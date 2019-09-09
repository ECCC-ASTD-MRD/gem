#!/usr/bin/perl
#
# s.dependencies.pl

#use warnings;        # Avertissement des messages d'erreurs
use strict;          # V�rification des d�clarations
use File::Spec::Functions;
use URI::file;
use Cwd "realpath";

#########################################################################
#                                                                       #
#                      List of global variables                         #
#                                                                       #
#########################################################################

my $msg =  2;
my %msgl = ("q", 0 , "e", 1 ,"w", 2 ,"v", 3 ,"vv", 4,"vvv", 5 ) ;
my $items_per_line = 4 ;   # number of items per Makefile line
my $item_count = 0;
my $ext = undef;
my @listfile;
my $use_strict = undef;
my $deep_include = undef;
my $soft_restriction = undef;
my $export_list = undef;
my @current_dependencies_list = ();
my @outside_tree_list = ();
my %module_missing = ();
my $bool = 0;
my %LISTOBJECT = ( ); # Hash of SRCFile object with filename, path and extension as the key

#########################################################################
#                                                                       #
#              List of function and object definition                   #
#                                                                       #
#########################################################################

{ package SRCFile;
    {
    
    # List of the type of file associated with the extension
    our %TYPE_LIST = (
        f     => "COMPILABLE",
        ftn   => "COMPILABLE",
        ptn   => "INCLUDE",
        f90   => "COMPILABLE",
        ftn90 => "COMPILABLE",
        ptn90 => "INCLUDE",
        cdk   => "INCLUDE",
        cdk90 => "COMPILABLE",
        c     => "COMPILABLE",
        h     => "INCLUDE",
        hf    => "INCLUDE",
        fh    => "INCLUDE",
	hc    => "INCLUDE",
	ch    => "INCLUDE",
        inc   => "INCLUDE",
	tmpl90 => "COMPILABLE",
    );

    # @new: Constructor of the SRCFile Object
    # 
    # input:
    #   $1 = {
    #           path => 
    #           filename =>
    #           extension =>
    #        }
    #
    # output:
    #   pointer to the object

    sub new 
    {
        my ( $class, $ref_arguments ) = @_;
        
        $class = ref($class) || $class;
        my $this = {};
        bless( $this, $class );
            
        $this->{FULLPATH_SRC}     = " ";
        $this->{FILENAME}         = " ";
        $this->{EXTENSION}        = " ";

        $this->{FULLPATH_SRC}     = $ref_arguments->{path}; 
        $this->{FILENAME}         = $ref_arguments->{filename};
        $this->{EXTENSION}        = $ref_arguments->{extension};
      %{$this->{DEPENDENCIES}}    = ();   # A list to the required file
        $this->{TYPE}             = $TYPE_LIST{lc $this->{EXTENSION}};   # Type based on the extension
        $this->{STATUS}           = undef;
      @{$this->{UNSOLVED_MODULE}} = ();
      @{$this->{MODULE_LIST}}     = ();
      %{$this->{UNKNOWN_MODULE}}  = ();
      %{$this->{UNKOWN_USE}}	  = ();
            
        return $this;
    }

    # @displ: Display file information
    #
    # input:
    #   none
    #
    # output:
    #   none

    sub displ
    {
        my $self = $_[0];
        my $key;
    
        print "\t$_[0]->{FULLPATH_SRC} \t$_[0]->{FILENAME} \t$_[0]->{EXTENSION} \t$_[0]->{TYPE}\t\t$_[0]->{STATUS}\n";
    
        foreach $key (keys %{$self->{DEPENDENCIES}})
        {
            print "\t\tDEPENDENCIES:\t $key\n";
        }
    
        foreach $key (@{$self->{MODULE_LIST}})
        {
            print "\t\tMODULE:\t $key\n";
        }
    
        foreach $key (@{$self->{UNSOLVED_MODULE}})
        {
            print "\t\tMISSING MODULE:\t $key\n";
        }
    } 

    # @getFilename: Get full filename with path and extension
    # 
    # input: 
    #   none
    #
    # output:
    #   full path

    sub getFilename
    {
        return "$_[0]->{FULLPATH_SRC}$_[0]->{FILENAME}.$_[0]->{EXTENSION}";
    }

    # @set_status: Set status of object to true
    #
    # input: none
    #
    # output: none

    sub set_status
    {
        $_[0]->{STATUS} = 1; #true
    }

    # @get_status: get status the object
    #
    # input: none
    #
    # output: 
    #   status

    sub get_status
    {
        return $_[0]->{STATUS};
    }

    # @reset_status: Set status of object to false
    #
    # input: none
    #
    # output: none

    sub reset_status
    {
        $_[0]->{STATUS} = undef; #false
    }

    # @as_module: find if the object defined the module 
    #
    # intput: 
    #   $1 = Module name to find
    #
    # output:
    #   true (1) if the module as been found, false (undef) otherwise

    sub as_module
    {
        my $module_name = $_[1];
        my $module;
    
        foreach $module (@{$_[0]->{MODULE_LIST}})
        {
            return 1 if $module eq $module_name;
        }
    
        return undef;
    }

    # @as_unsolved_module: find if the object has the module in his unsolved module list
    #
    # input: 
    #   $1 = Module name to find
    #
    # output:
    #   true (1) if the module as been found, false (undef) otherwise

    sub as_unsolved_module
    {
        my $module_name = $_[1];
        my $module = "";
        
	foreach $module (@{$_[0]->{UNSOLVED_MODULE}})
        {
            
            return 1 if ($module eq $module_name);
        }
    
        return undef;
    }

    # @remove_unsolved_module: delete the module in the unsolved module list
    #
    # input:
    #   $1 = Module name to delete
    #
    # output: none

    sub remove_unsolved_module
    {
        my $module_name = "$_[1]";
        my $idx = 0;
        my $module = "";

	#print "MODULE? $_[0]->{FULLPATH_SRC} $_[0]->{FILENAME}.$_[0]->{EXTENSION}\n" if $module_name eq "";

	foreach $module (@{$_[0]->{UNSOLVED_MODULE}})
        {
	 #   print "TST: $module vs $module_name\n";
            if( $module eq $module_name )
	    { 
		delete ${$_[0]->{UNSOLVED_MODULE}}[$idx]; 
	#	print "MODULE DELETED\n"; 
	    } 
            $idx++; 
        }

    }

    # @find_depedencies: find if the object has the filename in his depedencie list
    #
    # input:
    #   $1 = filename to search
    # 
    # output
    #   true (1) if the module as been found, false (undef) otherwise

    sub find_depedencies
    {
        my $search_depedencies = "$_[1]";
	#print "SEARCHING IN : $_[0]->{FULLPATH_SRC} $_[0]->{FILENAME} $_[0]->{EXTENSION}\n";        

        while( my($dep_filename, $dep_ptr) = each(%{$_[0]->{DEPENDENCIES}}) )
        {
            return 1 if( ($search_depedencies eq $dep_filename) and ($_[0]->{FULLPATH_SRC} ne $dep_ptr->{FULLPATH_SRC} ) );
	    return undef if( ($search_depedencies eq $dep_filename) and ($_[0]->{FULLPATH_SRC} eq $dep_ptr->{FULLPATH_SRC} ) );
	    return 1 if( $dep_ptr->SRCFile::find_depedencies($search_depedencies) );
        }
        
	return undef;
    }

    } #end package SRCFile
}

# @reset_all_file: Reset status of all object
#
# input:
#   %0 = Hash of object
#
# output: none

sub reset_all_file
{
  keys %{$_[0]}; #reset current location of the hash

  while(my($key, $file) = each(%{$_[0]}))
  {
    $file->reset_status();
  }
}

# @search_undone_file: Find the key of the object that hasn't been processed yet (based on its status)
#
# input: 
#   %0 = Hash of object
#
# output:
#   key of the object, undef otherwise.

sub search_undone_file
{
  keys %{$_[0]}; #reset current location of the hash

  while(my($key, $file) = each(%{$_[0]}))
  {
    return $key if !$file->get_status(); 
  }

  return undef;
}

# @search_module: Find the key of the object that own the module name
#
# input: 
#   %0 = Hash of objects
#   $1 = Module name
#
# output:
#   key of the object where the module has been found, undef otherwise.

sub search_module
{
    keys %{$_[0]}; #reset current location of the hash
    my $module_name = $_[1];

    while(my($key, $file) = each(%{$_[0]}))
    {
        return $key if( $file->as_module($module_name) ); 
    }
    
    return undef;
}

# @search_unsolved_module: Find the key of the first object that has the module as one of his unsolved module list
#
# input:
#   %0 = Hash of object
#   $1 = Module name to be search
#
# output:
#   key of the object, undef otherwise.

sub search_unsolved_module
{
  keys %{$_[0]}; #reset current location of the hash
  my $module_name = $_[1];

  while((my($key, $file)) = each(%{$_[0]}))
  {
    return $key if( $file->SRCFile::as_unsolved_module($module_name) ); 
  }

  return undef;
}

# @print_header: Print the first line of a dependency rule or variable list
#
# input:
#   $0 = First word(s) of line
#   $1 = Seperator
#   $2 = Word/item right after seperator, empty if no word is needed
#
# output: none

sub print_header {
  $item_count = 0 ;
  my($item1,$separ,$item2) = @_ ;
  print STDOUT "$item1$separ" ;
  print STDOUT "\t$item2" if ( "$item1" ne "$item2" && "$item2" ne "" ) ;
}

# @print_item: print each item of dependency rule or variable list (items_per_line items per line)
#
# input:
#   $0 = Item to print
#
# output: none

sub print_item {
  my($item1) = @_ ;
  if( "$item1" ne "" )
  {
    if ( $item_count == 0 ) { print STDOUT " \\\n\t" ; }
    print STDOUT "$item1  " ;
    $item_count = 0 if ( $item_count++ >= $items_per_line ) ;
  }
}

# @find_same_filename: Look if the filename is already used somewhere else.
#
# input: 
#   %0 = Hash of objects
#   $1 = filename to compare with.
#
# output:
#    key (filename) of the object if the file already exist in the list, false (undef) otherwise.

sub find_same_filename {
    keys %{$_[0]}; #reset current location of the hash
    my $cmp_file = ${$_[0]}{$_[1]};
    
    # Return if the soft_restriction is activated and if the file isn't compilable. 
    # soft_restriction = disable warning if 2 headers files have the same name
    return undef if ( $soft_restriction and ($cmp_file->{TYPE} ne "COMPILABLE") );

    while(my($key, $file) = each(%{$_[0]}))
    {
        return $key if( ("$file->{FILENAME}.$file->{EXTENSION}" eq "$cmp_file->{FILENAME}.$cmp_file->{EXTENSION}") and ($file->{FULLPATH_SRC} ne $cmp_file->{FULLPATH_SRC}) ); 
    }
    
    return undef;
}

# @find_same_output: Look if the filename is already used in the Object list.
#
# input: 
#   %0 = Hash of objects
#   $1 = Object to compare with.
#
# output:
#   key (filename) of the object if the file already exist in the list, false (undef) otherwise.

sub find_same_output {
    keys %{$_[0]}; #reset current location of the hash
    my $cmp_file = ${$_[0]}{$_[1]};

    return undef if $cmp_file->{TYPE} ne "COMPILABLE";
    
    while(my($key, $file) = each(%{$_[0]}))
    {
        return $key if( ($file->{FILENAME} eq $cmp_file->{FILENAME}) and ($file->{TYPE} eq "COMPILABLE") and ($key ne $_[1]) ); 
    }
    
    return undef;
}

# @find_string_in_array:
# 
# input: 
#   @0 = Array to search in.
#   $1 = String to search.
#
# output:
#   1 if the string as been found, undef otherwise.

sub find_string_in_array {
    my @myArray = @{$_[0]}; 
    my $string = $_[1];

    foreach my $string_tmp (@myArray)
    {
	return 1 if ($string eq $string_tmp);
    }

    return undef;
}

# @rec_print_dependencies: 
#
# input:
#   %0 = Hash of objects
#   $1 = Filename to print dependencies
#
# output:
#   none

sub rec_print_dependencies {
    my $file = ${$_[0]}{$_[1]};
    
    while( my($dep_filename, $dep_ptr) = each(%{$file->{DEPENDENCIES}}) )
    {
        my $tmp_filename = $dep_filename;
        $tmp_filename = "$dep_ptr->{FULLPATH_SRC}$dep_ptr->{FILENAME}.o" if ($dep_ptr->{TYPE} eq "COMPILABLE");
    
	next if( ($_[1] eq $dep_filename) or find_string_in_array(\@current_dependencies_list, $tmp_filename) );
	
        print_item($tmp_filename);
	push @current_dependencies_list, $tmp_filename;

        # Recursively call the function to print all depedencies
        rec_print_dependencies(\%{$_[0]}, $dep_filename) if( $dep_ptr->{TYPE} ne "COMPILABLE" );
    }
}

# @as_legal_extension: 
#
# input: 
#   $0 = Extension to search
#
# output:
#   1 if the extension is valid, undef otherwise.

sub as_legal_extension
{
    my $search_extension = lc  $_[0];
    
    foreach my $extension (keys(%SRCFile::TYPE_LIST))
    {
	return 1 if $extension eq $search_extension;
    }
    
    return undef;
}

#########################################################################
#                                                                       #
#                        Main program beginning                         #
#                                                                       #
#########################################################################

# Si pas d'arguments, on prend STDIN

if ( $#ARGV > -1 ){
 @listfile = (@ARGV) ;
 $listfile[$#listfile] =~ s/\n/ /;
}

#
# Process command line arguments
#

foreach my $target (@listfile){
#  process options
    if( $target =~ /^-[-]*h.*$/ ) { print STDERR "usage: s.dependencies.pl -q|e|w|v|vv|vvv -strict|deep-include|soft-restriction -EXPoutput_of_produced_file -Ooutfile list_of_files\n" ; exit ; }
    if( $target =~ /^-(q|e|w|v|vv|vvv)$/ ) { $msg = $msgl{"$1"} ; next ; }
    if( $target =~ /^-strict$/ ) { $use_strict = 1; next; }
    if( $target =~ /^-deep-include$/ ) { $deep_include = 1; next; }
    if( $target =~ /^-soft-restriction$/ ) { $soft_restriction = 1; next; }
    
    #if( $target =~ /^-lib$/ ) { $prefix1="\$(EC_ARCH)/lib\$(MALIB)(" ; $prefix2=")" ; next ; }
    #if( $target =~ /^-dir$/ ) { $prefix1="\$(EC_ARCH)/" ; $prefix2="" ; next ; }
    #if( $target =~ /^-obj$/ ) { $prefix1="" ; $prefix2="" ; next ; }
    if( $target =~ /^-EXP([a-zA-Z0-9.\/]*)$/ ) { $export_list = $1; next ; }
    if( $target =~ /^-O([a-zA-Z0-9.\/]*)$/ ) { open(STDOUT,">", "$1") or die "ERROR: Can't redirect STDOUT\n" ; next ; }
#    process element of list_of_files_and_directories
    print STDERR "target = '$target'\n" if( $msg >= 5 );
}


while(my $target = <STDIN>){
    foreach my $entry (glob $target) {

	if( ! -f "$entry" ) { print STDERR "INFO: entry '$entry' is not a file\n" if( $msg >= 5 ); next; }
        my $file = "$entry" ;
        $file =~ s/,v$// ;
        $file =~ s/[\s]+// ;
        $file = File::Spec->abs2rel( realpath($file), "./");
        
        next if $file !~  /(.*\/)*(.*)[.]([^.]*$)/ ;  # ,v and path trimmed filename must be optional/path/root_file.extension

        next if exists $LISTOBJECT{$file}; 

	#print $1;
        my $path = ($1 ? $1 : "");
        my $filn = ($2 ? $2 : "");
        my $exte = ($3 ? $3 : "");
        my $duplicated_filename = "";

	next if !as_legal_extension($exte);

	$LISTOBJECT{"$path$filn.$exte"} = new SRCFile({path => $path, filename => $filn, extension => $exte});
    
	# Error handler
        die "ERR: using 2 files with the same name $duplicated_filename with $path$filn.$exte" if( $duplicated_filename = find_same_filename(\%LISTOBJECT, "$path$filn.$exte"));
        die "ERR: using 2 files ($duplicated_filename and $path$filn.$exte) that will produce the same object file ($filn.o)\n" if ( $duplicated_filename = find_same_output(\%LISTOBJECT, "$path$filn.$exte") );

    }  # foreach $entry

}  #  foreach $target

reset_all_file(\%LISTOBJECT);

while( my $filename = search_undone_file(\%LISTOBJECT) )
{
    open(INPUT,"<", $filename) or print STDERR "ERROR: Can't open file '".$filename."\n"; #  if( $msg >= 1 )
    my $file = $LISTOBJECT{$filename};
    my $line_number = 0;

    while (<INPUT>) 
    {
        if ($_ =~ /^[@]*[\s]*#[\s]*include[\s]*[<'"\s]([\w.\/\.]+)[>"'\s][\s]*/)
        {
            my $include_path = "";
            my $tmp_dir = $1;
        
	    if ($tmp_dir =~ /^\.\.\//)
            {
                $include_path = File::Spec->abs2rel( realpath("$file->{FULLPATH_SRC}/$tmp_dir"), "./"); # Convert file path position relatively to the base path
            }
            elsif (-f realpath("$file->{FULLPATH_SRC}/$tmp_dir") )
            {
                $include_path = File::Spec->abs2rel( realpath("$file->{FULLPATH_SRC}/$tmp_dir"), "./");
            }
            else
            {
                $include_path = File::Spec->abs2rel( realpath($tmp_dir), "./");
            }

	    next if $include_path !~  /(.*\/)*(.*)[.]([^.]*$)/ ;  # ,v and path trimmed filename must be optional/path/root_file.extension

	    my $path = ($1 ? $1 : "");
            my $filn = ($2 ? $2 : "");
            my $exte = ($3 ? $3 : "");
            my $duplicated_filename = "";

	    next if !as_legal_extension($exte);

	    if( ! -f "$path$filn.$exte")
	    {
		push @outside_tree_list, "$path$filn.$exte" if(!find_string_in_array(\@outside_tree_list, "$path$filn.$exte"));
		next;
	    }

            # Add file in the database if it's not in yet and if the file really exists.
            $LISTOBJECT{"$path$filn.$exte"} = new SRCFile({path => $path, filename => $filn, extension => $exte}) if( !exists $LISTOBJECT{"$path$filn.$exte"});

            # Force the file to not be analysed.
            $LISTOBJECT{"$path$filn.$exte"}->set_status() if $deep_include;

	    # Error handler
            die "ERR: using 2 files with the same name $duplicated_filename with $path$filn.$exte\n" if( $duplicated_filename = find_same_filename(\%LISTOBJECT, "$path$filn.$exte"));
            die "ERR: using 2 files ($duplicated_filename and $path$filn.$exte) that will produce the same object file ($filn.o)\n" if ( $duplicated_filename = find_same_output(\%LISTOBJECT, "$path$filn.$exte") );
            die "ERR: cannot include compilable file ($tmp_dir) in $filename while using strict mode\n" if ( $use_strict and $LISTOBJECT{"$path$filn.$exte"}->{TYPE} eq "COMPILABLE" );

            # Add to dependencies, if not already there
            ${$file->{ DEPENDENCIES }}{"$path$filn.$exte"} = $LISTOBJECT{"$path$filn.$exte"} if( !exists ${$file->{ DEPENDENCIES }}{"$path$filn.$exte"} );

        }
        next if ( $file->{EXTENSION} =~ /(c|cc|CC)$/);
        
        # FORTRAN include statement : include "..."    include ',,,"
        if ($_ =~ /^[@]*[\s]*include[\s]*[<'"\s]([\w.\/\.]+)[>"'\s][\s]*/i)
        {
            my $include_path = "";
            my $tmp_dir = $1;
    
            if ($tmp_dir =~ /^\.\.\//)
            {
                $include_path = File::Spec->abs2rel( realpath("$file->{FULLPATH_SRC}/$tmp_dir"), "./"); # Convert file path position relatively to the base path
            }
            elsif (-f realpath("$file->{FULLPATH_SRC}/$tmp_dir"))
            {
                $include_path = File::Spec->abs2rel( realpath("$file->{FULLPATH_SRC}/$tmp_dir"), "./");
            }
            else 
            {
                $include_path = File::Spec->abs2rel( realpath($tmp_dir), "./");
            }

            next if $include_path !~  /(.*\/)*(.*)[.]([^.]*$)/ ;  # ,v and path trimmed filename must be optional/path/root_file.extension

            my $path = ($1 ? $1 : "");
            my $filn = ($2 ? $2 : "");
            my $exte = ($3 ? $3 : "");
            my $duplicated_filename = "";

	    next if !as_legal_extension($exte);

	    if( ! -f "$path$filn.$exte")
	    {
		push @outside_tree_list, "$path$filn.$exte" if(!find_string_in_array(\@outside_tree_list, "$path$filn.$exte"));
		next;
	    }

            # Add file in the database if it's not in yet and if the file really exists.
            $LISTOBJECT{"$path$filn.$exte"} = new SRCFile({path => $path, filename => $filn, extension => $exte}) if( !exists $LISTOBJECT{"$path$filn.$exte"});

            # Force the file to not be analysed.
            $LISTOBJECT{"$path$filn.$exte"}->set_status() if $deep_include;

            # Error handler
            die "ERR: using 2 files with the same name $duplicated_filename with $path$filn.$exte" if( $duplicated_filename = find_same_filename(\%LISTOBJECT, "$path$filn.$exte"));
            die "ERR: using 2 files ($duplicated_filename and $path$filn.$exte) that will produce the same object file ($filn.o)\n" if ( $duplicated_filename = find_same_output(\%LISTOBJECT, "$path$filn.$exte") );
            die "ERR: cannot include compilable file ($tmp_dir) in $filename while using strict mode" if ( $use_strict and $LISTOBJECT{"$path$filn.$exte"}->{TYPE} eq "COMPILABLE" );

            # Add to dependencies, if not already there
            ${$file->{ DEPENDENCIES }}{"$path$filn.$exte"} = $LISTOBJECT{"$path$filn.$exte"} if( !exists ${$file->{ DEPENDENCIES }}{"$path$filn.$exte"} );
        }
        # FORTRAN use statement : use yyy 
        if ($_ =~ /^[@]*[\s]*\buse[\s]+([a-z][\w]*)(,|\t| |$)/i)
        {
            my $modname = $1 ; $modname =~ tr/A-Z/a-z/ ; # modules names are case insensitive

            if( my $include_filename = search_module(\%LISTOBJECT, $modname)) # If the module can be found, add the file to dependencies
            {
              ${$file->{ DEPENDENCIES }}{$include_filename} = $LISTOBJECT{$include_filename} if( !exists ${$file->{ DEPENDENCIES }}{$include_filename} );
            }
            else #module not found yet!
            {
              push @{$file->{UNSOLVED_MODULE}}, $modname; 
            }

        }

	elsif($_ =~ /^[@]*[\s]*\buse[\s]+([a-z][\w]*)/i)
	{
	    ${$file->{ UNKOWN_USE }}{$line_number} = $_;
	}

        # FORTRAN module declaration : module xxx
        if ($_ =~ /^[@]*[\s]*\bmodule[\s]+([a-z][\w]*)(,|\t| |$)/i)
        {
            my $modname = $1 ; $modname =~ tr/A-Z/a-z/ ; # modules names are case insensitive
            my $search_filename = '';

	    next if $modname eq "procedure";

	    #print "SEARCHING MODULE: ";
	    #print "$modname\n";

            # Verifier que le nom du module n'existe pas dans un autre fichier
            if($search_filename = search_module(\%LISTOBJECT, $modname))
            { 
                print STDERR "Module ".$modname." (".$filename.") already defined in ".$search_filename."\n"; 
                next; 
            }

            # Ajouter le module dans la liste des modules associer au fichier.
            push @{$file->{ MODULE_LIST }}, $modname;

            # Recherche tous les fichiers analyser precedemment qui avait besoin de ce module la
            while( my $key = search_unsolved_module(\%LISTOBJECT, $modname) )
            {
                #print "unsolved module: $key".${$LISTOBJECT{$key}->{ DEPENDENCIES }}{$filename}."\n";
                # Ajouter a la liste des dependence, le fichier en cours
                ${$LISTOBJECT{$key}->{ DEPENDENCIES }}{$filename} = $file if( !exists ${$LISTOBJECT{$key}->{ DEPENDENCIES }}{$filename} );

                # Enlever le module de la liste des unsolved modules 
                $LISTOBJECT{$key}->remove_unsolved_module($modname);
            }
        }
	elsif ($_ =~ /^[@]*[\s]*\bmodule[\s]+/i)
	{
	    ${$file->{ UNKOWN_MODULE }}{$line_number} = $_;
	    #print STDERR "Unknown module statement: $filename: $_\n";
	}
	$line_number++;
    }
    
    $file->set_status();

    close INPUT;

}

while( my($filename, $file) = each(%LISTOBJECT) )
{
    my $result;
    #print "TEST: $filename \n";

    # Loop sur tous les fichiers pour voir les dependances circulaires
    if($result = $file->find_depedencies($filename))
    { 
        print STDERR "ERR: Circular dependencies in $filename FAILED\n";
        exit 1; 
    }

    #print "($result)\n";
}

#print "DONE\n";

#
#  lists of file types FDECKS, CDECKS, ...
#

reset_all_file(\%LISTOBJECT);

for $ext (keys %SRCFile::TYPE_LIST)
{
    print_header(uc $ext."DECKS", "=", "");

    foreach my $filename (sort keys %LISTOBJECT)
    {
	my $file = $LISTOBJECT{$filename};
        print_item($filename) if( $file->{EXTENSION} eq $ext );
    }

    print STDOUT "\n";
}

#
#  OBJECTS LIST
#

print_header("OBJECTS","=","");

foreach my $filename (sort keys %LISTOBJECT)
{
    my $file = $LISTOBJECT{$filename};
    print_item("$file->{FULLPATH_SRC}$file->{FILENAME}.o") if( $file->{TYPE} eq "COMPILABLE" );
}

print STDOUT "\n";

#
#   Build depedencie rules
#

foreach my $filename (sort keys %LISTOBJECT)
{
    my $file = $LISTOBJECT{$filename};
    @current_dependencies_list = ();

    if( $file->{TYPE} eq "COMPILABLE" )
    {
        print_header("$file->{FULLPATH_SRC}$file->{FILENAME}.o",":","$filename");

        rec_print_dependencies(\%LISTOBJECT, $filename);

        print "\n";
    }
}

#
#    Print the missing module(s) / file(s) from the current tree
#

print STDERR "Includes missing from the current tree: " if( $#outside_tree_list );

foreach my $filename (@outside_tree_list)
{
    print STDERR "$filename ";
}
print STDERR "\n" if( $#outside_tree_list );

%module_missing = ();
$bool = 0;
while( my($filename, $file) = each(%LISTOBJECT) )
{
    foreach my $module (@{$file->{UNSOLVED_MODULE}})
    {
	next if ($module eq "");

	if (!$bool)
	{
	    $bool = 1;
	    print STDERR "Modules missing from the current tree: ";
	}
	$module_missing{$module} = $filename if( !exists $module_missing{$module});
    }
}

while( my($module,$filename) = each(%module_missing))
{
    print STDERR "$module ($filename) ";
}

print STDERR "\n" if $bool;

#
#   Unknown module statement 
#

$bool = 0;

while( my($filename, $file) = each(%LISTOBJECT) )
{
    while( my($line_number,$text_line) = each(%{$file->{UNKNOWN_MODULE}}))
    {
	if (!$bool)
	{
	    $bool = 1;
	    print STDERR "Unknown module statement: \n";
	}

	print STDERR "\t($filename) $line_number: $text_line\n";
    }
}

#
#   Unknown use statement 
#

$bool = 0;

while( my($filename, $file) = each(%LISTOBJECT) )
{
    while( my($line_number,$text_line) = each(%{$file->{UNKOWN_USE}}))
    {
	if (!$bool)
	{
	    $bool = 1;
	    print STDERR "Unknown use statement: \n";
	}

	print STDERR "\t($filename) $line_number: $text_line\n";
    }
}

#
#   Export a list of produced files (.o and .mod)
#

if($export_list)
{
    open(my $EXPOUT,'>',$export_list);

    my @list_of_modules = ();

    foreach my $filename (sort keys %LISTOBJECT)
    {
	my $file = $LISTOBJECT{$filename};
	print $EXPOUT "$file->{FULLPATH_SRC}$file->{FILENAME}.o\n" if $file->{TYPE} eq "COMPILABLE";

	foreach my $module (@{$file->{MODULE_LIST}})
	{
	    push @list_of_modules, $module if $module ne "";
	}
    }

    foreach my $module (sort @list_of_modules)
    {
	print $EXPOUT "$module.mod\n";
    }

    close($EXPOUT);
}

## output list of produced files (.mod .o) in specified file. Else do nothing! 

# DEBUG tree
# while( my($filename, $file) = each(%LISTOBJECT) )
# {
#     $file->displ();
# }

