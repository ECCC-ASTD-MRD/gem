/* RMNLIB - Library of useful routines for C and FORTRAN programming 
 * Copyright (C) 1975-2005  Environnement Canada
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation,
 * version 2.1 of the License.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */
#include <stdlib.h>
#include <unistd.h>
#include <malloc.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include <rpnmacros.h>

#define C_SOURCE_CODE
#include <WhiteBoard.h>

#define TRIM(string,Lstring) { while( Lstring>0 && string[Lstring-1]==' ' ) Lstring-- ; }

/* default verbosity level, print warnings and up */
static int message_level = WB_MSG_WARN;

/* set verbosity level (C callable) */
int c_wb_verbosity(int level)
{
   int old_level = message_level;
   message_level = level;
   return(old_level);
}

/* set verbosity level (FORTRAN callable) */
wordint f77_name(f_wb_verbosity)(wordint *level)
{
   wordint old_level = message_level;
   message_level = *level;
   return(old_level);
}

static char *linebuffer=NULL;              /* line buffer used to process directive files */
static char *current_char=NULL;            /* pointer to current character in directive file */

static DEFTBL *definition_table=NULL;      /* pointer to definition table */
static int definition_table_entries=0;     /* number of valid entries in definition table */
static int max_definition_table_entries=0;

/* manage definition_table (used by whiteboard read routine) */
/* define_mode=0 : processing a define */
/* define_mode=1 : processing key=value in strict dictionary mode */
/* define_mode=2 : processing key=value in normal mode */
static int wb_define_check(LINE *line, int define_mode){
   int i;

   for( i=0 ; i<definition_table_entries ; i++) {
      if( definition_table[i].line == line ) {       /* a match has been found */

         if( define_mode==1 ) return(WB_ERR_REDEFINE) ; /* found when it should not (we are processing a define) */

         if( definition_table[i].assigned ) {  /* already assigned to, should not be assigned again */
            if( message_level<=WB_MSG_ERROR )
               fprintf(stderr,"ERROR: key value already assigned in this directive file \n");
            return(WB_ERR_REDEFINE) ;
            }
         definition_table[i].assigned = 1 ;    /* flag assignation */

         return(i);
         }
   }
   if( define_mode == 0 ) return (WB_ERR_NOTFOUND);  /* not in define mode, strict dictionary mode, not found, fail */

   if(definition_table_entries >= max_definition_table_entries) {
      if( message_level<=WB_MSG_ERROR )
         fprintf(stderr,"ERROR: too many keys appear in this directive file \n");
      return (WB_ERROR);  /* OOPS, definition table is full */
      }

   definition_table[definition_table_entries].line=line;   /* create and initialize entry */
   definition_table[definition_table_entries].defined = ( define_mode==1 )? 1 : 0 ;   /* define(...) */
   definition_table[definition_table_entries].assigned= ( define_mode!=1 )? 1 : 0 ;   /* key=...     */
   return(definition_table_entries++); /* return position*/
}

static int default_dict_option = 0;              /* default options when creating variables (used by Whiteboard read) */

static SYMTAB dict_options[] = {                 /* table of variable creation options (used by Whiteboard read) */
                { WB_IS_LOCAL,"WB_IS_LOCAL"},
                { WB_REWRITE_AT_RESTART,"WB_REWRITE_AT_RESTART"},
                { WB_REWRITE_NONE,"WB_REWRITE_NONE"},
                { WB_REWRITE_UNTIL,"WB_REWRITE_UNTIL"},
                { WB_REWRITE_MANY,"WB_REWRITE_MANY"},
                { WB_READ_ONLY_ON_RESTART,"WB_READ_ONLY_ON_RESTART"},
                { 0, NULL}        } ;

static SYMTAB verb_options[] = {                  /* table for verbosity options (used by Whiteboard read) */
                { WB_MSG_DEBUG,"WB_MSG_DEBUG"},
                { WB_MSG_INFO,"WB_MSG_INFO"},
                { WB_MSG_WARN,"WB_MSG_WARN"},
                { WB_MSG_ERROR,"WB_MSG_ERROR"},
                { WB_MSG_SEVERE,"WB_MSG_SEVERE"},
                { WB_MSG_FATAL,"WB_MSG_FATAL"},
                { 0, NULL}        } ;

typedef PAGE *PAGEPTR;
typedef struct{                      /* a whiteboard instance */
   PAGE *PageChain;
   int validpages;
   }WhiteBoard;
typedef WhiteBoard *WhiteBoardPtr;

static WhiteBoard DummyWhiteboard;    /* dummy whiteboard instance, returned when freeing a whiteboard */
static WhiteBoardPtr DummyWhiteboardPtr=&DummyWhiteboard;
static WhiteBoard BaseWhiteboard;     /* permanent whiteboard instance, the only one that can be checkpointed */
static WhiteBoardPtr BaseWhiteboardPtr=&BaseWhiteboard;

static int new_page(WhiteBoard *WB, int nlines);

/* create a new whiteboard instance */
 WhiteBoard *c_wb_new(){
   int status;
   WhiteBoard *New_WhiteBoard=(WhiteBoard *)calloc(1,sizeof(WhiteBoard));

   if( New_WhiteBoard == NULL ) return(NULL);   /* OOPS, too many whiteboard instances  */
   New_WhiteBoard->PageChain=NULL;
   New_WhiteBoard->validpages=0;
   status=new_page(New_WhiteBoard,WB_MAXLINESPERPAGE);
   if(status<0) return(NULL);
   return(New_WhiteBoard);
}
wordint f77_name(f_wb_new)(WhiteBoard **WB){
   WhiteBoard *New_WhiteBoard=c_wb_new();
   *WB=New_WhiteBoard;
   return( New_WhiteBoard==NULL ? WB_ERROR : WB_OK);
}

/* get rid of a whiteboard instance */
int c_wb_free(WhiteBoard *WB) {
   PAGE *temp, *temp2;

   if( WB==NULL ) return(WB_OK);      /* global whiteboard page chain is never freed */

   if( WB-> validpages <= 0 ) return(WB_ERROR) ; /* invalid whiteboard, negative or zero number of valid pages */
   WB->validpages = 0;
   temp=WB->PageChain;
   while( temp != NULL ) {    /* free all pages for this whiteboard */
      temp2=temp;             /* keep copy of pointer to current page */
      temp = temp->next;      /* pointer to next page in page chain */
      free(temp2);            /* free current page */
      }
   free(WB);
   return(WB_OK);
}
wordint f77_name(f_wb_free)(WhiteBoard **WB){
   wordint status=c_wb_free(*WB);
   *WB = DummyWhiteboardPtr ;  /* set whiteboard address to address of dummy whiteboard, so that we can trap it */
   return status;
}

static SYMTAB errors[] = {      /* error message table */
                { WB_ERR_NAMETOOLONG , "key name is too long" },
                { WB_ERR_NOTFOUND, "requested item not found" },
                { WB_ERR_READONLY, "key value cannot be redefined"},
                { WB_ERR_WRONGTYPE, "key object and value have different types/length combinations" },
                { WB_ERR_WRONGDIMENSION, "key object has dimensions smaller than value assigned" },
                { WB_ERR_WRONGSTRING, "target string element is too short" },
                { WB_ERR_ALLOC, "could not allocate memory with malloc" },
                { WB_ERR_NOTYPE, "requested type/length combination is not valid" },
                { WB_ERR_NOMEM, "not enough space available to allocate value" },
                { WB_ERR_NOVAL, "key has not been initialized yet" },
                { WB_ERR_BADVAL, "key value is dubious due to a failed put" },
                { WB_ERR_CKPT, "problem opening/reading/writing checkpoint file" },
                { WB_ERR_REDEFINE, "attempt to redefine an entry" },
                { WB_ERR_BIG, "token too large in directive " },
                { WB_ERR_SYNTAX, "syntax error in directive " },
                { WB_ERR_OPTION, "invalid option or combination of options" },
                { WB_ERR_READ, "open/read error in file" },
                { WB_ERROR, "" },
                { 0, NULL}        } ;

static char *datatypes[] = { "    ", "Real", "Int ", "Char", "Bool"  }; /* THIS MUST REMAIN CONSISTENT WITH Whiteboard.h */

/* get integer value associated with a string(symbol) from a symbol table */
static int wb_value(char *text, SYMTAB *table){

   if( message_level<=WB_MSG_DEBUG )
      fprintf(stderr,"looking for value of '%s', ",text);
   while( table->text != NULL ){
      if(strcmp(table->text,text)==0) {
         if( message_level<=WB_MSG_DEBUG ) fprintf(stderr,"found %d\n",table->code);
         return(table->code);
         }
      table++;
      }
   if( message_level<=WB_MSG_DEBUG ) fprintf(stderr,"found NONE\n");
   return(0);  /* no value found in table */
}

/* default user error handler, do nothing routine */
static void null_error_handler()
{
return ;
}

/* pointer to the user's FORTRAN handler routine, defaults to internal null_error_handler */
/* the error handler receives two pointers to wordint variables */
static void (*ErrorHandler)() = &null_error_handler;

/* FORTRAN callable routine to define user error handler */
void f77_name(f_wb_error_handler)(  void (*UserErrorHandler)() )
{
   ErrorHandler = UserErrorHandler;
}

static char text_of_last_error[256];
static char *extra_error_message=NULL;
static char extra_error_buffer[256];
/* copy a possibly non null terminated string to the extra error buffer */
static void set_extra_error_message(char *name,int Lname){
   if(Lname > 255) Lname = 255;
   extra_error_buffer[Lname]='\0';
   while(Lname > 0) { extra_error_buffer[Lname-1]=name[Lname-1]; Lname-- ; }
   extra_error_message=&(extra_error_buffer[0]);
}
/*
 error processing routine, returns the error code it was given
 a message will be printed on standard error if severity is >= message level
 the user error handler will be called if it was defined
*/
static int wb_error(int severity, int code){
   SYMTAB *message=&errors[0];
   wordint Severity=severity;
   wordint Code=code;

   if( code == WB_OK || severity<message_level ) return(code);
   while( message->text != NULL ){
      if( message->code == code ){
         snprintf(text_of_last_error,sizeof(text_of_last_error)-1,"%s - %s",
                  message->text,(extra_error_message!=NULL)?extra_error_message:"");
         fprintf(stderr,"ERROR: %s \n", text_of_last_error);
         extra_error_message=NULL;
         break;
         }
      message++;
      }
   (*ErrorHandler)(&Severity,&Code);  /* call user's FORTRAN error handler */
   return(code);
}

/* error exit macro */
#define WB_ERR_EXIT(level,code)  return( wb_error(level,code) )

/*
  get the symbolic type code (WB_FORTRAN_...) , validate the type/length combination,
  return the symbolic type code if OK  (WB_FORTRAN_...)
  WB_FORTRAN_REAL must be 4 or 8 bytes long
  WB_FORTRAN_INT  must be 4 or 8 bytes long
  WB_FORTRAN_BOOL must be 4 bytes long
  WB_FORTRAN_CHAR must be 1 to WB_MAXSTRINGLENGTH bytes long
  Ltype_char is a way to permit a possibly excessive length for a character string on a get operation
*/
static int get_typecode(unsigned char typecode, int Ltype, int Ltype_char){
   if( typecode == WB_FORTRAN_REAL && (Ltype==4 || Ltype==8) ) return(WB_FORTRAN_REAL);
   if( typecode == WB_FORTRAN_INT && (Ltype==4 || Ltype==8) ) return(WB_FORTRAN_INT);
   if( typecode == WB_FORTRAN_BOOL && Ltype==4 ) return(WB_FORTRAN_BOOL);
   if( typecode == WB_FORTRAN_CHAR && (Ltype>0 && Ltype_char<=WB_MAXSTRINGLENGTH) ) return(WB_FORTRAN_CHAR);
   WB_ERR_EXIT(WB_MSG_ERROR,WB_ERR_NOTYPE);   /* type / length combination is not valid */
}

/*
  FORTRAN/C style trimmed source to blank padded destination copy, return trimmed source length
  the src string can be either null byte (C style) terminated or up to lsrc bytes long (FORTRAN style)
  trailing blanks are ignored in both cases
  if destination pointer is NULL, no copy takes place, just return trimmed source length 
*/
static int c_fortran_string_copy(unsigned char *src, unsigned char *dest, int lsrc, int ldest, unsigned char pad)
{
   int i;

   /* if there is a null before lsrc characters in the source string, act as if it was terminated at null */
   for (i=0 ; src[i]!=0 && i<lsrc ; i++);
   lsrc = i;
   while( lsrc > 0 && src[lsrc-1] == ' ' ) lsrc--;  /* ignore trailing blanks in source string */
   if( dest == NULL ) return(lsrc);     /* no destination, just return trimmed source length */
   if(lsrc > ldest) return(-1);  /* OOPS, source longer than destination */
   if(pad==0 && lsrc==ldest) return(-1);  /* OOPS, not enough space for padding */
   for( i=0 ; i<lsrc ; i++ ) dest[i] = src[i] ;  /* copy src to dest */
   if(pad)
      while( i < ldest ) dest[i++] = pad ;          /* pad destination */
   else
      dest[i] = pad;
   return(lsrc) ; /* return number of significant characters copied */
}

#ifdef NOT_USED
/* FORTRAN callable version of above routine */
wordint f77_name(fortran_string_copy)(unsigned char *src, unsigned char *dest, F2Cl lsrc, F2Cl ldest)
{
   int Ldest = ldest;
   int Lsrc = lsrc;
   return(c_fortran_string_copy(src, dest, Lsrc, Ldest, ' '));
}
#endif

/* initialize a line with a name coming from a WB_FORTRAN string */
static void fill_line(LINE *line, unsigned char *name, int lname)  
{
   int i;
   unsigned char c;

   for( i=0 ; i < WB_MAXNAMELENGTH ; i++){
      c = (i<lname) ? name[i] : ' ';   /* blank fill to WB_MAXNAMELENGTH (FORTRAN style) */
      line->m.name.c[i] = toupper(c) ; /* force uppercase */
      }
   line->m.name.c[WB_MAXNAMELENGTH] = 0;  /* force NULL terminator */
}

/* short match for names, stop at first NULL encountered, match at most nc characters */
static int wb_match_name(unsigned char *c1,unsigned char *c2, int nc){
   int j;
   if(nc == 0) return(1);  /* everything matches 0 characters */
   for( j=0 ; j<nc && *c2!=0 && *c1!=0 ; j++ ) { if( toupper(*c1) != toupper(*c2) ) break ; c1++ ; c2++ ; }
   if( (j==nc)  || (*c1==0) || (*c2==0) ) return(1);
/* fprintf(stderr,"c1='%c',c2='%c'\n",*c1,*c2); */
/*   if( (toupper(*c1)==toupper(*c2)) || (*c1==0) || (*c2==0) || (j==nc) ) return(1); */
   return(0);
}

/*
 find a line in a page using the key name, return line number if found, error if not
 if err_not_found is false, it is a silent error if line is not found
*/
static int lookup_line(LINE *line, PAGE *page, int err_not_found) 
{
   int i=0;
   int linelimit;
   if(page != NULL) {
      linelimit=page->NextFreeEntry - 1;
      while( i <= linelimit ){
         LINE *pageline=&(page->line[i]);
         unsigned char *c1=&(pageline->m.name.c[0]), *c2=&(line->m.name.c[0]);
         if( wb_match_name(c1,c2,WB_MAXNAMELENGTH) ) return(i);      /* target found, return it's position in page */

         if(pageline->m.flags.lines == 0) break;   /* only one entry in this page in this case */
         i += pageline->m.flags.lines;
         }
      }
    if(err_not_found) WB_ERR_EXIT(WB_MSG_ERROR,WB_ERR_NOTFOUND);
    return(WB_ERR_NOTFOUND);
}

/* allocate new whiteboard page containing nlines lines */
static int new_page(WhiteBoard *WB, int nlines){
   PAGE *page;
   int pagesize;

   if(WB == NULL) WB=BaseWhiteboardPtr;
   if( message_level<=WB_MSG_DEBUG )
      fprintf(stderr,"allocating new page for %d lines \n",nlines);
   /* note: sizeof(PAGE) gives the size of a page containing WB_MAXLINESPERPAGE lines */
   /* adjustment for actual number of lines in page has to be done to allocate correct size */
   pagesize = sizeof(PAGE) + sizeof(LINE)*(nlines-WB_MAXLINESPERPAGE);
   page=(PAGE *)malloc( pagesize );
   if(page==NULL) WB_ERR_EXIT(WB_MSG_FATAL,WB_ERR_ALLOC);        /* OOPS, not enough memory, malloc failed */
   memset(page,0,pagesize);  /* zero page contents */

   page->NextFreeEntry = 0;       /* next free entry is first entry in page */
   page->LinesInPage = nlines;    /* the page can contain nlines data/metadata lines */
   page->next = WB->PageChain;        /* put current page chain after this page */
   WB->PageChain = page;              /* put new page at beginning of page chain*/
   WB->validpages++;
   return (0);
}

static char *WhiteBoardCheckpointFile="Whiteboard.ckpt";

int c_wb_reload();
/* whiteboard initialization routine, will try to reload the checkpoint file if one is found */
static int wb_init(){
   int status;

   if( BaseWhiteboardPtr->PageChain == NULL ) {   
      int fd;
      fd=open(WhiteBoardCheckpointFile,O_RDONLY);
      if(fd >= 0) {
         close(fd);
         if( message_level<=WB_MSG_INFO )
            fprintf(stderr,"whiteboard checkpoint file found, loading it\n");
         return( c_wb_reload() );
         }
      status=new_page(BaseWhiteboardPtr, WB_MAXLINESPERPAGE ); /* first time through, allocate first page */
      if ( status < 0 ) return(status);        /* allocation failed, OUCH !! */
      }
   return(0);
}
int c_wb_checkpoint_name(char *filename){
   extra_error_message="Setting checkpoint file name";
   WhiteBoardCheckpointFile=(char *)malloc(strlen(filename));
   if(WhiteBoardCheckpointFile==NULL) WB_ERR_EXIT(WB_MSG_FATAL,WB_ERR_ALLOC);
   strncpy(WhiteBoardCheckpointFile,filename,strlen(filename));
   return(wb_init());
}
wordint f77_name(f_wb_checkpoint_name)(char *filename, F2Cl lfilename){
   int Lfilename=lfilename;
   extra_error_message="Setting checkpoint file name";
   WhiteBoardCheckpointFile=(char *)malloc(Lfilename+1);
   if(WhiteBoardCheckpointFile==NULL) WB_ERR_EXIT(WB_MSG_FATAL,WB_ERR_ALLOC);
   c_fortran_string_copy(filename,WhiteBoardCheckpointFile,Lfilename,Lfilename,'\0');
   return(wb_init());
}
int c_wb_checkpoint_get_name(char *filename,int Lfilename){
   extra_error_message="NO checkpoint file found while getting checkpoint file name";
   if(WhiteBoardCheckpointFile==NULL) WB_ERR_EXIT(WB_MSG_FATAL,WB_ERR_ALLOC);
   strncpy(filename,WhiteBoardCheckpointFile,Lfilename);
   return(wb_init());
}

wordint f77_name(f_wb_checkpoint_get_name)(char *filename, F2Cl lfilename){
   int Lfilename=lfilename;
   extra_error_message="NO checkpoint file found while getting checkpoint file name";
   if(WhiteBoardCheckpointFile==NULL) WB_ERR_EXIT(WB_MSG_FATAL,WB_ERR_ALLOC);
   c_fortran_string_copy(WhiteBoardCheckpointFile,filename,Lfilename,Lfilename,' ');
   return(wb_init());
}

/* return pointer to first line of nlines consecutives lines in a whiteboard page */
/* as well as pointer to page containing these lines */
/* whiteboard may get extended by one page as a result of this call */
static LINE *new_line(WhiteBoard *WB, PAGE **page, int nlines){
   PAGE *lookuppage;
   LINE *result;
   int status;

   if( BaseWhiteboardPtr->PageChain == NULL ) { /* first time through */
      status=wb_init();
      if ( status < 0 ) return(NULL);           /* init failed, OUCH !! */
      }
   if(WB == NULL) WB=BaseWhiteboardPtr;
   lookuppage = WB->PageChain; /* start at beginning of pagechain */
   while( lookuppage->NextFreeEntry+nlines > lookuppage->LinesInPage ) {  /* while page is not full */
      if( lookuppage->next == NULL ){                               /* already last page, allocate a new page */
         status=new_page(WB, ( nlines > WB_MAXLINESPERPAGE ) ? nlines : WB_MAXLINESPERPAGE ); /* make sure allocated page is big enough */
         if ( status < 0 ) return(NULL);                            /* allocation of new page failed, OUCH !! */
         lookuppage=WB->PageChain;                                  /* potential target is new allocated page */
         break;
         }
      else
         lookuppage = lookuppage->next;                             /* look into next page */
      }
   result = &(lookuppage->line[lookuppage->NextFreeEntry]);
   *page = lookuppage;
   return( result ) ;
}

/* lookup name in whiteboard, return elementtype, elementsize, elements, pointers to line and page */
/* if err_not_found is non zero, issue an error message, otherwise be silent and only return error code */
static int c_wb_lookup(WhiteBoard *WB,unsigned char *name, int *elementtype, int *elementsize, int *elements,
                        LINE **result_line, PAGE **result_page, int err_not_found, int lname){
   LINE line;
   int target_page=0;
   int target_line=-1;
   PAGE *lookuppage=NULL;

   if(WB == NULL) WB=BaseWhiteboardPtr;

   TRIM(name,lname)
   if(lname > WB_MAXNAMELENGTH) WB_ERR_EXIT(WB_MSG_ERROR,WB_ERR_NAMETOOLONG);

   lookuppage=WB->PageChain;
   fill_line(&line,name,lname);    /* put key name into local line */
#ifdef DEBUG
printf("target name = :%s:\n",&(line.m.name.c[0]));
#endif
   while( lookuppage != NULL ){
      target_line=lookup_line(&line,lookuppage,0);   /* no screaming if not found, full length match */
      if( target_line >= 0 ){       /* a matching key has been found in this page */
            *elementtype = lookuppage->line[target_line].m.flags.type;
            *elementsize = lookuppage->line[target_line].m.flags.elementsize;
         if( lookuppage->line[target_line].m.flags.array == 1 ) 
            *elements = lookuppage->line[target_line].m.data.a.MaxElements;
         else
            *elements = 0;
         *result_line = &(lookuppage->line[target_line]);
         *result_page = lookuppage;
         return( WB_MAXLINESPERPAGE*target_page + target_line );
         }
      target_page++;
      lookuppage=lookuppage->next;
      }
   if(err_not_found) WB_ERR_EXIT(WB_MSG_ERROR,WB_ERR_NOTFOUND);
   return(WB_ERR_NOTFOUND);
}

static int options_to_flags(LINE *line, int Options){
   line->m.flags.array = (Options & WB_IS_ARRAY) ? 1 : 0 ;
   line->m.flags.fromrestart = (Options & WB_CREATED_BY_RESTART) ? 1 : 0  ;     /* not created by a restart */
   line->m.flags.readonly = (Options & WB_REWRITE_NONE) ? 1 : 0 ;
   line->m.flags.badval = (Options & WB_BADVAL) ? 1 : 0 ;   /* mark as bad value in case it fails unless it is a create only call*/
   line->m.flags.islocal = (Options & WB_IS_LOCAL) ? 1 : 0 ;     /* variable has the local attribute, this will be used when checkpointing */
   line->m.flags.resetuntil = (Options & WB_REWRITE_UNTIL) ? 1 : 0 ;
   line->m.flags.resetmany = (Options & WB_REWRITE_MANY) ? 1 : 0 ;
   line->m.flags.noresetafterrestart = (Options & WB_READ_ONLY_ON_RESTART) ? 1 : 0 ;
   line->m.flags.initialized = (Options & WB_INITIALIZED) ? 1 : 0;
   if( line->m.flags.readonly + line->m.flags.resetmany + line->m.flags.resetuntil > 1 ) WB_ERR_EXIT(WB_MSG_ERROR,WB_ERR_OPTION);
   return(WB_OK);
}

static int flags_to_options(LINE *line){
   int Options=0;
   if(line->m.flags.array) Options += WB_IS_ARRAY;
   if(line->m.flags.fromrestart) Options += WB_CREATED_BY_RESTART;
   if(line->m.flags.readonly) Options += WB_REWRITE_NONE;
   if(line->m.flags.badval) Options += WB_BADVAL;
   if(line->m.flags.islocal) Options += WB_IS_LOCAL;
   if(line->m.flags.resetuntil) Options += WB_REWRITE_UNTIL;
   if(line->m.flags.resetmany) Options += WB_REWRITE_MANY;
   if(line->m.flags.noresetafterrestart) Options += WB_READ_ONLY_ON_RESTART;
   if(line->m.flags.initialized) Options += WB_INITIALIZED;
   else                          Options += WB_NOTINITIALIZED;
   return(Options);
}

/* FORTRAN callable subroutine to get emtadata associated with a whiteboard name */
wordint f77_name(f_wb_get_meta)(WhiteBoard **WB, unsigned char *name, wordint *elementtype, wordint *elementsize,
                                       wordint *elements, wordint *options, F2Cl lname){
   int Elementtype, Elementsize, Elements;
   LINE *line;
   PAGE *page;
   int Lname=lname;
   int status;
/*
   int Options=0;
*/
   status = c_wb_lookup(*WB,name,&Elementtype,&Elementsize,&Elements,&line,&page,0,Lname); /* no screaming if not found */
   if(status < 0)return(status);

   *elementtype = Elementtype;
   *elementsize = Elementsize;
   *elements = Elements;
   *options = flags_to_options(line);
   return(status);
}

void f77_name(f_logical_move)(void *,void *,wordint *);
void f77_name(f_logical2int)(void *,void *,wordint *);
void f77_name(f_int2logical)(void *,void *,wordint *);

/*
  get the data associated with a whiteboard entry 
  name   : pointer to character string containing name of key (length MUST be supplied in Lname)
  Type   : character value R/I/L/C , key type  real/integer/logical/character
  Ltype  : length in bytes of each key element 4/8 for R/I/L, 1->WB_MAXSTRINGLENGTH for character strings
  value  : pointer to where data is to be returned (everything considered as unsigned bytes)
  Nval   : number of elements that can be stored into value
  Lname  : length of key pointed to by name (FORTRAN style)
*/
int c_wb_get(WhiteBoard *WB, unsigned char *name, char Type, int Ltype,unsigned char *value, int Nval, int Lname){
   LINE *line;
   PAGE *page;
   int Elementtype, Elementsize, Elements;
   int status,typecode,array,Result;
   int i;
   unsigned char *cdst=value;
   unsigned char *csrc;

   extra_error_message = " invalid whiteboard instance";
   if(WB == DummyWhiteboardPtr) WB_ERR_EXIT(WB_MSG_ERROR,WB_ERR_NOTFOUND);
   if(WB == NULL) WB=BaseWhiteboardPtr;

   TRIM(name,Lname)
   set_extra_error_message(name,Lname);

   if(Type==WB_FORTRAN_BOOL && Ltype==1) Ltype=4;  /* special case, force 4 byte container for FORTRAN logical type */

   typecode = get_typecode(Type,Ltype,WB_MAXSTRINGLENGTH);  /* check validity of type / ltype combination */
   if(typecode < 0 ) return(typecode);    /* bad type / ltype combination */
   if(Nval < 0)  WB_ERR_EXIT(WB_MSG_ERROR,WB_ERR_WRONGDIMENSION);  /* a negative number of values is a bad idea   */
   array = (Nval > 0) ? 1 : 0;
   if(Nval == 0) Nval = 1;                /* scalars */
   status = c_wb_lookup(WB,name,&Elementtype,&Elementsize,&Elements,&line,&page,0,Lname); /* no screaming if not found */
   if(status < 0) WB_ERR_EXIT(WB_MSG_INFO,WB_ERR_NOTFOUND);
   if( line->m.flags.badval == 1 ) WB_ERR_EXIT(WB_MSG_ERROR,WB_ERR_BADVAL);
   if( line->m.flags.initialized != 1 ) WB_ERR_EXIT(WB_MSG_ERROR,WB_ERR_NOVAL);
   if(status < 0) WB_ERR_EXIT(WB_MSG_ERROR,WB_ERR_NOTFOUND);
   if( Elementtype != typecode ) WB_ERR_EXIT(WB_MSG_ERROR,WB_ERR_WRONGTYPE);  /* types MUST match */
   /* for non character variables, element length must match */
   /* for character variables, the length verification will be performed when copying values */
   if( Elementtype != WB_FORTRAN_CHAR && Elementsize != Ltype ) WB_ERR_EXIT(WB_MSG_ERROR,WB_ERR_WRONGTYPE);
   if(array && line->m.flags.array != 1) WB_ERR_EXIT(WB_MSG_ERROR,WB_ERR_WRONGDIMENSION);
   /* we can now proceed to copy the value(s) from whiteboard entry */
   if( array ) {   /* one dimensional array */
      if(Nval > line->m.data.a.MaxElements) Nval = line->m.data.a.MaxElements ;   /* will do a short get */
      Result = line->m.data.a.MaxElements;   /* if successful return max dimension of array */
      csrc=&((line+1)->d.data[0]);  /* array data is always stored in lines following metadata */
      }
   else {         /* scalar */
      Result = 0;
      csrc = &(line->m.data.c[0]); /* scalar data is stored starting in metadata line */
      }
   if( Elementtype != WB_FORTRAN_CHAR ) {  /* not a character string */
      wordint Nvalf=Nval;
      if ( Elementtype == WB_FORTRAN_BOOL )   /* special case for FORTRAN logicals */
         /*f 77name(f_logical_move)( cdst , csrc, &Nvalf );    FORTRAN helper to move logicals */
         f77_name(f_int2logical)( cdst , csrc, &Nvalf );   /*  FORTRAN helper to move ints into logicals */
      else {
         memcpy( cdst , csrc , Nval*Ltype );   /* straight memory move */
         }
      }
   else {  /* character strings,  use trimmed to padded copy */
      for( i=0 ; i<Nval ; i++ ) {   /* loop over values */
         int tempstat=c_fortran_string_copy(csrc,cdst,Elementsize,Ltype,' ');
         if(tempstat < 0) WB_ERR_EXIT(WB_MSG_ERROR,WB_ERR_WRONGSTRING);  /* OOPS, not enough space to store string */
         cdst += Ltype;
         csrc += Elementsize;
         }
      }
   return(Result);  /* dimension of whiteboard array if array, 0 if scalar, <0 if error */
}
/* FORTRAN callable version of above routine */
wordint f77_name(f_wb_get)(WhiteBoard **WB, unsigned char *name,wordint *type, wordint *ltype,void *value, wordint *nval, F2Cl lname){
   int Lname=lname;
   int Ltype=*ltype;
   int Nval=*nval;
   char Type=*type;
   return( c_wb_get(*WB,name, Type, Ltype, value, Nval, Lname) );
}

static LINE *LastPutLine=NULL;  /* address of the line that was the target of the last put/create operation */
                              /*  used by wb_read to avoid a costly lookup */
/*
  put entry into whiteboard
  name   : pointer to character string containing name of key (length MUST be supplied in Lname)
  Type   : character value R/I/L/C , key type  real/integer/logical/character
  Ltype  : length in bytes of each key element 4/8 for R/I/L, 1->WB_MAXSTRINGLENGTH for character strings
  value  : pointer to data associated with key (everything considered as unsigned bytes)
  Nval   : number of elements (0 means a scalar) (1 or more means an array)
  Options: options associated with name
  Lname  : length of key pointed to by name (FORTRAN style)
*/
int c_wb_put(WhiteBoard *WB, unsigned char *name,char Type,int Ltype,unsigned char *value,int Nval,int Options,int Lname){
   LINE *line;
   PAGE *page;
   int Elementtype, Elementsize, Elements;
   int status,typecode,lines,array,LinesAvailable,Result;
   int i;
   unsigned char *csrc=value;
   unsigned char *cdst;

   extra_error_message = " invalid whiteboard instance";
   if(WB == DummyWhiteboardPtr) WB_ERR_EXIT(WB_MSG_ERROR,WB_ERR_NOTFOUND);
   if(WB == NULL) WB=BaseWhiteboardPtr;
   extra_error_message = "";

   TRIM(name,Lname)
   set_extra_error_message(name,Lname);
   LastPutLine=NULL;
   if(Type==WB_FORTRAN_BOOL && Ltype==1) Ltype=4;  /* special case, force 4 byte container for FORTRAN logical type */

   typecode = get_typecode(Type,Ltype,Ltype);  /* check validity of requested Type/Ttype combination */
   if(typecode < 0) return(typecode);    /* bad Type/Ltype combination */
   if(Nval < 0)  WB_ERR_EXIT(WB_MSG_ERROR,WB_ERR_WRONGDIMENSION);  /* a negative number of values is a bad idea   */
   array = (Nval > 0) ? 1 : 0;       /* array if Nval > 0 , scalar if Nval == 0 */
   if(Nval == 0) Nval = 1;           /* scalars */
   status = c_wb_lookup(WB,name,&Elementtype,&Elementsize,&Elements,&line,&page,0,Lname); /* no screaming if not found */

   if(status >= 0) {  /* entry has been found, we want to rewrite it, the options argument is ignored,
                         whiteboard options are checked for consistency and permissions */
      if(Options & WB_CREATE_ONLY)  WB_ERR_EXIT(WB_MSG_ERROR,WB_ERR_REDEFINE); /* redefinition not allowed */
/* above error has to be ignored if it was coming from a restart and create from restart flag then has to be erased */
      if(line->m.flags.readonly && line->m.flags.initialized )
          WB_ERR_EXIT(WB_MSG_ERROR,WB_ERR_READONLY);        /* entry is READONLY, OOPS */
      line->m.flags.badval = 1;  /* mark as bad value in case it fails */
      if( line->m.flags.array==1 && array!=1 )
         WB_ERR_EXIT(WB_MSG_ERROR,WB_ERR_WRONGDIMENSION);  /* array to scalar or scalar to array is a NO NO */
      if( Elementtype != typecode ) WB_ERR_EXIT(WB_MSG_ERROR,WB_ERR_WRONGTYPE);  /* types MUST match */
      /* for non character variables, element length must match */
      /* for character variables, the length verification will be performed when copying values */
      if( Elementtype != WB_FORTRAN_CHAR && Elementsize != Ltype ) WB_ERR_EXIT(WB_MSG_ERROR,WB_ERR_WRONGTYPE);
      if( array &&  line->m.data.a.MaxElements < Nval )
         WB_ERR_EXIT(WB_MSG_ERROR,WB_ERR_WRONGDIMENSION);  /* array in whiteboard is too small */
      }
   else{   /* there was an error during lookup */
      if(status != WB_ERR_NOTFOUND) return(status); /* if ERROR = not found, create entry, otherwise, return error code */
      lines = 1;                        /* number of storage lines is 1 for scalars using less than 9 bytes */
      if(array || typecode==WB_FORTRAN_CHAR )    /* arrays and strings need supplementary storage space */
         lines += ((Nval * Ltype) + sizeof(DATALINE) - 1 ) / sizeof(DATALINE) ;
      line = new_line(WB,&page,lines);              /* create new entry with space for lines  */
      if(line == NULL) return(WB_ERR_ALLOC);     /* no need to use error macro, it will have already been called */
      fill_line(line, name, Lname);              /* put key name into line */
      Elementtype = typecode;
      Elementsize = Ltype;
      if( Elementtype==WB_FORTRAN_CHAR && array==0 ) {   /* scalar string, round size up to DATALINE size */
         Elementsize = ( (Elementsize+sizeof(DATALINE)-1) / sizeof(DATALINE) ) * sizeof(DATALINE) ; 
/*         Ltype = Elementsize;  */
         }
      line->m.flags.type = Elementtype;           /* element type and size */
      line->m.flags.elementsize = Elementsize;
      line->m.flags.lines = ( lines > WB_MAXLINESPERPAGE ) ? 0 : lines ;   /* only one entry in page if lines>WB_MAXLINESPERPAGE */
      line->m.flags.readonly = (Options & WB_REWRITE_NONE) ? 1 : 0 ;
      line->m.flags.resetmany = (Options & WB_REWRITE_MANY) ? 1 : 0 ;
      line->m.flags.resetuntil = (Options & WB_REWRITE_UNTIL) ? 1 : 0 ;
      if( line->m.flags.readonly + line->m.flags.resetmany + line->m.flags.resetuntil > 1 ) WB_ERR_EXIT(WB_MSG_ERROR,WB_ERR_OPTION);
      line->m.flags.initialized = 0;
      line->m.flags.islocal = (Options & WB_IS_LOCAL) ? 1 : 0 ;     /* variable has the local attribute, this will be used when checkpointing */
      line->m.flags.badval = (Options & WB_CREATE_ONLY) ? 0 : 1 ;   /* mark as bad value in case it fails unless it is a create only call*/
      line->m.flags.fromrestart = 0 ;     /* not created by a restart */
      line->m.flags.noresetafterrestart = (Options & WB_READ_ONLY_ON_RESTART) ? 1 : 0 ;
      line->m.flags.array = array ;
      if(array) {                 /* arrays are not explicitly initialized for the time being */
         line->m.data.a.UsedElements = Nval;
         line->m.data.a.MaxElements = Nval;
         }
      else {                      /* scalar, set data to 0 */
         line->m.data.i[0] = 0;
         line->m.data.i[1] = 0;
         }
      LinesAvailable = page->LinesInPage - page->NextFreeEntry ;             /* check if there is enough space in page */
      if( LinesAvailable < lines ) WB_ERR_EXIT(WB_MSG_ERROR,WB_ERR_NOMEM);   /* this should never happen */
      page->NextFreeEntry += lines;   /* allocate space for data, adjust page pointer to next available entry */
      }
   LastPutLine=line;
   if(Options & WB_CREATE_ONLY) {   /* create entry only, do not copy value(s)  */
      if(array)
         Result = line->m.data.a.MaxElements;
      else
         Result = 0;
      }
   else{   /* we can now proceed to copy the value(s) into (old or new) whiteboard entry */
      if( array ) {   /* one dimensional array */
         if(Nval > line->m.data.a.MaxElements) WB_ERR_EXIT(WB_MSG_ERROR,WB_ERR_WRONGDIMENSION);  /* OOPS, not enough space to store values */
         Result = line->m.data.a.MaxElements;   /* if successful return max dimension of array */
         cdst=&((line+1)->d.data[0]);  /* array data is always stored in lines following metadata */
         }
      else {         /* scalar */
         Result = 0;
         cdst = &(line->m.data.c[0]); /* scalar data is stored starting in metadata line */
         }
      if( Elementtype != WB_FORTRAN_CHAR ) {  /* not a character string */
         wordint Nvalf=Nval;
         if ( Elementtype == WB_FORTRAN_BOOL )   /* special case for FORTRAN logicals */
            /* f77_name(f_logical_move)( cdst , csrc, &Nval );   FORTRAN helper to move logicals */
            f77_name(f_logical2int)( cdst , csrc, &Nvalf );   /* FORTRAN helper to move logicals into ints */
         else {
            memcpy( cdst , csrc , Nval*Ltype );   /* straight memory move */
            }
         }
      else {  /* character strings,  use trimmed to padded copy */
         for( i=0 ; i<Nval ; i++ ) {   /* loop over values */
            Result=c_fortran_string_copy(csrc,cdst,Ltype,Elementsize,' ');
            if(Result < 0) WB_ERR_EXIT(WB_MSG_ERROR,WB_ERR_WRONGSTRING);  /* OOPS, not enough space to store string */
            csrc += Ltype;
            cdst += Elementsize;
            }
         }
      line->m.flags.initialized = 1 ;  /* mark entry as initialized */
      line->m.flags.badval = 0;        /* and good */
      }
   /* return <0 error code if there was an error during the copy of value(s) into whiteboard */
   /* return 0 for non character scalars, whiteboard string length if character scalar, max size of array if array */
   return (Result);
}
/* FORTRAN version of c_wb_put */
wordint f77_name(f_wb_put)(WhiteBoard **WB, unsigned char *name,wordint *type, wordint *ltype,void *value, wordint *nval, wordint *options, F2Cl lname){
   char Type=*type;
   int Lname=lname;
   int Ltype=*ltype;
   int Nval = *nval;
   int Options = *options;
   return(c_wb_put(*WB,name,Type,Ltype,value,Nval,Options,Lname));
}

/* write checkpoint file */
int c_wb_checkpoint(){
   PAGE *lookuppage;
   int pageno=0;
   int zero=0;
   int status;
   int fd=open(WhiteBoardCheckpointFile,O_WRONLY|O_CREAT,0777);

   if(fd < 0) WB_ERR_EXIT(WB_MSG_ERROR,WB_ERR_CKPT);  /* OOPS, cannot open checkpoint file */
   lookuppage = BaseWhiteboardPtr->PageChain;

   status = write(fd,"WBckp100",8);                  /* write signature */
   if(status < 0) { close(fd); WB_ERR_EXIT(WB_MSG_ERROR,WB_ERR_CKPT); }      /* write error */

   while( lookuppage != NULL ) {
      status = write(fd,&(lookuppage->LinesInPage),4); /* write number of lines in page */
      if(status < 0) { close(fd); WB_ERR_EXIT(WB_MSG_ERROR,WB_ERR_CKPT); }      /* write error */
      status = write(fd,&(lookuppage->NextFreeEntry),4); /* write number of next entry */
      if(status < 0) { close(fd); WB_ERR_EXIT(WB_MSG_ERROR,WB_ERR_CKPT); }      /* write error */
      status = write(fd,&(lookuppage->line[0]),sizeof(LINE)*lookuppage->LinesInPage);  /* write page */
      if(status < 0) { close(fd); WB_ERR_EXIT(WB_MSG_ERROR,WB_ERR_CKPT); }      /* write error */
      if( message_level <= WB_MSG_DEBUG )
         fprintf(stderr,"wb_checkpoint: Page %d, length %d lines, Next entry %d \n",pageno,lookuppage->LinesInPage,lookuppage->NextFreeEntry);
      lookuppage=lookuppage->next;
      pageno++;
      }
   status = write(fd,&zero,4) ;   /* 0 length page is the end marker */
   close(fd);
   return(0);
}
/* FORTRAN version of above */
wordint f77_name(f_wb_checkpoint)(WhiteBoard **WB){
   return(c_wb_checkpoint());
}

/*
 basic whiteboard check/action routine, the Whiteboard "swiss knife"
  name   : pointer to character string containing name of key (length MUST be supplied in Lname)
  OptionsMask: options mask to be tested, call Action if OptionsMask&options  is non zero
  Lname  : length of key pointed to by name (FORTRAN style)
  printflag: print messages if nonzero
  Action:  routine to be called if name and OptionsMask &options  is non zero
  blinddata: pointer to be passed to Action routine as a second argument (first argument is matching line)
*/

typedef int (*ACTION)(LINE *,void *);
int c_wb_check(WhiteBoard *WB, unsigned char *name, int OptionsMask, int Lname, int printflag, ACTION Action, void *blinddata ){
   PAGE *lookuppage;
   int pageno=0;
   int printcount=0;
   int status;

   extra_error_message = " invalid whiteboard instance";
   if(WB == DummyWhiteboardPtr) WB_ERR_EXIT(WB_MSG_ERROR,WB_ERR_NOTFOUND);
   if(WB == NULL) WB=BaseWhiteboardPtr;

   TRIM(name,Lname)
   lookuppage = WB->PageChain;

   while( lookuppage != NULL ){
      int i,linelimit;
      i=0;
      linelimit=lookuppage->NextFreeEntry - 1;
      while( i <= linelimit ){
         int Options;
         LINE *line=&(lookuppage->line[i]);
         unsigned char *c1=&(line->m.name.c[0]);
         Options = flags_to_options(line);
         if( (Options & OptionsMask) && wb_match_name(c1,name,Lname) ) {
            printcount--;
            if(printflag){
               fprintf(stderr,"Page %2d, Line %3d, KEY=%s Datatype=%s[%4d] Size=%4d %s %s %s %s %s %s %s %s\n",
                     pageno,i,c1,
                     datatypes[line->m.flags.type],
                     line->m.flags.elementsize,
                     (line->m.flags.array) ? line->m.data.a.MaxElements : 0,
                     (line->m.flags.resetuntil) ? "LOCKABLE" : "        ",
                     (line->m.flags.initialized) ? "SET  " : "UNSET",
                     (line->m.flags.badval) ? "BAD " : "GOOD",
                     (line->m.flags.resetmany) ? "NOTLOCKABLE" : "           ",
                     (line->m.flags.noresetafterrestart) ? "ROonRst" : "RWonRst",
                     (line->m.flags.readonly) ? "LOCKED  " : "WRITABLE",
                     (line->m.flags.islocal) ? "LOCAL " : "GLOBAL",
                     (line->m.flags.fromrestart) ? "FromRestart" : "FromPut"
                     );
               }
            if(Action != NULL) {
               status = (*Action)(line,blinddata);
               if(status < 0) return(status);
               }
            }
         if(line->m.flags.lines == 0) break;   /* only one entry in this page in this case */
         i += line->m.flags.lines;
         }
      lookuppage=lookuppage->next;
      pageno++;
      }
return(-printcount);
}
wordint f77_name(f_wb_check)(WhiteBoard **WB, unsigned char *name, wordint *optionsmask, F2Cl lname){
   int OptionsMask=*optionsmask;
   int Lname=lname;
   /* print flag on, no action routine is supplied, no blind data pointer is supplied */
   return( c_wb_check(*WB, name, OptionsMask, Lname, 1, NULL, NULL) );
}

/* action routine to make a variable read-only */
static int MakeReadOnly(LINE *line, void *blinddata){
   if(line->m.flags.initialized == 0) WB_ERR_EXIT(WB_MSG_ERROR,WB_ERR_NOVAL); /* making read-only a non initialized variable ???!!! */
   if( message_level <= WB_MSG_DEBUG ) 
      fprintf(stderr,"key '%s' is now read-only\n",&(line->m.name.c[0]));
   line->m.flags.readonly = 1; line->m.flags.resetmany = 0;
   return(0);
}

/* action routine to make a variable no longer read-only */
static int MakeNotReadOnly(LINE *line, void *blinddata){
   if( message_level <= WB_MSG_DEBUG ) 
      fprintf(stderr,"key '%s' is now no longer read-only\n",&(line->m.name.c[0]));
   line->m.flags.readonly = 0;
   return(0);
}

/* action routine to make a variable undefined (no longer initialized) */
static int MakeUndefined(LINE *line, void *blinddata){
   if( message_level <= WB_MSG_DEBUG ) 
      fprintf(stderr,"key '%s' is now undefined\n",&(line->m.name.c[0]));
   line->m.flags.initialized = 0;
   line->m.flags.badval = 0;
   return(0);
}

/* action routine to make a variable undefined (no longer initialized) */
static int MakeCreatedByRestart(LINE *line, void *blinddata){
   if( message_level <= WB_MSG_DEBUG ) 
      fprintf(stderr,"key '%s' is now marked as created by a restart\n",&(line->m.name.c[0]));
   line->m.flags.fromrestart = 1;
   return(0);
}

/* framework for broadcasting values across MPI tiles */
typedef struct{
   wordint pe_root;
   wordint pe_me;
   char *domain;
   F2Cl ldomain;
   void (*broadcast_function)();
   void (*allreduce_function)();
   } BROADCAST;
static BROADCAST broadcast_table;

void f77_name(f_wb_bcst_init)(wordint *pe_root, wordint *pe_me, char *domain, void (*broadcast_function)(), void (*allreduce_function)() ){
   broadcast_table.pe_root = *pe_root;
   broadcast_table.pe_me = *pe_me;
   broadcast_table.domain = domain;
   broadcast_table.broadcast_function = broadcast_function;
   broadcast_table.allreduce_function = allreduce_function;
}

static int BroadcastLine(LINE *line, void *blinddata){
fprintf(stderr,"Broadcasting %s :-)\n",line->m.name.c);
   return(0);
}

wordint f77_name(f_wb_bcst)(WhiteBoard **WB, unsigned char *name, wordint *wildcard, F2Cl lname){
   int Lname = lname;
   LINE line;
   int i,status;

   fill_line(&line,name,Lname);
   if(*wildcard != 0) {   /* wildcard match, get rid of trailing blanks */
      for( i=0 ; i<WB_MAXNAMELENGTH-1 ; i++ ) {
         if(line.m.name.c[i]==' ') line.m.name.c[i]=0;
         }
      }
   status = c_wb_check(*WB, name, -1, Lname, message_level<=WB_MSG_INFO, BroadcastLine, &broadcast_table) ;
   return(status);
}

typedef struct{
   unsigned char *UserKeyArray;
   int UserKeyLength;
   int UserMaxLabels;
   } KEYS;

/* action routine used to copy keys into user table when collecting a list of keys */
static int CopyKeyName(LINE *line, void *blinddata){
   unsigned char *key=&(line->m.name.c[0]);
   KEYS *keys=(KEYS *)blinddata;
   if(keys->UserMaxLabels <= 0) return(WB_ERROR);   /* too many keys for user array */
   c_fortran_string_copy(key, keys->UserKeyArray, WB_MAXNAMELENGTH, keys->UserKeyLength, ' ');
   (keys->UserMaxLabels)--;
   keys->UserKeyArray += keys->UserKeyLength;
   return(0);
}

/* get a list of keys matching a name */
wordint f77_name(f_wb_get_keys)(WhiteBoard **WB, unsigned char *labels, wordint *nlabels, unsigned char *name, F2Cl llabels, F2Cl lname){
   int Lname = lname;
   int status;
   KEYS keys;

   keys.UserKeyArray = labels;
   keys.UserMaxLabels = *nlabels;
   keys.UserKeyLength = llabels;
   status = c_wb_check(*WB, name, -1, Lname, message_level<=WB_MSG_INFO, CopyKeyName, &keys) ;
   return (status);
}

/* make readonly all entries that have the WB_REWRITE_UNTIL attribute */
int c_wb_lock(WhiteBoard *WB, unsigned char *name, int Lname){
   unsigned char localname[WB_MAXSTRINGLENGTH];
   int llname;
   llname=c_fortran_string_copy(name,localname,Lname,sizeof(localname),'\0');
   if( message_level <= WB_MSG_DEBUG )
      fprintf(stderr,"locking variables with name beginning with '%s'\n",localname);
   return( c_wb_check(WB, localname, WB_REWRITE_UNTIL, llname, 1, MakeReadOnly, NULL) );
}
wordint f77_name(f_wb_lock)(WhiteBoard **WB, unsigned char *name, F2Cl lname){
   int Lname=lname;
   int status = c_wb_lock(*WB, name,Lname);
   return(status);
}

/* write whiteboard checkpoint file */
int c_wb_reload(){
   int pageno=0;
   int pagelen;
   int status;
   char signature[9];
   int fd=open(WhiteBoardCheckpointFile,O_RDONLY);

   if(fd < 0) WB_ERR_EXIT(WB_MSG_ERROR,WB_ERR_CKPT);  /* OOPS, cannot open checkpoint file */
   if(BaseWhiteboardPtr->PageChain != NULL) WB_ERR_EXIT(WB_MSG_ERROR,WB_ERROR);
   status=read(fd,signature,8) ; signature[8]=0;
   if( status!=8 || 0!=strncmp(signature,"WBckp100",8) ) { close(fd); WB_ERR_EXIT(WB_MSG_ERROR,WB_ERR_CKPT); }

   status=read(fd,&pagelen,4);     /* page size */
   while(pagelen>0){
      status=new_page(BaseWhiteboardPtr, pagelen);                    /* allocate new page */
      status=read(fd,&(BaseWhiteboardPtr->PageChain->NextFreeEntry),4);  /* next usable entry in page */
      if( message_level <= WB_MSG_DEBUG )
         fprintf(stderr,"wb_reload: Page %d, length=%d lines, next entry=%d\n",
                 pageno,BaseWhiteboardPtr->PageChain->LinesInPage,BaseWhiteboardPtr->PageChain->NextFreeEntry);
      status=read(fd,&(BaseWhiteboardPtr->PageChain->line[0]),sizeof(LINE)*pagelen);  /* read page */
      status=read(fd,&pagelen,4);     /* page size of next page, 0 means no more */
      pageno++;
      }
   close(fd);
/* if variable is WB_REWRITE_AT_RESTART, erase the read-only flag */
   status=c_wb_check(BaseWhiteboardPtr, (unsigned char*)"", WB_REWRITE_AT_RESTART, 0, 0, MakeNotReadOnly, NULL);
/* if variable is WB_READ_ONLY_ON_RESTART, make sure it is now marked as read-only */
   status=c_wb_check(BaseWhiteboardPtr, (unsigned char*)"", WB_READ_ONLY_ON_RESTART, 0, 0, MakeReadOnly, NULL);
/* if variable is local (not the same on ALL MPI whiteboards, mark it as not initialized */
   status=c_wb_check(BaseWhiteboardPtr, (unsigned char*)"", WB_IS_LOCAL, 0, 0, MakeUndefined, NULL);
/* mark all variables as having been created by a restart */
   status=c_wb_check(BaseWhiteboardPtr, (unsigned char*)"", -1, 0, 0, MakeCreatedByRestart, NULL);
   if( message_level<=WB_MSG_INFO )
      fprintf(stderr,"whiteboard has been reloaded, variables read only after restart have been locked\n");
   return(WB_OK);
}
wordint f77_name(f_wb_reload)(WhiteBoard **WB){
   int status=c_wb_reload();
   return(status);
}

static int Action1(LINE *line, void *blinddata){   /* DUMMY ACTION ROUTINE */
   fprintf(stderr,"Action1 has been called, key=%s\n",&(line->m.name.c[0]));
   return(0);
}

#define MISC_BUFSZ 1024
/* get a line from file, reset current character pointer */
static int wb_get_line(FILE *infile){
   current_char=fgets(linebuffer, MISC_BUFSZ, infile);
   if( message_level<=WB_MSG_INFO && current_char )
      fprintf(stderr,">>%s",linebuffer);
   if(current_char) return(*current_char);
   else             return(EOF);
}

/* push character back onto input buffer */
static int wb_ungetc(int c){
   if( current_char==NULL || linebuffer==NULL ) return(WB_ERROR);
   if( current_char > linebuffer ) {
      current_char--;
      *current_char = c;
      }
   else {
      return(WB_ERROR);
      }
   return(WB_OK);
}

/* get next character from input stream, take care of \newline sequence */
static int wb_getc(FILE *infile){
   int c;

   if( current_char==NULL ) {                  /* no buffer pointer, fill buffer */
      c=wb_get_line(infile);
      if( current_char==NULL ) return(EOF);    /* End Of File or error */
      }
s: if( *current_char==0 ) {
      c=wb_get_line(infile);           /* try to refill empty buffer */
      if( current_char==NULL ) return(EOF);    /* End Of File or error */
      }
   c=*current_char++ ;                         /* get next character and bump pointer */
   if( c=='\\' && *current_char=='\n' ) {      /* \newline, get rid of it */
      current_char++;                          /* point to character after newline */
      goto s;                                  /* and start all over again */
      }
   return(c);
}

/* flush input until newline or EOF , return what was found */
static int wb_flush_line(FILE *infile){
   int c=wb_getc(infile);
   while( c!='\n' && c!=EOF && c!=';' ) c=wb_getc(infile);
   return (c);
}

/* skip blanks and tabs from input, return next character, newline, or EOF ,
    take care of comments, recognize ; as newline    */
static int wb_get_nonblank(FILE *infile){
   int c=wb_getc(infile);
  if( c==EOF ) return(EOF);
   while( (c==' ' || c=='\t' || c=='#' || c=='!' ) && c!=EOF ) {  /* space, tab, comment character */
      if( c=='#' || c=='!' ) {             /* comment, get rid of the rest of the line */
         c=wb_flush_line(infile);
         if( c==EOF ) return(EOF);
         }
      c=wb_getc(infile);           /* get next character */
      }
   if( c==';' ) c='\n' ;      /* treat a non quoted ; as a newline */
   return (c);
}

/* special error print when reading directives, print directive line up to error, add ^^^^ error marker, print rest of line */
static void wb_read_error(){
   int temp;
   if(current_char==NULL || linebuffer==NULL) return;
   temp=*current_char;
   *current_char=0;
   if( message_level<=WB_MSG_ERROR )
      fprintf(stderr,"ERROR in directives, offending line:\n%s^^^^",linebuffer);
   *current_char=temp;
   if( message_level<=WB_MSG_ERROR )
      fprintf(stderr,"%s\n",current_char);
}

/* get next token (alphanum_, number, delimited string, single character, force non quoted tokens to uppercase */
static int wb_get_token(unsigned char *token,FILE *infile,int maxtoken, int noskip){
   int ntoken=0;
   int c;

   if(noskip)  /* do not skip spaces */
      c=wb_getc(infile);
   else
      c=wb_get_nonblank(infile);
   if( c==EOF ) return(WB_ERROR);
   token[ntoken++] = toupper(c);  /* first character */
   maxtoken--; if(maxtoken<0) goto error_big;
   if(isalpha(c)) {   /* collect alphanum _ token */
      for( c=wb_getc(infile) ; isalnum(c) || c=='_' ; c=wb_getc(infile) ) {
         maxtoken--; if(maxtoken<0) goto error_big;
         token[ntoken++]=toupper(c) ;
         }
      if( c==EOF ) return(WB_ERROR);
      if(wb_ungetc(c)==EOF ) return(WB_ERROR);  /* push back extra input character */
      }
   else if( c=='\'' || c=='"' ) {         /* collect ' or " delimited string */
      int quote=c;
      c=wb_getc(infile);
      while( c!=quote && c!='\n' && c!=EOF) {   /* look for matching quote , error end if newline/EOF  encountered */
         if(c==EOF || c=='\n') break;
         maxtoken--; if(maxtoken<0) goto error_big;
         token[ntoken++]=c;
         c=wb_getc(infile);
         }
      if( c=='\n' ) {
         if( message_level<=WB_MSG_ERROR ) fprintf(stderr,"ERROR: improperly terminated string\n");
         if(wb_ungetc(c)==EOF ) return(WB_ERROR); /* push back newline that has been erroneously swallowed */
         }
      if(c==EOF || c=='\n') return(WB_ERROR);
      maxtoken--; if(maxtoken<0) goto error_big;
      token[ntoken++]=c;                       /* store end delimiter in token */
      }
   else if(isdigit(c) || c=='.' || c=='+' || c=='-' ) {  /* digit, point, sign, potential number or boolean */
      c=wb_getc(infile) ;
      if( isalpha(c) && token[ntoken-1]=='.' ) {    /* collect .true. , .false. , etc ... */
         while( isalpha(c) || c=='.' ) {
            maxtoken--; if(maxtoken<0) goto error_big;
            token[ntoken++]=toupper(c) ;\
            c=wb_getc(infile);
            }
         }
      else {    /* collect a potential number */
         while( isdigit(c) || c=='.' || c=='-' || c=='+' || c=='E' || c=='e' ) {
            maxtoken--; if(maxtoken<0) goto error_big;
            token[ntoken++]=toupper(c) ;
            c=wb_getc(infile);
            }
         }
      if( c==EOF ) return(WB_ERROR);
      if(wb_ungetc(c)==EOF ) return(WB_ERROR);  /* push back extra input character */
      }
   /* none of the above means single character */
   maxtoken--; if(maxtoken<0) goto error_big;
   token[ntoken]=0;  /* null terminate token */

   if( message_level<=WB_MSG_DEBUG )
      fprintf(stderr,"GetToken, maxtoken=%d, ntoken=%d, nospkip=%d, token='%s'\n",maxtoken,ntoken,noskip,token);
   return(ntoken);
error_big:    /* token is too big to be stored in supplied array */
   wb_read_error() ;
   WB_ERR_EXIT(WB_MSG_ERROR,WB_ERR_BIG);
}

/* process options set [option,option,...] or (option,option,...) */
static int wb_options(FILE *infile,char start_delim,char end_delim,SYMTAB *option_table){
   unsigned char Token[WB_MAXNAMELENGTH+1];
   int options=0;
   int ntoken,newoption;

   ntoken=wb_get_token(Token,infile,2,0) ; /* get starting delimiter ( or [  */
   if(Token[0]!=start_delim || ntoken!=1) goto error_syntax ;

   while(Token[0] != end_delim ) {
      ntoken=wb_get_token(Token,infile,sizeof(Token),0) ; /* expect keyname  */
      if( !isalpha(Token[0]) || ntoken<=0 ) goto error_syntax ;
      newoption = wb_value((char *)Token,option_table);
      if( newoption == 0 ) goto error_option;
      options += newoption;
      if( message_level<=WB_MSG_DEBUG )
         fprintf(stderr,"newoption=%d, options=%d\n",newoption,options);
      ntoken=wb_get_token(Token,infile,2,0) ; /* expect end delimiter ( or ]  or comma ,  */
      if( (Token[0]!=end_delim && Token[0]!=',') || ntoken!=1) goto error_syntax ;
      }
   return(options);
error_option:
   wb_read_error() ;
   wb_flush_line(infile);
   WB_ERR_EXIT(WB_MSG_ERROR,WB_ERR_OPTION);
error_syntax:
   wb_read_error() ;
   wb_flush_line(infile);
   WB_ERR_EXIT(WB_MSG_ERROR,WB_ERR_SYNTAX);
}

/* process define directive  define(key_name[TYPE,array_size], ... , ... )  TYPE=R4/R8/I4/I8/L1/Cn  */
static int wb_define(WhiteBoard *WB, FILE *infile, char *package){
   unsigned char name[WB_MAXNAMELENGTH+1];
   unsigned char type[WB_MAXNAMELENGTH+1];
   unsigned char length[WB_MAXNAMELENGTH+1];
   unsigned char Token[MISC_BUFSZ];
   int ntoken;
   int key_type,len_type,nread,array_length,options=default_dict_option;
   char key_c;
   int opt_cmd=0;
   int desc_cmd=0;
   int status=WB_OK;

   ntoken=wb_get_token(Token,infile,2,0) ; /* get (  */
   if(Token[0]!='(' || ntoken!=1) goto error_syntax ;

   strncpy((char *)name,package,strlen(package));
   ntoken=wb_get_token(name+strlen(package), infile, sizeof(name)-strlen(package), 0) ; /* expect keyname  */
   if( !isalpha(name[strlen(package)]) || ntoken<=0 ) goto error_syntax ;

   ntoken=wb_get_token(Token,infile,2,0) ; /* expect [  */
   if(Token[0]!='[' || ntoken!=1) goto error_syntax ;

   ntoken=wb_get_token(type,infile,sizeof(type),0) ; /* expect entry type and length  */
   /*if( !isalpha(type[0]) || ntoken<=0 ) goto error_syntax ;   validate type and length */
   sscanf((char *)type,"%c%d%n",&key_c,&len_type,&nread);
   switch(key_c) {          /* get type code */
      case 'R' : key_type=WB_FORTRAN_REAL ; break ;
      case 'I' : key_type=WB_FORTRAN_INT ; break ;
      case 'C' : key_type=WB_FORTRAN_CHAR ; break ;
      case 'L' : key_type=WB_FORTRAN_BOOL ; break ;
      default: goto error_syntax ;
      }

   if(nread != ntoken ) goto error_syntax ;
   if( (key_type=get_typecode(key_type,len_type,len_type)) < 0 ) goto error_syntax ;  /* check type and length combination */

   ntoken=wb_get_token(Token,infile,2,0) ; /* expect , */
   if(Token[0]!=',' || ntoken!=1) goto error_syntax ;

   ntoken=wb_get_token(length,infile,sizeof(length),0) ; /* expect array length  */
   if( !isdigit(length[0]) || ntoken<=0 ) goto error_syntax ;
   sscanf((char *)length,"%d%n",&array_length,&nread);
   if(nread != ntoken ) goto error_syntax ;

   ntoken=wb_get_token(Token,infile,2,0) ; /* expect ]  */
   if(Token[0]!=']' || ntoken!=1) goto error_syntax ;

   ntoken=wb_get_token(Token,infile,2,0) ; /* expect ) or ,  */
   if( (Token[0]!=')' && Token[0]!=',') || ntoken!=1) goto error_syntax ;

   if( message_level<=WB_MSG_DEBUG )
      fprintf(stderr,"type=%c,key_type=%d,len=%d,array_length=%d\n",key_c,key_type,len_type,array_length);
   while(Token[0] != ')' ) {
      ntoken=wb_get_token(Token,infile,WB_MAXNAMELENGTH,0) ; /* expect subcommand name OPT or DESC  */
      if( strcmp("OPT",(char *)Token)==0 && opt_cmd==0 ) {   /* OPT= subcommand */
         ntoken=wb_get_token(Token,infile,2,0) ; /* expect =  */
         if(Token[0]!='=' || ntoken!=1) goto error_syntax ;
         options =  wb_options(infile,'[',']',dict_options);
         opt_cmd++;
         }
      else if( strcmp("DESC",(char *)Token)==0 && desc_cmd==0 ) {   /*DESC= subcommand */
         ntoken=wb_get_token(Token,infile,2,0) ; /* expect =  */
         if(Token[0]!='=' || ntoken!=1) goto error_syntax ;
         ntoken=wb_get_token(Token,infile,MISC_BUFSZ,0) ; /* expect quoted string  */
         if(Token[0]!='\'' && Token[0]!='"') goto error_syntax ;
         desc_cmd++;  /* DESC is ignored for now, just counted */
         }
      else {
         if( message_level<=WB_MSG_ERROR )
            fprintf(stderr,"invalid/undefined define subcommand\n");
         goto error_syntax ;
         }
      ntoken=wb_get_token(Token,infile,WB_MAXNAMELENGTH,0) ; /* expect ) or ,  */
      if( (Token[0]!=')' && Token[0]!=',') || ntoken!=1) goto error_syntax ;
      }
   status=c_wb_put(WB,name,key_type,len_type,Token,array_length,options|WB_CREATE_ONLY,strlen((char *)name));
   if(status<0) goto error_syntax;                 /* put failed for some reason */
   status=wb_define_check(LastPutLine,1);  /* must not already be in table */
   if(status<0) goto error_syntax;                 /* already defined in table, OOPS */
   wb_flush_line(infile);
   return(status);
error_syntax:
   wb_read_error() ;
   wb_flush_line(infile);
   WB_ERR_EXIT(WB_MSG_ERROR,WB_ERR_SYNTAX);
}

/* process key= directive from file Token contains key name, package contains package prefix string */
static int wb_key(WhiteBoard *WB, FILE *infile,unsigned char *Token, char *package, int Options){
   int status=WB_OK;
   int err_not_found=1;
   int elementtype,elementsize,elements,ntoken,nread,items;
   LINE *line;
   PAGE *page;
   unsigned char name[WB_MAXNAMELENGTH+1];
   unsigned char separator[3];
   unsigned char buffer[WB_MAXSTRINGLENGTH+3];  /* add the two delimiters and the null terminator */
   unsigned char *dataptr;

   strncpy((char *)name,package,strlen(package));   /* concatenate package name with keyname */
   strncpy((char *)name+strlen(package),(char *)Token,sizeof(name)-strlen(package));
   name[WB_MAXNAMELENGTH] = 0;

   if(message_level<=WB_MSG_INFO) fprintf(stderr,"Assigning to '%s'\n",name);
   extra_error_message=name;
   ntoken=wb_get_token(Token,infile,2,0) ; /* expect =  */
   if(Token[0]!='=' || ntoken!=1) {
      if( message_level<=WB_MSG_ERROR )
         fprintf(stderr,"= sign not found where expected \n");
      goto error_syntax ;
      }

   status = c_wb_lookup(WB,name,&elementtype,&elementsize,&elements,&line,&page,err_not_found,strlen((char *)name));
   if(status<0) goto error_syntax;  /* key MUST be defined */

   if( Options==WB_STRICT_DICTIONARY )
      status = wb_define_check(line,0);    /* must be found in table if in strict dictionary mode */
   else
      status = wb_define_check(line,2);    /* must be inserted in table if in normal mode */
   if( status<0 ) goto error_syntax;               /* other source of error : key already assigned a value in this file */

   if(elements == 0) elements = 1;           /* scalar, 1 item */
      if(line->m.flags.array)                /* array data */
         dataptr = &((line+1)->d.data[0]);
      else                                   /* scalar data */
         dataptr = &(line->m.data.c[0]);
   while(elements--){                        /* get values for up to max items for key */
      ntoken=wb_get_token(buffer,infile,sizeof(buffer)-1,0);
      if( ntoken <= 0 ) goto error_syntax ;

      switch(elementtype){
         case(WB_FORTRAN_INT):
            if(elementsize==4){              /* 4 byte integer */
               int *target = (int *)dataptr;
               items = sscanf((char *)buffer,"%d%n",target,&nread);
               if(nread != ntoken || items != 1) goto error_syntax ;
               }
            else if(elementsize==8){         /* 8 byte integer */
               long long *target = (long long *)dataptr;
               items = sscanf((char *)buffer,"%lld%n",target,&nread);
               if(nread != ntoken || items != 1) goto error_type ;
               }
            else
               goto error_syntax;
            break;
         case(WB_FORTRAN_REAL):
            if(elementsize==4){              /* 4 byte real */
               float *target = (float *)dataptr;
               items = sscanf((char *)buffer,"%E%n",target,&nread);
               if(nread != ntoken || items != 1) goto error_syntax ;
               }
            else if(elementsize==8) {        /* 8 byte real */
               double *target = (double *)dataptr;
               items = sscanf((char *)buffer,"%lE%n",target,&nread);
               if(nread != ntoken || items != 1) goto error_syntax ;
               }
            else
               goto error_syntax;
            break;
         case(WB_FORTRAN_BOOL):
            if(elementsize==4){              /* FORTRAN logical */
               if( strncmp(".FALSE.",(char *)buffer,7) == 0 || strncmp(".F.",(char *)buffer,3) == 0 ) *dataptr = 0;
               else if( strncmp(".TRUE.",(char *)buffer,6) == 0 || strncmp(".T.",(char *)buffer,3) == 0 ) *dataptr = 1;
               else goto error_type;
               }
            else
               goto error_syntax;
            break;
         case(WB_FORTRAN_CHAR): /* copy string in buffer into target, omit first and last character, the quotes */
            nread = c_fortran_string_copy(buffer+1, dataptr, strlen((char *)buffer)-2, elementsize, ' ');
            if(nread < 0) goto error_string;
            break;
         default:
            goto error_syntax;
         }
      dataptr += elementsize;

      ntoken=wb_get_token(separator,infile,sizeof(separator)-1,0);   /* separator must be either comma or newline */
      if( ntoken != 1 ) goto error_syntax;
      if( separator[0] == ',' ) continue ;
      if( separator[0] == '\n' ) break ;
      else                      goto error_syntax;
      }
   if( separator[0] != '\n' ) goto error_toomany;
   line->m.flags.initialized = 1 ;  /* mark entry as initialized */
   line->m.flags.badval = 0 ;
   return(status);
error_toomany:
   if(message_level<=WB_MSG_ERROR) {
      fprintf(stderr,"attempting to assign too many values to %s (type %s)\n",name,datatypes[elementtype]);
      goto error_syntax;
      }
error_string:
   if(message_level<=WB_MSG_ERROR) {
      fprintf(stderr,"string %s is too long to be assigned to %s (length=%d>%d)\n",buffer,name,strlen((char *)buffer)-2,elementsize);
      goto error_syntax;
      }
error_type:
   if(message_level<=WB_MSG_ERROR) {
      fprintf(stderr,"type mismatch error while assigning value %s to %s (type %s)\n",buffer,name,datatypes[elementtype]);
      goto error_syntax;
      }
error_syntax:
   wb_read_error() ;
   wb_flush_line(infile);
   WB_ERR_EXIT(WB_MSG_ERROR,WB_ERR_SYNTAX);
}

/* read a dictionary or user directive file */

int c_wb_read(WhiteBoard *WB, char *filename, char *package, char *section, int Options, int Lfilename, int Lpackage, int Lsection){
   DEFTBL mytable[MISC_BUFSZ];
   char localbuffer[MISC_BUFSZ];
   char localfname[MISC_BUFSZ];
   unsigned char Token[MISC_BUFSZ];
   char Package[WB_MAXNAMELENGTH];
   char Section[WB_MAXNAMELENGTH];
   int i,status;
   FILE *infile;
   int ntoken;
   int temp;
   int version=999,items,nread;
   int errors=0;

   extra_error_message = " invalid whiteboard instance";
   if(WB == DummyWhiteboardPtr) WB_ERR_EXIT(WB_MSG_ERROR,WB_ERR_NOTFOUND);
   if(WB == NULL) WB=BaseWhiteboardPtr;

   TRIM(filename,Lfilename)
   TRIM(package,Lpackage)
   TRIM(section,Lsection)
   linebuffer=localbuffer;
   definition_table=mytable; definition_table_entries=0; max_definition_table_entries=MISC_BUFSZ;  /* initialize definition table control */
   for( i=0 ; i<MISC_BUFSZ ; i++) { definition_table[i].line = NULL ; definition_table[i].assigned = 0 ; definition_table[i].defined = 0 ; }

   for( i=0 ; i<Lfilename && i<MISC_BUFSZ-1 && filename[i]!=0 ; i++) localfname[i] = filename[i];
   localfname[i] = 0;  /* filename collected */
   if(message_level<=WB_MSG_DEBUG) fprintf(stderr,"localfname='%s'\n",localfname);

   for( i=0 ; i<Lpackage && i<WB_MAXNAMELENGTH-1 && package[i]!=0 ; i++) Package[i] = toupper(package[i]);
   Package[i] = 0;     /* package name collected */
   if(message_level<=WB_MSG_DEBUG) fprintf(stderr,"Package='%s'\n",Package);

   for( i=0 ; i<Lsection && i<WB_MAXNAMELENGTH-1 && section[i]!=0 ; i++) Section[i] = toupper(section[i]);
   Section[i] = 0;     /* section_name  collected */
   if(message_level<=WB_MSG_DEBUG) fprintf(stderr,"Section='%s'\n",Section);

   extra_error_message = localfname;  /* add directive file name to error message */
   current_char = NULL;               /* make sure that no previous input is left in buffers */
   infile=fopen(localfname,"r");      /* try to open file */
   if(infile==NULL) WB_ERR_EXIT(WB_MSG_ERROR,WB_ERR_READ);  /* OOPS, cannot open/read file */

   extra_error_message = NULL;

   while(1){   /* loop until desired section "@xxxx" is found or end of directive file */
      temp = wb_get_nonblank(infile);
      while( temp!= '@' && temp!= EOF ) {   /* look for @ as first nonblank character in line (or End of file) */
         temp = wb_get_line(infile);
         temp = wb_get_nonblank(infile);
         }
      if(temp==EOF) {                    /* end of file or error , requested section not found */
         fclose(infile);
         extra_error_message = Section;  /* add section name to error message */
         WB_ERR_EXIT(WB_MSG_ERROR,WB_ERR_NOTFOUND);
         }
      ntoken = wb_get_token(Token,infile,MISC_BUFSZ-1,1);
      if( strncmp((char *)Token,Section,strlen(Section))==0 ) break;   /* exit loop if requested section is found */
      }
   if( message_level<=WB_MSG_INFO )fprintf(stderr,"INFO: directive section %s found\n",Section);

   /* we have found our section, let's process it */
   temp = wb_get_line(infile);     /* get rid of rest of @section line */
   default_dict_option = 0;                /* default options = none unless there is a directive */
   ntoken = wb_get_token(Token,infile,MISC_BUFSZ-1,0);
   while( strncmp((char *)Token,"@",1) && ntoken!=WB_ERROR ) {  /* loop until end of section (beginning of next section or EOF) */
      if( strncmp((char *)Token,"OPTIONS",7)==0 ) {   /* default options directive */
         default_dict_option = wb_options(infile,'(',')',dict_options);
         }
      else if( strncmp((char *)Token,"MESSAGES",8)==0 ) {   /* verbosity control directive */
         message_level = wb_options(infile,'(',')',verb_options);
         }
      else if( strncmp((char *)Token,"DEFINE",6)==0 && Options!=WB_FORBID_DEFINE ) {  /* define directive (if allowed) */
         status = wb_define(WB,infile,Package);
         if(status<0) errors++;
         }
      else if( isalpha(Token[0]) ){   /* must be key= */
         status = wb_key(WB,infile,Token,Package,Options);
         if(status<0) errors++;
         }
      else if( Token[0] == '\n' ) {  /* ignore newlines */
         temp = '\n';
         }
      else{
         if( message_level<=WB_MSG_ERROR )
            fprintf(stderr,"Unexpected token found at beginning of directive: '%s'\n",Token);
         errors++;
         wb_flush_line(infile);
         }
      ntoken = wb_get_token(Token,infile,MISC_BUFSZ-1,0);   /* get next token */
      }
   fclose(infile);    /* we are done, close file and return success */
   if( message_level<=WB_MSG_INFO || ( message_level<=WB_MSG_ERROR && errors>0 ) )
      fprintf(stderr,"INFO: %d error(s) detected\n",errors);
   return( errors==0 ? WB_OK : WB_ERROR);
}

/* read a dictionary or user directive file (FORTRAN version) */
wordint f77_name(f_wb_read)(WhiteBoard **WB, char *package, char *filename, char *section, wordint *options, F2Cl lpackage, F2Cl lfilename,  F2Cl lsection){
   int Lfilename=lfilename;
   int Lpackage=lpackage;
   int Lsection=lsection;
   int Options=*options;

/*
   if(*WB != NULL) fprintf(stderr,"INFO: files can only be read into MAIN WhiteBoard, option IGNORED\n");
*/
   return(c_wb_read(*WB, filename, package, section, Options, Lfilename, Lpackage, Lsection));
}

void f77_name(c_wb_test)() {
   int status,myint,myint2;
   long long myll,myll2;
   float myreal,myreal2;
   double mydouble,mydouble2;
   int integer_array[]={1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25};
   int integer_array2[25];
   char *string1="this is string 1";
   char string2[32];
   WhiteBoard *WB=NULL;

   message_level = WB_MSG_INFO;
   status=WB_PUT_C(WB,"string1",string1,strlen(string1),0+WB_REWRITE_UNTIL);
   printf("status=%d\n",status);
   myint = 12; myll=1212;
   myreal=12.12; mydouble=2424.68;
   status=WB_PUT_I4(WB,"valeur1",&myint,WB_CREATE_ONLY);
   printf("status=%d\n",status);
   status=WB_PUT_I8(WB,"ll1",&myll,0+WB_READ_ONLY_ON_RESTART);
   printf("status=%d\n",status);
   status=WB_PUT_R4(WB,"reel1",&myreal,0+WB_READ_ONLY_ON_RESTART);
   printf("status=%d\n",status);
   status=WB_PUT_R8(WB,"dble1",&mydouble,0+WB_REWRITE_UNTIL);
   printf("status=%d\n",status);
   myint = -134;myll=-134431;
   myreal=-13.45; mydouble=-12345.6789;
   status=WB_PUT_I4(WB,"valeur2",&myint,0+WB_REWRITE_UNTIL);
   printf("status=%d\n",status);
   status=WB_PUT_I8(WB,"ll2",&myll,0);
   printf("status=%d\n",status);
   status=WB_PUT_R4(WB,"reel2",&myreal,0);
   printf("status=%d\n",status);
   status=WB_PUT_R8(WB,"dble2",&mydouble,0);
   printf("status=%d\n",status);
printf("========\n");
   status=c_wb_check(WB,(unsigned char *)"", -1, 0, 1, NULL,NULL);
   printf("c_wb_check printed %d entries\n",status);
printf("========\n");
   printf("before c_wb_lock\n");
   c_wb_lock(WB,(unsigned char *)"D",1);

   c_wb_checkpoint();

   BaseWhiteboardPtr->PageChain=NULL;
   c_wb_reload();

   status=WB_GET_C(WB,"string1",string2,32);
   string2[31]=0;
   printf("status=%d, string2='%s'\n",status,string2);
   myint2 = -1; myll2=-1;
   myreal2=-1.0; mydouble2=-1.1111111111;
   status=WB_GET_I4(WB,"valeur1",&myint2);
   printf("status=%d, myint2=%d\n",status,myint2);
   status=WB_GET_I8(WB,"ll1",&myll2);
   printf("status=%d, myll2=%lld\n",status,myll2);
   status=WB_GET_R4(WB,"reel1",&myreal2);
   printf("status=%d, myreal2=%lf\n",status,myreal2);
   status=WB_GET_R8(WB,"dble1",&mydouble2);
   printf("status=%d, mydouble2=%f\n",status,mydouble2);
   myint2 = -1; myll2=-1;
   myreal2=-1.0; mydouble2=-1.1111111111;
   status=WB_GET_I4(WB,"valeur2",&myint2);
   printf("status=%d, myint2=%d\n",status,myint2);
   status=WB_GET_I8(WB,"ll2",&myll2);
   printf("status=%d, myll2=%lld\n",status,myll2);
   status=WB_GET_R4(WB,"reel2",&myreal2);
   printf("status=%d, myreal2=%f\n",status,myreal2);
   status=WB_GET_R8(WB,"dble2",&mydouble2);
   printf("status=%d, mydouble2=%lf\n",status,mydouble2);
printf("========\n");
   status=WB_PUT_I4V(WB,"intarray1",integer_array,25,0);
   printf("status=%d\n",status);
printf("========\n");
   status=WB_GET_I4V(WB,"intarray1",integer_array2,5);
   printf("status=%d, integer_array2[4]=%d\n",status,integer_array2[4]);

   c_wb_checkpoint();

   BaseWhiteboardPtr->PageChain=NULL;
   c_wb_reload();

printf("========\n");
   status=WB_PUT_I4V(WB,"intarray1",integer_array,31,0);
   printf("status=%d\n",status);
printf("========\n");
   status=WB_GET_I4V(WB,"intarray1",integer_array2,5);
   printf("status=%d, integer_array2[4]=%d\n",status,integer_array2[4]);
   status=c_wb_check(WB,(unsigned char *)"", -1, 0, 1, NULL,NULL);
   printf("c_wb_check printed %d entries\n",status);
printf("========\n");
   status=WB_PUT_I4V(WB,"intarray1",integer_array,10,0);
   printf("status=%d\n",status);
printf("========\n");
   status=WB_GET_I4V(WB,"intarray1",integer_array2,5);
   printf("status=%d, integer_array2[4]=%d\n",status,integer_array2[4]);
printf("========\n");
   status=WB_GET_I4V(WB,"intarray1",integer_array2,15);
   printf("status=%d, integer_array2[14]=%d\n",status,integer_array2[14]);
printf("==== Locking Test ====\n");
   status=c_wb_check(WB,(unsigned char *)"V", -1, 1, 1, Action1,NULL);
   printf("c_wb_check printed %d entries\n",status);
   c_wb_lock(WB,"V",1);
   status=c_wb_check(WB,(unsigned char *)"V", -1, 1, 1, Action1,NULL);
   printf("c_wb_check printed %d entries\n",status);
printf("==== END of Locking Test ====\n");
}
