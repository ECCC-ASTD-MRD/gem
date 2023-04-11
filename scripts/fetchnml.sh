#!/bin/bash

# ie: fetchnml.sh grd_nj grid gem_settings.nml
# pas encore fini
string=`echo $1 | tr "[a-z]" "[A-Z]"`
namelist=`echo $2 | tr "[a-z]" "[A-Z]"`

cat > $TMPDIR/tmp$$ <<EOF
function ltrim(s) { sub(/^[ \t\r\n]+/, "", s); return s }
function rtrim(s) { sub(/[ \t\r\n]+$/, "", s); return s }
function trim(s)  { return rtrim(ltrim(s)); }
BEGIN{
key=0 
nml=0
}
/&${namelist}$/ {nml=1; nmlfound=1}
/\// {nml=0}
{
if (key == 1 && nml == 1) {
    if (index(trim(\$0),"'") > 1) { #quote found in line not on col 1
       printf("\n")
       key =2
    }
    else if (index(trim(\$0),"'") == 1 ) { #quote found on col 1
       inquote = 0
       nq = split(trim(\$0),qfield,"'")
       for (k=1; k<=nq; k++){
            if (inquote==1) {
               if (nval==1){
               printf(",\n%s",trim(qfield[k]))
               nval=0
               }
               else{printf("%s",trim(qfield[k]))}
               }
            else {
                 if (index(qfield[k],"=") !=0) {key=2;printf("\n")}
                 else { printf("%s",trim(qfield[k])) }
                 }
           if (inquote==0) {
               inquote=1
               }
           else {
              inquote=0
              }
           }
         if (key==1) {printf("\n")}
         }
    else {
    if (length(\$0) !=0) {
         if (index(\$0,"=") ==0){
             printf("%s\n",trim(\$0))
             }
         else {
            nf = split(\$0,field,"=")
            val=trim(substr(field[1],1,(index(field[1],",")-1)))
            printf("%s\n",trim(val))
            key=2
              }
          }
      }
}
if (nml == 1 && key==0) {
  if (\$0~/${string}/) {
    if (index(\$0,"'")) { #quotes in line found
       nq = split(trim(\$0),qfield,"'")
       inquote=0
       nval=0
       for (i=1; i<=nq; i++){
           if (key ==1) {
               if (inquote==1) {
               printf("%s",trim(qfield[i]))
               nval=1
               }
               else {
                 if (index(qfield[i],"=") !=0) {key=2;printf("\n")}
                 if (i != nq) { printf("%s",trim(qfield[i])) }
               }
           }
           if (inquote==0 && key == 0 ) {
               key=findkey("${string}",qfield[i])
               }
           if (inquote==0) {
               inquote=1
               }
           else {
              inquote=0
              }
           }
           }
    else { #no quotes in line found
        key=findkey("${string}",\$0)
         }
        }
    }
}
#END {
#if (nmlfound !=1) { printf("WARNING: Namelist ${namelist}  not found\n")}
#if (nmlfound ==1 && mykey ==0) { printf("WARNING: Key ${string} in namelist ${namelist} not found\n")}
#}

func findkey(name,mystring)
{
# look for namelist key
    nv = split(mystring,field,"=")
    mykey=0
    for (j=1; j<=nv; j++){
         ii=index(field[j],",")
         if (trim(substr(field[j],ii+1,length(field[j]))) == trim(name)) {
            mykey=1
           if (index(field[j+1],",") !=0){
             val=trim(substr(field[j+1],1,index(field[j+1],",")-1))
             mykey=2
             printf("%s\n",trim(val))
             }
         else if (length(trim(field[j+1])) > 0) {
             printf("%s\n",trim(field[j+1]))
             mykey=2
             }
            }
            }
         return(mykey)
}
EOF
cat $3 | tr "\"" "\'" | tr "[a-z]" "[A-Z]" | awk -f $TMPDIR/tmp$$ 
/bin/rm -f $TMPDIR/tmp$$
