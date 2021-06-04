#!/bin/bash

string=`echo $1 | tr "[a-z]" "[A-Z]"`
nml=`   echo $2 | tr "[a-z]" "[A-Z]"`

cat > tmp$$ <<EOF
BEGIN{pr=0}
{ if (pr==1) print }
/&${nml}$/ {pr=1}
/\// {pr=0}
EOF
cat $3 | tr "[a-z]" "[A-Z]" | awk -f tmp$$ |\
  grep -i $string | tr "[a-z]" "[A-Z]" | sed "s-\(.*\)\(${string} *=.*\)-\2-" \
                  | sed "s/,/ /g" | sed "s/=/ /g" | awk '{print $2}'          \
                  | sed "s/'//g" | sed 's/"//g'
/bin/rm -f tmp$$







