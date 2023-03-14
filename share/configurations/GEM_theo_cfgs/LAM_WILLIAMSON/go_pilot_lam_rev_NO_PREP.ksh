set -x

name_exp=$1

run1=$2
run2=$4

name1=$3
name2=$5

exp1=${name_exp}/${name1}
exp2=${name_exp}/${name2}

if [ -z "${run1}" -a -z "${run2}"  ] ; then
echo "PILOT run1=0/1 AND LAM run2=0/1"
exit
fi

######################################
#CAUTION: CHANGE both REP_pilot and REP_lam
######################################

#----------------------------
#run GEM PILOT
#----------------------------
if test "${run1}" = "1"; then

REP_pilot=/users/dor/armn/mta/raid/${TRUE_HOST}/resultats_nlm/williamson_5.1.a5/${name_exp}/${name1}
mkdir -p $REP_pilot

/bin/rm -rf RUN*/*
/bin/rm -rf PREP/* 

runprep -dircfg ${exp1} 
###runmod  -dircfg ${exp1} -ptopo 6x6x1 -inorder > std1_mod_${name1}
runmod  -dircfg ${exp1} -ptopo 3x3x1 -inorder > std1_mod_${name1}

if [ "$_status" == "ABORT" ] ; then
   printf "\n\n   PROBLEM with Um_runmod -ABORT \n\n"
   exit 1
fi

/bin/rm -rf $REP_pilot/*
cp RUNMOD/output/cfg_0000/laststep_000000*/0*/* $REP_pilot
mkdir $REP_pilot/casc_PILOT
scp $REP_pilot/casc_* $REP_pilot/casc_PILOT/. 
\rm -rf $REP_pilot/casc_*00* 
cp std1_mod_${name1} $REP_pilot/.

cat > ${exp2}/cfg_0000/configexp.cfg <<EOF
GEM_inrep=$REP_pilot/casc_PILOT
EOF

#----------------------------
fi
#----------------------------

#----------------------------
#run GEM LAM 
#----------------------------
if test "${run2}" = "1"; then

REP_lam=/users/dor/armn/mta/raid/${TRUE_HOST}/resultats_nlm/williamson_5.1.a5/${name_exp}/${name2}
mkdir -p $REP_lam

#################################################################################################################################
### COMMENTER LES 3 LINES SUIVANTES si on a deja tourne RUNPREP pour un LAM et que le transfert des fichiers de pilotage est long 
#################################################################################################################################
/bin/rm -rf RUN*/*
###/bin/rm -rf PREP/* 
###runprep -dircfg ${exp2} 

runmod  -dircfg ${exp2} -ptopo 4x2x1 > std2_mod_${name2}
###runmod  -dircfg ${exp2} -ptopo 1x1x1 > std2_mod_${name2}

mkdir -p $REP_lam/casc_${name2}
cp RUNMOD/output/cfg_0000/laststep_000000*/0*/* $REP_lam/casc_${name2}
cp std2_mod_${name2} $REP_lam/.

###xrec -imflds $REP/casc/dm* $REP/dm*

#----------------------------
fi
#----------------------------
