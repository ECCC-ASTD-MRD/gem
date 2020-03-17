#!/bin/bash

# Script used to define your own grid and visualise it with SPI First, copy
# one of the gem_settings.nml files located in the different configurations
# available, edit it, and then run the command grille.sh -spi
# If you want geophysical fields and historical meteorological data for the
# region you defined in that new grid, contact us.  

export PATH=${PWD}:$PATH
eval `cclargs_lite0 -D $0 \
  -xrec     "0"     "1"     "[Visualize grid with xrec                ]"\
  -spi      "0"     "1"     "[Visualize grid with spi                 ]"\
  ++ $*`

set -ex

printf "\n######################### \n"
printf " GRILLE utility for GEM \n"
printf "######################### \n\n"

GEM_ON_A_STICK_WORK_DIR=`(cd ..; pwd)`
GEMGRID=${GEM_ON_A_STICK_WORK_DIR}/bin/maingemgrid
ATM_MODEL_DFILES=${GEM_ON_A_STICK_WORK_DIR}/../../gem_dbase

if [ ! -s gem_settings.nml ] ; then
  echo "Must provide file gem_settings.nml containing a valid &grid namelist"
  exit
fi

export TMPDIR=${TMPDIR-/tmp/$USER}
mkdir -p ${TMPDIR}

ROOT_WORK=${TMPDIR}/grid$$
WORKDIR=${ROOT_WORK}/cfg_0001
/bin/rm -rf ${ROOT_WORK} ; mkdir -p ${WORKDIR}/CACHE

OUTPUT=${PWD}/grid$$
mkdir -p ${OUTPUT}
cp gem_settings.nml ${WORKDIR}/model_settings.nml

ngrids=1
GRDTYP=$(./Um_fetchnml2.sh Grd_typ_s grid gem_settings.nml)
if [ "$GRDTYP" == "GY" ] ; then 
   export GEM_YINYANG=YES
   ngrids=2
   mkdir -p ${WORKDIR}/YIN/000-000 ${WORKDIR}/YAN/000-000
   ln -s ${ATM_MODEL_DFILES}/datafiles/constants/thermoconsts ${WORKDIR}/YIN/000-000/constantes
   ln -s ${ATM_MODEL_DFILES}/datafiles/constants/thermoconsts ${WORKDIR}/YAN/000-000/constantes
else
   mkdir -p ${WORKDIR}/000-000
   ln -s ${ATM_MODEL_DFILES}/datafiles/constants/thermoconsts ${WORKDIR}/000-000/constantes
fi

export TASK_INPUT=${ROOT_WORK}
export TASK_WORK=${ROOT_WORK}
export DOMAIN_start=1
export DOMAIN_end=1
export DOMAIN_total=1
export GEM_NDOMAINS=1:1

# Produce Grid Pos Rec
/bin/rm -f tape1 tape1_core tape1_free tape2 tape2_core tape2_free 2> /dev/null

mpirun -np 1 ${GEMGRID}

find ${WORKDIR} -type f -name "*tape*" -exec mv {} ${OUTPUT}/ \;
find ${WORKDIR} -type f -name "eigenv_v1_*" -exec mv {} ${OUTPUT}/ \;

/bin/rm -rf $ROOT_WORK

# View the grid(s)
function launch_spi {
   if ! which SPI; then
      echo "==================================================="
      echo "   Error!"
      echo "   Please load a version of SPI."
      echo "   Search the wiki for SPI to learn how to load it."
      echo "   Aborting"
      echo "==================================================="
      rm -rf $TMPDIR
      exit 1
   fi
   SPI -field $*
}

viewer=""
if [ $xrec -eq 1 ]; then
   viewer="xrec -imflds"
elif [ $spi -eq 1 ]; then
   viewer="launch_spi"
fi

if [ -n "$viewer" ] ; then
   cat > ${OUTPUT}/directives.pgsm << EOF
 sortie (std,1000,a)
 compac=-12
 etiksrt='GRID'
 LAGRILLE
*
 heure(0)
 outlalo(-1,-1,-1)
*
 setintx(cubique)
 champ('ME')
end
EOF
   tmp1=''
   if [ -s ${OUTPUT}/tape1 ] ; then
      ./r.fstliste.sh -izfst ${OUTPUT}/tape1 -nomvar ">>" | sed 's/\:/ /g' > ${OUTPUT}/liste
      ip1=`cat ${OUTPUT}/liste | awk '{print $3}'`
      ip2=`cat ${OUTPUT}/liste | awk '{print $4}'`
      ip3=`cat ${OUTPUT}/liste | awk '{print $5}'`
      ETK=`cat ${OUTPUT}/liste | awk '{print $9}'`
      cat ${OUTPUT}/directives.pgsm | sed "s/LAGRILLE/grille (tape2,$ip1,$ip2,$ip3)/" | sed "s/GRID/$ETK/" > ${OUTPUT}/p1.dir
      tmp1=${OUTPUT}/tmp1$$
      cp ${OUTPUT}/tape1 ${tmp1}
      pgsm -iment $ATM_MODEL_DFILES/bcmk/geophy/Gem_geophy.fst -ozsrt ${tmp1} -i ${OUTPUT}/p1.dir
   fi
   tmp2=''
   if [ -s ${OUTPUT}/tape2 ] ; then
      ./r.fstliste.sh -izfst  ${OUTPUT}/tape2 -nomvar ">>" | sed 's/\:/ /g' > ${OUTPUT}/liste
      ip1=`cat ${OUTPUT}/liste | awk '{print $3}'`
      ip2=`cat ${OUTPUT}/liste | awk '{print $4}'`
      ip3=`cat ${OUTPUT}/liste | awk '{print $5}'`
      ETK=`cat ${OUTPUT}/liste | awk '{print $9}'`
      cat ${OUTPUT}/directives.pgsm | sed "s/LAGRILLE/grille (tape2,$ip1,$ip2,$ip3)/" | sed "s/GRID/$ETK/" > ${OUTPUT}/p2.dir
      tmp2=${OUTPUT}/tmp2$$
      cp ${OUTPUT}/tape2 ${tmp2}
      pgsm -iment $ATM_MODEL_DFILES/bcmk/geophy/Gem_geophy.fst -ozsrt ${tmp2} -i ${OUTPUT}/p2.dir
   fi
   $viewer ${tmp1} ${tmp2}
   /bin/rm -f ${OUTPUT}/directives.pgsm ${OUTPUT}/p1.dir ${OUTPUT}/p2.dir ${OUTPUT}/liste ${tmp1} ${tmp2}
fi

rm -rf $TMPDIR
