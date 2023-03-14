#!/bin/ksh -ex
# run this scripts in [dev_dir]/GEM_cfgs_acid directory

if [ -z "${ATM_MODEL_BNDL}" ] ; then
  echo MUST set GEM env first -- ABORT --
  exit 1
fi

ici=$(basename $PWD)
REP=${HOME}/raid/\${TRUE_HOST}/acid_test_results
/bin/rm -rf ../GEM_cfgs_acid1 ../GEM_cfgs_acid2 ../go.ksh
mkdir -p ../GEM_cfgs_acid1 ../GEM_cfgs_acid2

cp -r cfg_0000 ../GEM_cfgs_acid1
cp -r cfg_0001 ../GEM_cfgs_acid2/cfg_0000 

cp configexp.cfg ../GEM_cfgs_acid1/cfg_0000
cat >> ../GEM_cfgs_acid1/cfg_0000/configexp.cfg <<EOF
GEM_inrep=$ATM_MODEL_DFILES/bcmk
GEM_anal=\$GEM_inrep/2009042700_000
GEM_xfer=${REP}
EOF

cp configexp.cfg ../GEM_cfgs_acid2/cfg_0000
cat >> ../GEM_cfgs_acid2/cfg_0000/configexp.cfg <<EOF
GEM_inrep=${REP}/casc
GEM_anal=$ATM_MODEL_DFILES/bcmk/2009042700_000
GEM_xfer=${REP}/casc
EOF

cd ../
rm -f go.ksh
ln -s ${ici}/go.ksh

