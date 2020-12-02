set -x

REP=${HOME}/raid/${TRUE_HOST}/acid_test_results
/bin/rm -rf $REP std1_mod std2_mod ; mkdir -p $REP/casc

linkit
/bin/rm -rf RUNMOD/* PREP/*
runprep -dircfg GEM_cfgs_acid1
runmod  -dircfg GEM_cfgs_acid1 -ptopo 4x4 -inorder > std1_mod
cp RUNMOD/output/cfg_0000/laststep_000000*/0*/* $REP
mv $REP/casc_* $REP/casc

/bin/rm -rf RUNMOD/* PREP/*
runprep -dircfg GEM_cfgs_acid2
runmod  -dircfg GEM_cfgs_acid2 -ptopo 4x4 -inorder > std2_mod

cp RUNMOD/output/cfg_0000/laststep_000000*/0*/* $REP/casc

xrec -imflds $REP/casc/dm* $REP/dm*
