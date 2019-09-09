TMPL90DECKS=
FTN90DECKS=
PTNDECKS=
CDKDECKS=
CDECKS= \
	vgrid.c  
HDECKS= \
	armnlib.h  vgrid.h  
FHDECKS=
FTNDECKS=
HCDECKS= \
	BODY_C_compute_heights_0001.hc  BODY_C_compute_heights_21001.hc  BODY_C_compute_heights_4001.hc  BODY_C_compute_pressure_1001_1002.hc  BODY_C_compute_pressure_1003_5001.hc   \
	BODY_C_compute_pressure_2001.hc  BODY_C_compute_pressure_5002_5003_5004_5005.hc  BODY_C_compute_pressure_5100.hc  BODY_Cvgd_diag_withref_2ref.hc  vgrid_version.hc  
CDK90DECKS=
INCDECKS=
HFDECKS= \
	vgrid_descriptors.hf  
F90DECKS=
PTN90DECKS=
FDECKS=
CHDECKS=
OBJECTS= \
	vgrid.o  vgrid_descriptors.o  vgrid_utils.o  
vgrid.o:	vgrid.c \
	BODY_C_compute_pressure_2001.hc  BODY_C_compute_pressure_5100.hc  BODY_C_compute_heights_4001.hc  BODY_Cvgd_diag_withref_2ref.hc  armnlib.h   \
	BODY_C_compute_pressure_5002_5003_5004_5005.hc  BODY_C_compute_heights_0001.hc  BODY_C_compute_pressure_1001_1002.hc  BODY_C_compute_pressure_1003_5001.hc  BODY_C_compute_heights_21001.hc   \
	vgrid.h  
vgrid_descriptors.o:	vgrid_descriptors.F90 \
	vgrid_utils.o  vgrid_descriptors.hf  
vgrid_utils.o:	vgrid_utils.F90 \
	vgrid_descriptors.hf  
