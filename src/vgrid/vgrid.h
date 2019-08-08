#include <stdint.h>
#ifndef VGRID_H
#define VGRID_H

#define VGD_OK       0
#define VGD_ERROR    -1
#define VGD_MISSING  -9999.
#define VGD_LEN_NAME 5
#define VGD_LEN_RFLD 5
#define VGD_LEN_RFLS 5
#define VGD_LEN_ETIK 13
#define VGD_LEN_TYPVAR 3
#define VGD_LEN_GRTYP  2
#define VGD_NO_REF_NOMVAR "    "
extern float VGD_STDA76_SFC_T;
extern float VGD_STDA76_SFC_P;


typedef struct VGD_TFSTD {
   int   dateo;                 // date d'origine du champs
   int   deet;                  // duree d'un pas de temps
   int   npas;                  // pas de temps
   int   nbits;                 // nombre de bits du champs
   int   datyp;                 // type de donnees
   int   ip1,ip2,ip3;           // specificateur du champs
   int   ig1,ig2,ig3,ig4;       // descripteur de grille
   char  typvar[VGD_LEN_TYPVAR]; // type de variable
   char  nomvar[VGD_LEN_NAME]; // nom de la variable
   char  etiket[VGD_LEN_ETIK]; // etiquette du champs
   char  grtyp[VGD_LEN_GRTYP];   // type de grilles
   char  fstd_initialized;      // if the fstd struct is initialized
} VGD_TFSTD;

typedef struct VGD_TFSTD_ext {
   int   dateo;               // date d'origine du champs
   int   datev;               // date de validitee du champs
   int   deet;                // duree d'un pas de temps
   int   npas;                // pas de temps
   int   nbits;               // nombre de bits du champs
   int   datyp;               // type de donnees
   int   ip1,ip2,ip3;         // specificateur du champs
   int   ni,nj,nk;            // dimensions
   int   ig1,ig2,ig3,ig4;     // descripteur de grille
   int   swa;
   int   lng;
   int   dltf;
   int   ubc;
   int   extra1,extra2,extra3;
   char  typvar[VGD_LEN_TYPVAR]; // type de variable
   char  nomvar[VGD_LEN_NAME]; // nom de la variable
   char  etiket[VGD_LEN_ETIK]; // etiquette du champs
   char  grtyp[VGD_LEN_GRTYP];   // type de grilles
} VGD_TFSTD_ext;

typedef struct vgrid_descriptor {
  VGD_TFSTD rec;          // RPN standard file header
  double   ptop_8;        // Top level pressure (Pa)
  double   pref_8;        // Reference pressure (Pa)
  double   *table;        // Complete grid descriptor record
  int      table_ni;      //    ni size of table
  int      table_nj;      //    nj size of table
  int      table_nk;      //    nk size of table
  double   *a_m_8;        // A-coefficients for momentum levels  
  double   *b_m_8;        // B-coefficients for momentum levels
  double   *c_m_8;        // C-coefficients for momentum levels
  double   *a_t_8;        // A-coefficients for thermodynamic levels
  double   *b_t_8;        // B-coefficients for thermodynamic levels
  double   *c_t_8;        // C-coefficients for thermodynamic levels
  double   *a_w_8;        // A-coefficients for vertical-velocity levels
  double   *b_w_8;        // B-coefficients for vertical-velocity levels
  double   *c_w_8;        // C-coefficients for vertical-velocity levels
  int      *ip1_m;        // ip1 values for momentum levels
  int      *ip1_t;        // ip1 values for momentum levels
  int      *ip1_w;        // ip1 values for momentum levels
  int      nl_m;          // Number of momentum      level
  int      nl_t;          // Number ot thermodynamic level
  int      nl_w;          // Number ot vertical volocity level
  float    dhm;           // Diag level Height (m) for Momentum variables UU,VV
  float    dht;           // Diag level Height (t) for Thermo variables TT,HU, etc
  float    dhw;           // Diag level Height (m) for vertical-velocity variables WT1 ZZ
  char*    ref_name;      // Reference field name
  char*    ref_namel;     // Reference field name
  float    rcoef1;        // Rectification coefficient
  float    rcoef2;        // Rectification coefficient
  float    rcoef3;        // Rectification coefficient
  float    rcoef4;        // Rectification coefficient
  int      nk;            // Number of momentum levels
  int      ip1;           // ip1 value given to the 3D descriptor
  int      ip2;           // ip2 value given to the 3D descriptor
  int      unit;          // file unit associated with this 3D descriptor
  int      vcode;         // Vertical coordinate code
  int      kind;          // Vertical coordinate code
  int      version;       // Vertical coordinate code
  char     match_ipig;    // do ip/ig matching for records
  char     valid;         // Validity of structure
} vgrid_descriptor;

int Cvgd_is_valid(vgrid_descriptor *self, char *valid_table_name);
void Cvgd_table_shape(vgrid_descriptor *self, int **tshape);
int Cvgd_print_desc(vgrid_descriptor *self, int sout, int convip);
int Cvgd_print_vcode_description(int vcode);
int Cvgd_levels_8(vgrid_descriptor *self, int ni, int nj, int nk, int *ip1_list, double *levels_8, double *sfc_field_8, int in_log);
int Cvgd_levels(vgrid_descriptor *self, int ni, int nj, int nk, int *ip1_list, float *levels, float *sfc_field, int in_log);
int Cvgd_levels_2ref_8(vgrid_descriptor *self, int ni, int nj, int nk, int *ip1_list, double *levels_8, double *sfc_field_8, double *sfc_field_ls_8, int in_log);
int Cvgd_levels_2ref(vgrid_descriptor *self, int ni, int nj, int nk, int *ip1_list, float *levels, float *sfc_field, float *sfc_field_ls, int in_log);
int Cvgd_diag_withref_8(vgrid_descriptor *self, int ni, int nj, int nk, int *ip1_list, double *levels_8, double *sfc_field_8, int in_log, int dpidpis);
int Cvgd_diag_withref(vgrid_descriptor *self, int ni, int nj, int nk, int *ip1_list, float *levels, float *sfc_field, int in_log, int dpidpis);
int Cvgd_diag_withref_2ref_8(vgrid_descriptor *self, int ni, int nj, int nk, int *ip1_list, double *levels_8, double *sfc_field_8, double *sfc_field_ls_8, int in_log, int dpidpis);
int Cvgd_diag_withref_2ref(vgrid_descriptor *self, int ni, int nj, int nk,
			   int *ip1_list, float *levels, float *sfc_field,
			   float *sfc_field_ls, int in_log, int dpidpis);
int Cvgd_set_vcode_i(vgrid_descriptor *VGrid,int Kind,int Version);
int Cvgd_set_vcode(vgrid_descriptor *VGrid);
int Cvgd_new_build_vert(vgrid_descriptor **self, int kind, int version, int nk, int ip1, int ip2, double *ptop_8, double *pref_8, float *rcoef1, float *rcoef2,
			double *a_m_8, double *b_m_8, double *a_t_8, double *b_t_8, int *ip1_m, int *ip1_t, int nl_m, int nl_t);
int Cvgd_new_build_vert_1001(vgrid_descriptor **self, int ip1, int ip2, 
			     double *a_m_8, double *b_m_8, int *ip1_m, int nk);
int Cvgd_new_build_vert_1002(vgrid_descriptor **self, int ip1, int ip2, double ptop_8,
			     double *a_m_8, double *b_m_8, int *ip1_m, int nk);
int Cvgd_new_build_vert_2001(vgrid_descriptor **self, int ip1, int ip2, 
			     double *a_m_8, double *b_m_8, int *ip1_m, int nk);
int Cvgd_new_build_vert_4001(vgrid_descriptor **self, int ip1, int ip2, 
			     double *a_m_8, double *b_m_8, int *ip1_m, int nk);
int Cvgd_new_build_vert_5001(vgrid_descriptor **self, int ip1, int ip2, double ptop_8, double pref_8, float rcoef1,
			     double *a_m_8, double *b_m_8, int *ip1_m, int nk);
int Cvgd_new_build_vert_5002(vgrid_descriptor **self, int ip1, int ip2, double ptop_8, double pref_8, float rcoef1, float rcoef2,
			     double *a_m_8, double *b_m_8, double *a_t_8, double *b_t_8, int *ip1_m, int *ip1_t, int nl_m, int nl_t);
int Cvgd_new_build_vert_5005(vgrid_descriptor **self, int ip1, int ip2, double pref_8, float rcoef1, float rcoef2,
			     double *a_m_8, double *b_m_8, double *a_t_8, double *b_t_8, int *ip1_m, int *ip1_t, int nl);
int Cvgd_new_build_vert_5100(vgrid_descriptor **self, int ip1, int ip2, double pref_8, float rcoef1, float rcoef2, float rcoef3, float rcoef4,
			     double *a_m_8, double *b_m_8, double *c_m_8, double *a_t_8, double *b_t_8, double *c_t_8, int *ip1_m, int *ip1_t, int nl);
int Cvgd_new_build_vert_5999(vgrid_descriptor **self, int ip1, int ip2, 
			     double *a_m_8, double *b_m_8, int *ip1_m, int nk);
int Cvgd_new_build_vert_21001(vgrid_descriptor **self, int ip1, int ip2, float rcoef1, float rcoef2, float rcoef3, float rcoef4, 
			      double *a_m_8, double *b_m_8, double *c_m_8, double *a_t_8, double *b_t_8, double *c_t_8, int *ip1_m, int *ip1_t, int nl);
int Cvgd_new_build_vert_21002(vgrid_descriptor **self, int ip1, int ip2, float rcoef1, float rcoef2, float rcoef3, float rcoef4, 
			      double *a_m_8, double *b_m_8, double *c_m_8,
			      double *a_t_8, double *b_t_8, double *c_t_8,
			      double *a_w_8, double *b_w_8, double *c_w_8,
			      int *ip1_m, int *ip1_t, int *ip1_w, int nl);
int Cvgd_new_build_vert2(vgrid_descriptor **self, int kind, int version, int nk, int ip1, int ip2, double *ptop_8, double *pref_8, float *rcoef1, float *rcoef2, float *rcoef3, float *rcoef4,
		     double *a_m_8, double *b_m_8, double *c_m_8, double *a_t_8, double *b_t_8, double *c_t_8, double *a_w_8, double *b_w_8, double *c_w_8, int *ip1_m, int *ip1_t, int *ip1_w, int nl_m, int nl_t, int nl_w);
int Cvgd_vgdcmp(vgrid_descriptor *vgd1, vgrid_descriptor *vgd2);
void Cvgd_free(vgrid_descriptor **self);
int Cvgd_getopt_int(char *key, int *value, int quiet);
int Cvgd_get_int(vgrid_descriptor *self, char *key, int *value, int quiet);
int Cvgd_get_int_1d(vgrid_descriptor *self, char *key, int **value, int *nk, int quiet);
int Cvgd_get_float(vgrid_descriptor *self, char *key, float *value, int quiet);
int Cvgd_get_float_1d(vgrid_descriptor *self, char *key, float **value, int *nk, int quiet);
int Cvgd_put_double(vgrid_descriptor **self, char *key, double value_put);
int Cvgd_get_double(vgrid_descriptor *self, char *key, double *value_get, int quiet);
int Cvgd_get_double_1d(vgrid_descriptor *self, char *key, double **value, int *nk, int quiet);
int Cvgd_get_double_3d(vgrid_descriptor *self, char *key, double **value, int *ni, int *nj, int *nk, int quiet);
int Cvgd_get_char(vgrid_descriptor *self, char *key, char out[], int quiet);
int Cvgd_put_char(vgrid_descriptor **self, char *key, char *value);
int Cvgd_putopt_int(char *key, int value);
int Cvgd_put_int(vgrid_descriptor **self, char *key, int value);
int Cvgd_new_gen2(vgrid_descriptor **self, int kind, int version, float *hyb, int size_hyb, float *rcoef1, float *rcoef2, float *rcoef3, float *rcoef4,
	      double *ptop_8, double *pref_8, double *ptop_out_8,
	      int ip1, int ip2, float *dhm, float *dht, float *dhw, int avg);
int Cvgd_new_gen(vgrid_descriptor **self, int kind, int version, float *hyb, int size_hyb, float *rcoef1, float *rcoef2,
	      double *ptop_8, double *pref_8, double *ptop_out_8,
		 int ip1, int ip2, float *dhm, float *dht, int avg);
int Cvgd_new_gen_1001(vgrid_descriptor **self, float *hyb, int size_hyb, int ip1, int ip2);
int Cvgd_new_gen_2001(vgrid_descriptor **self, float *hyb, int size_hyb, int ip1, int ip2);
int Cvgd_new_gen_5999(vgrid_descriptor **self, float *hyb, int size_hyb, int ip1, int ip2);
int Cvgd_new_gen_1002(vgrid_descriptor **self, float *hyb, int size_hyb, double ptop_8, int ip1, int ip2);
int Cvgd_new_gen_4001(vgrid_descriptor **self, float *hyb, int size_hyb, int ip1, int ip2);
int Cvgd_new_gen_5001(vgrid_descriptor **self, float *hyb, int size_hyb, double ptop_8, double pref_8, float rcoef1, int ip1, int ip2);
int Cvgd_new_gen_5002(vgrid_descriptor **self, float *hyb, int size_hyb, double ptop_8, double pref_8, float rcoef1, float rcoef2, int ip1, int ip2);
int Cvgd_new_gen_5005(vgrid_descriptor **self, float *hyb, int size_hyb, double pref_8, double *ptop_out_8, float rcoef1, float rcoef2, int ip1, int ip2, float dhm, float dht );
int Cvgd_new_gen_5100(vgrid_descriptor **self, float *hyb, int size_hyb, double pref_8, double *ptop_out_8, float rcoef1, float rcoef2, float rcoef3, float rcoef4, int ip1, int ip2, float dhm, float dht, int avg);
int Cvgd_new_gen_21001(vgrid_descriptor **self, float *hyb, int size_hyb, float rcoef1, float rcoef2, float rcoef3, float rcoef4, int ip1, int ip2, float dhm, float dht);
int Cvgd_new_gen_21002(vgrid_descriptor **self, float *hyb, int size_hyb, float rcoef1, float rcoef2, float rcoef3, float rcoef4, int ip1, int ip2, float dhm, float dht, float dhw);

int Cvgd_new_read(vgrid_descriptor **self, int unit, int ip1, int ip2, int kind, int version);
int Cvgd_write_desc (vgrid_descriptor *self, int unit);
int Cvgd_new_from_table(vgrid_descriptor **self, double *table, int ni, int nj, int nk);
int Cvgd_stda76_temp(vgrid_descriptor *self, int *i_val, int nl_t, float *temp);
int Cvgd_stda76_pres(vgrid_descriptor *self, int *i_val, int nl_t, float *pres, float *sfc_temp, float *sfc_pres);
int Cvgd_stda76_hgts_from_pres_list(float *hgts, float *pres, int nb);
int Cvgd_stda76_pres_from_hgts_list(float *pres, float *hgts, int nb);
#endif // VGRID_H
