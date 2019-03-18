#ifndef CONVERT_IP_DEFS
#define CONVERT_IP_DEFS

/* C++ needs to know that types and declarations are C, not C++.  */
#ifdef    __cplusplus
# define __DEB_DECL    extern "C" {
# define __FIN_DECL    }
#else
# define __DEB_DECL
# define __FIN_DECL
#endif

#define TO_IP 1
#define TO_RP -1
#define CONVERT_OK 0
#define CONVERT_GUESS 14
#define CONVERT_GOOD_GUESS 2
#define CONVERT_BAD_GUESS 4
#define CONVERT_TERRIBLE_GUESS 8
#define CONVERT_WARNING 32
#define CONVERT_ERROR 64

/*
KIND = 0, hauteur (m) par rapport au niveau de la mer (-20,000 -> 100,000)
KIND = 1, sigma   (0.0 -> 1.0)
KIND = 2, p est en pression (mb)  (0 -> 1100)
KIND = 3, code arbitraire  (-4.8e8 -> 1.0e10)
KIND = 4, hauteur (M) par rapport au niveau du sol    (-20,000 -> 100,000)
KIND = 5, coordonnee hybride        (0.0 -> 1.0)
KIND = 6, coordonnee theta (1 -> 200,000)
KIND =10, temps en heure    (0.0 -> 200,000.0)
KIND =15, reserve (entiers)        
KIND =17, indice x de la matrice de conversion (1.0 -> 1.0e10)
          (partage avec kind=1 a cause du range exclusif
KIND =21, p est en metres-pression  (partage avec kind=5 a cause du range exclusif)
                                                                       (0 -> 1,000,000) fact=1e4
*/
#define KIND_ABOVE_SEA 0
#define KIND_SIGMA 1
#define KIND_PRESSURE 2
#define KIND_ARBITRARY 3
#define KIND_ABOVE_GND 4
#define KIND_HYBRID 5
#define KIND_THETA 6
#define KIND_HOURS 10
#define KIND_SAMPLES 15
#define KIND_MTX_IND 17
#define KIND_M_PRES 21

typedef struct { /* if v1 == v2, it is not a range but a single value */
  float v1;      /* first value of range */
  float v2;      /* second value of range */
  int kind;      /* kind of value (see table above) */
} ip_info;

static ip_info invalid_ip_info={0.0,0.0,-1};
#define NULL_ip_info &invalid_ip_info
#define INIT_ip_info(a) {(a).v1 = 0.0; (a).v2 = 0.0; (a).kind =-1 ;};

/* see fortran module convert_ip123.f90 for quick documentation of arguments */

__DEB_DECL void ConvertIp(int *ip, float *p, int *kind, int mode); __FIN_DECL

__DEB_DECL int EncodeIp( int *ip1, int *ip2, int *ip3, ip_info *p1, ip_info *p2, ip_info *p3); __FIN_DECL
__DEB_DECL int DecodeIp(ip_info *p1, ip_info *p2, ip_info *p3, int ip1, int ip2, int ip3); __FIN_DECL

__DEB_DECL int EncodeIp_v(int ip[3],ip_info p[3]); __FIN_DECL
__DEB_DECL int DecodeIp_v(ip_info p[3],int ip[3]); __FIN_DECL

__DEB_DECL int ConvertPKtoIP(int *ip1, int *ip2, int *ip3, float p1, int kind1, float p2, int kind2, float p3, int kind3); __FIN_DECL
__DEB_DECL int ConvertIPtoPK(float *p1, int *kind1, float *p2, int *kind2, float *p3, int *kind3, int ip1, int ip2, int ip3); __FIN_DECL

__DEB_DECL int ConvertPKtoIP_v(int ip[3],float p[3],int kind[3]); __FIN_DECL
__DEB_DECL int ConvertIPtoPK_v(float p[3],int kind[3],int ip[3]); __FIN_DECL

__DEB_DECL void KindToString(int code, char *s1, char *s2) ; __FIN_DECL

#endif

