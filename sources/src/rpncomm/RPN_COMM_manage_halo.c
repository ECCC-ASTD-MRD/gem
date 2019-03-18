#include <stdlib.h>
#include <stdio.h>

#define MAX_TABLE 256
#define MAX_HASH 4
static struct {
  void *p;
  short Hx;
  short Hy;
  short Hz;  /* hz had to become Hz because of xlc */
  short spare;
  int hash_east[MAX_HASH];
  int hash_west[MAX_HASH];
  int hash_south[MAX_HASH];
  int hash_north[MAX_HASH];
} table[MAX_TABLE];

static int to_initialize=1;
static int table_full_error=0;

static void init(){
  int i;
  if(to_initialize){
    to_initialize = 0;
    for (i=0 ; i<MAX_TABLE ; i++) {
      table[i].p=NULL;
      table[i].Hx=-1;
      table[i].Hy=-1;
      table[i].Hz=-1;
    }
  }
}
/*
 * manage valid halo table
 * the array address is used as a tag
 * 
 * integer functions, return 0 if successful, -1 otherwise
 * 
 * from Fortran
 *   status = rpn_comm_get_valid_halo(array,halox,haloy,haloz)
 *     integer, intent(OUT) :: halox,haloy,haloz
 *   status = rpn_comm_set_valid_halo(array,halox,haloy,haloz)
 *     integer, intent(IN) :: halox,haloy,haloz
 *     a negative value in any of halox,haloy,haloz will invalidate the entry
 *   status = rpn_comm_inc_valid_halo(array,halox,haloy,haloz)
 *     integer, intent(IN) :: halox,haloy,haloz
 *     increment halos by halox,haloy,haloz (negative values will decrement)
 *
 * from C
 *   int rpn_comm_get_valid_halo(void *array,int *Hx, int*Hy, int *Hz)
 *   int rpn_comm_set_valid_halo(void *array,int *Hx, int*Hy, int *Hz)
 *   int rpn_comm_inc_valid_halo(void *array,int *Hx, int*Hy, int *Hz)
 */
#pragma weak rpn_comm_get_valid_halo_=rpn_comm_get_valid_halo
int rpn_comm_get_valid_halo_(void *array,int *Hx, int*Hy, int *Hz);
#pragma weak rpn_comm_get_valid_halo__=rpn_comm_get_valid_halo
int rpn_comm_get_valid_halo__(void *array,int *Hx, int*Hy, int *Hz);

int rpn_comm_get_valid_halo(void *array,int *Hx, int*Hy, int *Hz){
  int i;
  if(to_initialize)init();
  *Hx = -1;
  *Hy = -1;
  *Hz = -1;
  for (i=0 ; i<=MAX_TABLE ; i++){
    if(table[i].p==array) {
      *Hx = table[i].Hx;
      *Hy = table[i].Hy;
      *Hz = table[i].Hz;
      return 0;
    }
  }
  return -1;  /* array not found */
}

#pragma weak rpn_comm_inc_valid_halo_=rpn_comm_inc_valid_halo
int rpn_comm_inc_valid_halo_(void *array,int *Hx, int*Hy, int *Hz);
#pragma weak rpn_comm_inc_valid_halo__=rpn_comm_inc_valid_halo
int rpn_comm_inc_valid_halo__(void *array,int *Hx, int*Hy, int *Hz);

int rpn_comm_inc_valid_halo(void *array,int *Hx, int*Hy, int *Hz){
  int i;
  if(to_initialize)init();
  for (i=0 ; i<=MAX_TABLE ; i++){
    if(table[i].p==array) {
      table[i].Hx += *Hx; if(table[i].Hx<0) table[i].Hx=0;
      table[i].Hy += *Hy; if(table[i].Hy<0) table[i].Hy=0;
      table[i].Hz += *Hz; if(table[i].Hz<0) table[i].Hz=0;
      return 0;
    }
  }
  return -1;  /* array not found */
}

#pragma weak rpn_comm_set_valid_halo_=rpn_comm_set_valid_halo
int rpn_comm_set_valid_halo_(void *array,int *Hx, int*Hy, int *Hz);
#pragma weak rpn_comm_set_valid_halo__=rpn_comm_set_valid_halo
int rpn_comm_set_valid_halo__(void *array,int *Hx, int*Hy, int *Hz);

int rpn_comm_set_valid_halo(void *array,int *Hx, int*Hy, int *Hz){
  int i;
  if(to_initialize)init();
  for (i=0 ; i<=MAX_TABLE ; i++){
    if(table[i].p==array) {
      table[i].Hx = *Hx;
      table[i].Hy = *Hy;
      table[i].Hz = *Hz;
      if(*Hx<0 || *Hy<0 || *Hz<0) table[i].p=NULL;  /* cancel entry in table if halo not valid */
      return 0;
    }
  }
  if(*Hx<0 || *Hy<0 || *Hz<0) return 0;        /* array not found and it was an invalidate command */
  for (i=0 ; i<=MAX_TABLE ; i++){              /* array not found, see if we have a free slot */
    if(table[i].p==NULL) {
      table[i].p  = array;
      table[i].Hx = *Hx;
      table[i].Hy = *Hy;
      table[i].Hz = *Hz;
      return 0;
    }
  }
  if( ! table_full_error ) fprintf(stderr,"WARNING: valid_halo_table is full\n");
  table_full_error=1;
  return -1;  /* halo table overflow */
}

#if defined(SELF_TEST)
main() {
  int a[10];
  int i;
  int Hx,Hy,Hz;
  int status;

  for (i=0 ; i<10 ; i+=2){
    Hx = i+1;
    Hy = i+2;
    Hz = i+3;
    status = rpn_comm_set_valid_halo(&a[i],&Hx,&Hy,&Hz);
    Hx = +1;
    Hy = 0;
    Hz = -1;
    status = rpn_comm_inc_valid_halo(&a[i],&Hx,&Hy,&Hz);
  }
  for (i=0 ; i<10 ; i++){
    status = rpn_comm_get_valid_halo(&a[i],&Hx,&Hy,&Hz);
    if(status==-1){
      fprintf(stderr,"no halo valid for a[%d]\n",i);
    }else{
      fprintf(stderr,"halo for a[%d] = (%d,%d,%d)\n",i,Hx,Hy,Hz);
    }
  }
  for (i=0 ; i<10 ; i+=2){
    Hx = -1;
    Hy = -1;
    Hz = -1;
    status = rpn_comm_set_valid_halo(&a[i],&Hx,&Hy,&Hz);
  }
  for (i=0 ; i<10 ; i++){
    status = rpn_comm_get_valid_halo(&a[i],&Hx,&Hy,&Hz);
    if(status==-1){
      fprintf(stderr,"no halo valid for a[%d]\n",i);
    }else{
      fprintf(stderr,"halo for a[%d] = (%d,%d,%d)\n",i,Hx,Hy,Hz);
    }
  }
  for (i=1 ; i<10 ; i+=2){
    Hx = i;
    Hy = i+1;
    Hz = i+2;
    status = rpn_comm_set_valid_halo(&a[i],&Hx,&Hy,&Hz);
  }
  for (i=0 ; i<10 ; i++){
    status = rpn_comm_get_valid_halo(&a[i],&Hx,&Hy,&Hz);
    if(status==-1){
      fprintf(stderr,"no halo valid for a[%d]\n",i);
    }else{
      fprintf(stderr,"halo for a[%d] = (%d,%d,%d)\n",i,Hx,Hy,Hz);
    }
  }
}
#endif
