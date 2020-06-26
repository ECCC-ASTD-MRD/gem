//  Hopefully useful software for FORTRAN and C
//  Copyright (C) 2019  Division de Recherche en Prevision Numerique
//                      Environnement Canada
// 
//  This is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation,
//  version 2.1 of the License.
// 
//  This software is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <string.h>

#define NEED_PRIVATE
#include <time_trace.h>
#undef NEED_PRIVATE

long long WhatTimeIsIt(){           // get current time of day in microseconds
  struct timeval tv;
  long long time;
  gettimeofday(&tv, NULL);
  time = tv.tv_sec;                 // seconds
  time *= 1000000;                  // seconds to microseconds
  time += tv.tv_usec;               // add microseconds part
  return time;
}

// tp : pointer to a time context structure
void TimeTraceInit(time_context *tp){  // create initialized time context
  trace_table *tt = (trace_table*)TimeTraceCreate();  // create trace table
  tp->t = tt;                          // return time context to caller
}

void *TimeTraceCreate(){            // create and initialize a trace_table
  trace_table *tt = (trace_table*) malloc(sizeof(trace_table));
  time_context t;

  tt->initialized = 0;              // mark as not initialized
  tt->nbeads      = 0;              // 0 "beads" in table
  tt->nwords      = 0;              // nothing stored in table
  tt->major       = MAJORV;         // major version marker
  tt->minor       = MINORV;         // minor version marker
  tt->step        = 99999999;       // default "step" number
  tt->offset      = WhatTimeIsIt(); // initial offset
  tt->first       = NULL;           // pointer to first "bead"
  tt->last        = NULL;           // pointer to last(current) "bead"
  t.t            = tt;
  TimeTraceCreateNewBead(t);        // create first "bead" in table

  return tt;
}

// t : time context (obtained from TimeTraceInit)
void TimeTraceCreateNewBead(time_context t){   // add a new timing "bead"
  trace_table *tt = t.t;
  bead *temp = (bead *) malloc(sizeof(bead));  // allocate a new "bead"
  bead *last;

  if(tt->initialized == 0){      // not initialized yet, link first bead
    tt->first = temp;
    tt->initialized = 1;         // mark as initialized
  }else{
    last = (bead *) tt->last;    // link new "baed" after the last one
    last->next = temp;
  };
  temp->next  = NULL;            // no next "bead"
  temp->nbent = 0;               // "bead" is empty
  temp->mxent = MAX_TIMES;       // max number of items that can be stored in "bead"
  tt->last    = temp;            // last "bead" is the newly created one
  tt->nbeads += 1;               // increment "bead" counter 
}

// t   : time context (obtained from TimeTraceInit)
// val : 32 integer value
void TimeTraceInsert(time_context t, int val){  // insert a 32 bit value in table
  trace_table *tt = t.t;
  bead *last;

  last = (bead *) tt->last;
  if(last->nbent >= last->mxent) TimeTraceCreateNewBead(t);  // last "bead" full, add a new one
  last = (bead *) tt->last;
  last->t[last->nbent] = val;      // insert value
  tt->nwords  += 1;                // bump global insertion counter
  last->nbent += 1;                // bump "bead" insertion counter
}

// t   : time context (obtained from TimeTraceInit)
// tag : 32 integer value, must be >0 and <99999999
// t1  : 64 bit intgeger time value
// t2  : 64 bit intgeger time value
// if t2 < 0, it will be ignored
void TimeTraceNewTag(time_context t, int tag, long long t1, long long t2){  // insert a tag and one or two times
  trace_table *tt = t.t;
  int tm1, tm2, code;

  tm1 = t1 - tt->offset;        // remove current "step" offset from t1
  if(t2 >= 0){                  // t2 < 0 means t2 is to be ignored
    tm2 = t2 - tt->offset;      // remove current "step" offset from t2
    code = (tag << 3) | 2;      // 2 values will follow tag
  }else{
    code = (tag << 3) | 1;      // one value will follow tag
  }
  TimeTraceInsert(t, code);     // insert tag and value(s)
  TimeTraceInsert(t, tm1);
  if(t2 >= 0) TimeTraceInsert(t, tm2);
}

// t    : time context (obtained from TimeTraceInit)
// step : 32 integer value, must be >0 and <99999999
void TimeTraceNewStep(time_context t, int step){ // insert a new "step" (timing epoch) 
  trace_table *tt = t.t;
  long long t0;
  int code, thi, tlo;

  t0 = WhatTimeIsIt();            // get time of day in microseconds
  thi = t0 >> 32;
  tlo = t0 & 0xFFFFFFFF;
  tt->offset = t0;
  code = step << 3;
  TimeTraceInsert(t, code);       // coded "step" value and full 64 bit time
  TimeTraceInsert(t, thi);        // upper 32 bits of time
  TimeTraceInsert(t, tlo);        // lower 32 bits of time
}

// t   : time context (obtained from TimeTraceInit)
// tag : 32 integer value, must be >0 and <99999999
void TimeTraceBarr(time_context t, int tag, int barrier, int fcomm, void (*barrier_code)()){   // insert tag and one time value
  long long t0 = WhatTimeIsIt();
  long long t1 = -1;              //  make sure t1 will be ignored

  if(barrier != 0) t1 = WhatTimeIsIt();
  TimeTraceNewTag(t, tag, t0, t1);
}

// t   : time context (obtained from TimeTraceInit)
long long  TimeTraceSinceStep(time_context t){   // time since last step call
  long long t0 = WhatTimeIsIt();
  trace_table *tt = t.t;

  return (t0 - tt->offset);
}

// t   : time context (obtained from TimeTraceInit)
// tag : 32 integer value, must be >0 and <99999999
void TimeTrace(time_context t, int tag){   // insert tag and one time value
  long long t0 = WhatTimeIsIt();
  long long t1 = -1;              //  make sure t1 will be ignored

  TimeTraceNewTag(t, tag, t0, t1);
}

// t    : time context (obtained from TimeTraceInit)
// step : 32 integer value, must be >0 and <99999999
void TimeTraceStep(time_context t, int step){   // start a new "step" (timing epoch)
  TimeTraceNewStep(t, step);
}

// get addresses and lengths of timing "beads"
// addresses will be stored into array
// lengths will be stored into larray
// size of array and larray is n
// if n is larger than the number of "beads", NULL pointers and ZERO lengths will be stored into the end part of array and larray
// t : time context (obtained from TimeTraceInit)
void TimeTraceGetBuffers(time_context t, void **array, int *larray, int n){  
  trace_table *tt = t.t;
  bead *current;
  int i;

  current = (bead *)tt->first;         // head of chain
  for(i=0 ; i<n ; i++){
    if(current == NULL){               // no more "beads", store NULL and ZERO
      array[i]  = NULL;
      larray[i] = 0;
    }else{                             //store "bead" address and length
      array[i]  = current->t;
      larray[i] = current->nbent;
      current = (bead *)current->next;
    }
  }
}

// get number of time trace buffers and total number of entries (useful for call to TimeTraceGetBuffers)
// t           : time context (obtained from TimeTraceInit)
// nbuf        : number of buffers
// nent        : total number of entries
// consolidate : if non zero consolidate data into a single buffer (allocated locally)
void  *TimeTraceGetBufferData(time_context t, int *nbuf, int *nent, int consolidate){  
  trace_table *tt = t.t;
  bead *current;
  int nbent;
  unsigned int *data = NULL;
  unsigned int *temp;

  current = (bead *)tt->first;         // head of chain
  *nbuf = tt->nbeads;
  nbent = tt->nwords + HEADER;
  *nent = nbent;

  if(consolidate){
    data = (unsigned int *) malloc((nbent+HEADER) * sizeof(unsigned int));   // allocate single buffer
    data[0] = MAJORV;
    data[1] = MINORV;
    data[2] = tt->nbeads;
    data[3] = nbent - HEADER;
    temp = data + HEADER;
    if(data != NULL){                                               // malloc successful
      current = (bead *)tt->first;                                  // first bead
      while(current != NULL){                                       // loop over valid beads
        nbent = current->nbent;                                     // entries in this bead
        memcpy(temp, current->t, nbent * sizeof(unsigned int));     // copy to consolidation buffer
        temp += nbent;                                              // bump destination pointer
        current = (bead *)current->next;                            // next bead
      }
    }
  }
  return data;
}

// t        : time context (obtained from TimeTraceInit)
// filename : base part of file name
// ordinal  : integer (will be converted as 6 digit nnnnnn)
// file name will be filename_nnnnnn.dat
void TimeTraceDumpBinary(time_context t, char *filename, int ordinal){  // dump into file in binary form
  trace_table *tt = t.t;
  bead *current;
  char fname[4096];
  int fd, i, nw;
  ssize_t  nwr;

  snprintf(fname, sizeof(fname)-1, "%s_%6.6d.dat", filename, ordinal);
  fd = open(fname, O_RDWR | O_CREAT, 0644);
  current = (bead *)tt->first;
  for(i=0 ; i < tt->nbeads && current != NULL; i++){
    nw = current->nbent;
    nwr = write(fd, current->t, sizeof(int)*nw);
    if(nwr != sizeof(int)*nw) exit(1);
    current = (bead *)current->next;
  }
  close(fd);
}

// t        : time context (obtained from TimeTraceInit)
// filename : base part of file name
// ordinal  : integer (will be converted as 6 digit nnnnnn)
// file name will be filename_nnnnnn.txt
void TimeTraceDumpText(time_context t, char *filename, int ordinal){  // dump into file in text form
  trace_table *tt = t.t;
  char fname[4096];
  int cstep = 99999999;
  FILE *fd;
  bead *current;
  int i, tag, nval, j;
  unsigned long long tm[10];
  unsigned long long tm8;

  if(tt == NULL) return;

  snprintf(fname, sizeof(fname)-1, "%s_%6.6d.txt", filename, ordinal);
  fd = fopen(fname, "w");
  fprintf(fd,"%d %d %d %d\n",tt->major,tt->minor,tt->nbeads,tt->nwords);
  i = 0;
  current = (bead *)tt->first;
  while(current != NULL){
    tag = -2; 
    nval = 0; 
    for(j=0 ; j<10 ; j++) tm[j] = 0 ;
    for(j=0 ; j<4 ; j++){
      if(i >= current->nbent){
        current = (bead *)current->next;
        if(current == NULL) return;
        i = 0;
      }
      if(j == 0){       // tag or step
        tag = current->t[i]; nval = tag & 3 ; 
        if(nval == 0){
          cstep = tag >> 3 ; tag = -1; nval = 2; // fprintf(stderr,"step %d nval %d\n",cstep,nval);
        }else{
          tag = tag >> 3; //  fprintf(stderr,"tag   %d nval %d\n",tag,nval);
        }
      }else{           // j > 0, data
        tm[j-1] = current->t[i];
      }
      i++;
      if(j >= nval) break;  // done collecting
    }
    if(tag == -1) {                                          // step 
      tm8 = (tm[0] << 32) | tm[1] ;
      fprintf(fd,"%d %d %llu %d",cstep, tag, tm8, 0);
    }else{                                                   // tag with or without barrier
      fprintf(fd,"%d %d ",cstep, tag);
      if(nval == 1) { nval = 2 ; tm[1] = 0 ; }
      for(j=0 ; j<nval ; j++) fprintf(fd,"%llu ",tm[j]);
    }
    fprintf(fd,"\n");
  }
  fclose(fd);
}

// data     : data blob from TimeTraceGetBufferData
// nbent    : length of data
// filename : base part of file name
// ordinal  : integer (will be converted as 6 digit nnnnnn)
// file name will be filename_nnnnnn.txt
void TimeTraceSingleText(unsigned int *data, int nbent, char *filename, int ordinal){  // dump blob into file in text form
  char fname[4096];
  int cstep = 99999999;
  FILE *fd;
  int i, tag, nval, j;
  unsigned long long tm[10];
  unsigned long long tm8;

  nbent = nbent - HEADER;
  if(data == NULL || nbent <= 0) return;
  if(data[0] != MAJORV || data[1] != MINORV || data[3] != nbent) {
    fprintf(stderr,"ERROR: metadata is not consistent, expected %d %d %d, got %d, %d, %d\n",MAJORV,MINORV,nbent,data[0],data[1],data[3]);
    return;
  }

  snprintf(fname, sizeof(fname)-1, "%s_%6.6d.txt", filename, ordinal);
  fd = fopen(fname, "w");
  fprintf(fd,"%d %d %d %d\n",data[0], data[1], data[2], data[3]);
  data = data + HEADER;        // skip metadata at start of data buffer
//   fprintf(fd,"%d %d %d %d\n",tt->major,tt->minor,tt->nbeads,tt->nwords);

  i = 0;
  while(i < nbent){
    tag = -2; 
    nval = 0; 
    for(j=0 ; j<10 ; j++) tm[j] = 0 ;
    for(j=0 ; j<4 ; j++){
      if(i >= nbent) goto the_end;
      if(j == 0){       // tag or step
        tag = data[i]; 
        nval = tag & 3 ; 
        if(nval == 0){
          cstep = tag >> 3 ; tag = -1; nval = 2; // fprintf(stderr,"step %d nval %d\n",cstep,nval);
        }else{
          tag = tag >> 3; //  fprintf(stderr,"tag   %d nval %d\n",tag,nval);
        }
      }else{           // j > 0, data
        tm[j-1] = data[i];
      }
      i++;
      if(j >= nval) break;  // done collecting
    }
    if(tag == -1) {                                          // step 
      tm8 = (tm[0] << 32) | tm[1] ;
      fprintf(fd,"%d %d %llu %d",cstep, tag, tm8, 0);
    }else{                                                   // tag with or without barrier
      fprintf(fd,"%d %d ",cstep, tag);
      if(nval == 1) { nval = 2 ; tm[1] = 0 ; }
      for(j=0 ; j<nval ; j++) fprintf(fd,"%llu ",tm[j]);
    }
    fprintf(fd,"\n");
  }
the_end:
  fclose(fd);
}

// data     : data blob from TimeTraceGetBufferData
// nbent    : length of data
// ref      : integer array, [lines][linesize] in C,   (linesize,lines) in Fortran
// out      : integer array, [lines][linesize] in C,   (linesize,lines) in Fortran
// linesize : "line" size
// lines    : number of linesize integer elements items in out
// ver      : ref is used as input for consistency checking if ver not zero, as output if ver is zero
// the function returns the number of 'lines' used in out
// if out is too small, the functions return - (number of 'lines')
int TimeTraceExpand(unsigned int *data, int nbent, int *ref, int *out, int linesize, int lines, int ver){  // expand blob from TimeTraceGetBufferData
  int cstep = 99999999;
  int i, tag, nval, j, nl;
  unsigned long long tm[10];

  nbent = nbent - HEADER;
  if(data == NULL || nbent <= 0) return 0;
  if(data[0] != MAJORV || data[1] != MINORV || data[3] != nbent) {
    fprintf(stderr,"ERROR: metadata is not consistent, expected %d %d %d, got %d, %d, %d\n",MAJORV,MINORV,nbent,data[0],data[1],data[3]);
    return 0;
  }
  data = data + HEADER;        // skip metadata at start of data buffer

  i = 0;
  nl = 0;
  while(i < nbent){
    tag = -2; 
    nval = 0; 
    for(j=0 ; j<10 ; j++) tm[j] = 0 ;
    for(j=0 ; j<4 ; j++){
      if(i >= nbent) goto the_end;
      if(j == 0){       // tag or step
        tag = data[i]; 
        nval = tag & 3 ; 
        if(nval == 0){
          cstep = tag >> 3 ; tag = -1; nval = 2;
        }else{
          tag = tag >> 3;
        }
      }else{           // j > 0, data
        tm[j-1] = data[i];
      }
      i++;
      if(j >= nval) break;  // done collecting
    }
    nl++;
    if(nl > lines) return ( -(nl-1) ) ; // OOPS, out is too small
    if(ver == 0){         // ref is an output
      ref[0] = cstep;
      ref[1] = tag;
    }
    if(ref[0] == cstep){  // check that cstep is consistent with ref
      out[0] = tm[0];
      out[1] = tm[1];
    }else{                // inconsistent, flag as such
      out[0] = -1;
      out[1] = -1;
    }
    out += linesize;
    ref += linesize;
  }
the_end:
  return nl;
}

unsigned long long I8From2I4(unsigned int a, unsigned int b){
  unsigned long long t8;
  t8 = a;
  t8 <<= 32;
  t8 |= b;
  return t8;
}

#if defined(C_SELF_TEST)
int main(){
  time_context t;
  int i, tag;
  void *array[100];
  int larray[100];

  printf("C Test of time trace \n");
  TimeTraceInit(&t);
  printf("-------------------------\n");
  TimeTrace(t, 0);
  TimeTraceStep(t, 0);
  printf("=========================\n");
  TimeTraceBarr(t, 1, 1, 0, NULL);
  for(i=1 ; i<4 ; i++){
    TimeTraceStep(t, i);
    printf("++++ start step %d ++++\n",i);
    tag = 10*i+1;
    TimeTraceBarr(t, tag, 1, 0, NULL);
    printf("++++ end step %d ++++\n",i);
    tag = 10*i+2;
    TimeTrace(t, tag);
  }
  TimeTraceStep(t, 4);
  printf("=========================\n");
   TimeTrace(t, 77777777);
   TimeTraceDumpText(t, "time_list", 0);
   TimeTraceDumpBinary(t, "time_list", 0);
   TimeTraceGetBuffers(t, array, larray, 10);
   for(i=0 ; i<10 ; i++) printf("%d ",larray[i]);
   printf("\n");
  return 0;
}
#endif
