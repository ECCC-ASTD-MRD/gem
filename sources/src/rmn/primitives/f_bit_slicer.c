#include <rpnmacros.h>
/* FORTRAN callable bit slicer  
   call f_bits_get(bits,bits_per_slice,slices,nslices)
   call f_bits_put(bits,bits_per_slice,slices,nslices)
   integer bits(*)                    bit stream (IN for get) (OUT for put)
   integer bits_per_slice(nslices)    number of bits to extract for slice (IN) (IN)
   integer slices(nslices)            extracted slices (OUT for get) (IN for put)
   integer nslices                    number of slices to extract (IN) (IN)

   ex: call f_bits_get(bits,(/4,8,4,8,4,8,4,8,4,8/),slices,10)

   Note: 32*dimension of array bits must be >= sum of bits_per_slice(i),i=1,nslices
   
   esample test program:
program test_f_bit_slicer
integer bits(4)
integer slices(13)
real rbits(4)
rbits=.1234
slices=-1
bits=-2048
call f_bits_get(bits,(/32,4,0,8,4,8,4,8,4,8,4,8,4/),slices,13)
print *,slices
bits=1234
call f_bits_put(bits,(/32,4,0,8,4,8,4,8,4,8,4,8,4/),slices,13)
print *,bits
call f_bits_get(rbits,(/32,4,0,8,4,8,4,8,4,8,4,8,4/),slices,13)
rbits=1.234
call f_bits_put(rbits,(/32,4,0,8,4,8,4,8,4,8,4,8,4/),slices,13)
print *,rbits
stop
end

expected output:
        -2048           15            0          255           15          248 
            0           15           15          255           15          128 
            0
        -2048        -2048        -2048         1234
   0.1234000       0.1234000       0.1234000        1.234000
   
*/

void f77name(f_bits_get)(ftnword *bit_array, ftnword *bits_per_slice, ftnword *slices, ftnword *nslices)
{
  int number_of_slices = *nslices;
  int left = 0;  /* number of bits left in bit reservoir (mybits) */
  unsigned int mybits;  /* bit reservoir, always left aligned */

  while(--number_of_slices >= 0){
    int current = *bits_per_slice++;                        /* number of bits to extract, min 0, max 32 */
    register unsigned int temp_slice;

    if(left==0) { left = 32; mybits = *bit_array++ ;}       /* refill reservoir if empty */
    if(current > 32) current=32;
    if(current <= 0)                  /* 0 bit slice, the easy answer is 0 */
      temp_slice = 0;
    else if((current==32) && (left==32))
      temp_slice = *bit_array++ ;
    else if(left >= current){         /* enough bits in reservoir to satisfy request */
      register unsigned int temp = ( -1 << (32-current));   /* extraction mask */
      temp = (temp & mybits) >> (32-current);               /* extract and align right */
      temp_slice = temp;
      left -= current;                                      /* decrement reservoir bit count */
      mybits = mybits << current;                           /* left align reservoir */
    }else{                            /* not enough bits in reservoir to satisfy request, will do it in two steps */
      register unsigned int temp = ( -1 << (32-left));      /* extraction mask */
      temp = (temp & mybits) >> (32-left);                  /* extract and align right */
      temp_slice = temp;                                    /* partial slice, right aligned */
      current -= left;                                      /* decrement request count by what we already have */
      left = 32;                                            /* refill reservoir */
      mybits = *bit_array++ ;
      temp = ( -1 << (32-current));                         /* extraction mask */
      temp = (temp & mybits) >> (32-current);               /* extract and align right */
      temp_slice = (temp_slice << current) | temp;          /* add second subslice to partial slice */
      left -= current;                                      /* decrement reservoir bit count */
      mybits = mybits << current;                           /* left align reservoir */
    }
    *slices++ = temp_slice;                                 /* store result */
  }  /* end while */
}
void f77name(f_bits_put)(ftnword *bit_array, ftnword *bits_per_slice, ftnword *slices, ftnword *nslices)
{
  int number_of_slices = *nslices;
  int left = 32;          /* number of bits stored in bit reservoir (mybits) */
  unsigned int mybits=0;  /* bit reservoir, always right aligned */

  while(--number_of_slices >= 0){
    int current = *bits_per_slice++;                        /* number of bits to put, min0, max 32 */
    
    if(current > 32) current=32;
    if(current <= 0)                     /* 0 bit slice, do nothing */
      mybits = mybits;
    else if((current==32) && (left==32))  /* easy case, 32 bit token to insert into empty reservoir */
      *bit_array++ = *slices;
    else if(left >= current){            /* enough space in reservoir to satisfy request */
      register unsigned int temp = ( -1 << (32-current));   /* insertion mask */
      temp = (temp >> (32-current)) & *slices;              /* mask token to insert */
      left -= current;                                      /* decrement reservoir bit count */
      mybits = (mybits << current) | temp;                  /* insert token into right aligned reservoir */
    }else{                               /* not enough bits in reservoir to satisfy request, will do it in two steps */
      register unsigned int temp = ( -1 << (32-left));      /* insertion mask */
      temp = (temp>>(32-left)) & (*slices>>(current-left)); /* mask token to insert */
      mybits = (mybits << left) | temp;                     /* insert token into right aligned reservoir */
      current -= left;                                      /* decrement request count by what we have already inserted */
      left = 32;                                            /* refill reservoir */
      *bit_array++ = mybits;                                /* store bit reservoir */
      temp = ( -1 << (32-current));                         /* insertion mask */
      temp = (temp >> (32-current)) & *slices;              /* mask token to insert */
      left -= current;                                      /* decrement reservoir bit count */
      mybits = temp;                                        /* insert token into right aligned reservoir */
    }
    if(left==0) { *bit_array++ = mybits; left=32 ; mybits=0 ; }
    slices++;                                               /* bump source pointer */
  }  /* end while */
  if(left < 32) *bit_array = (mybits << left) ;              /* left align leftovers and store */
}
