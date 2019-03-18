#include <rpnmacros.h>
#include "ezscint.h"
#include "ez_funcdef.h"

int c_find_gdin(int gdin, int gdout)
  {
   int i, gdrow_out, gdcol_out, idx_gdin, found; 
   _gridset *gset;
    
   c_gdkey2rowcol(gdout, &gdrow_out, &gdcol_out);
   idx_gdin = Grille[gdrow_out][gdcol_out].idx_last_gdin;
   if (idx_gdin == -1)
      {
      c_ezdefset(gdout, gdin);
      idx_gdin = gdin % primes[Grille[gdrow_out][gdcol_out].log_chunk_gdin];
   
/*   idx_gdin = Grille[gdrow_out][gdcol_out].idx_last_gdin;*/
      }
     
   gset = Grille[gdrow_out][gdcol_out].gset;
   if (gset[idx_gdin].gdin == gdin) 
      {
      return idx_gdin;
      }
   
    idx_gdin = gdin % primes[cur_log_chunk];
    if (gset[idx_gdin].gdin == gdin) 
      {
      return idx_gdin;
      }
    
    i = idx_gdin+1;
    found = 0;   
    while (found == 0 && i != idx_gdin && gset[i].gdin != -1)
      {
      if (gset[i].gdin == gdin)
        {
        found = 1;
        idx_gdin = i;
        }
      i++;
      if (0 == (i % primes[Grille[gdrow_out][gdcol_out].log_chunk_gdin]))
        {
        i = 0;
        }
      }
    
    if (found == 0)
      {
      return -1;
      }
    else
      {
      Grille[gdrow_out][gdcol_out].idx_last_gdin = idx_gdin;
      return idx_gdin;
      }
     
  }
