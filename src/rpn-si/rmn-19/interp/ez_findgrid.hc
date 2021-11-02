#include <rpnmacros.h>
#include "ezscint.h"
#include "ez_funcdef.h"

#ifdef MUTEX
//JP
extern pthread_mutex_t EZ_MTX;
#endif

int c_ez_findgrid(int grid_index, _Grille *gr)
  {
  wordint gdrow, gdcol, index_found, gr_size, resax, resay;
  wordint found, end_reached, next_index, i,j;
  _Grille *refgd;
    
  if (grid_index == -1)
    {
    return -1;
    }
  

  refgd = gr_list[grid_index];
  
  found = -1;
  index_found = -1;
  end_reached = -1;

  if (gr == NULL)
     {
     fprintf(stderr, "gr = NULL!\n");
     found = -1;
     return -1;
     }

#ifdef MUTEX
// JP
   pthread_mutex_lock(&EZ_MTX);
#endif
     
  while (found == -1 && end_reached == -1)
    {
    if (gr->grtyp[0] == 'U')
        {
          if (gr->grtyp[0] == refgd->grtyp[0] &&
              gr->grref[0] == refgd->grref[0] &&
              gr->ni_ax == refgd->ni_ax &&  gr->nj_ay == refgd->nj_ay &&
       gr->fst.ig[IG1] == refgd->fst.ig[IG1] &&
       gr->fst.ig[IG2] == refgd->fst.ig[IG2] &&
       gr->fst.ig[IG3] == refgd->fst.ig[IG3] &&
       gr->fst.ig[IG4] == refgd->fst.ig[IG4])
              {
               found = 1;
               index_found = refgd->index;
               break;
              }
        }
     /*  printf("gr->grtyp=%c  ni=%d nj= %d\n",gr->grtyp[0],gr->ni,gr->nj);
         printf("refgd->grtyp=%c  ni=%d nj= %d\n",refgd->grtyp[0],refgd->ni,refgd->nj);
     */
  if (gr->grtyp[0] == refgd->grtyp[0] &&
      gr->ni == refgd->ni &&  gr->nj == refgd->nj &&
      gr->fst.ig[IG1] == refgd->fst.ig[IG1] && gr->fst.ig[IG2] == refgd->fst.ig[IG2] &&
      gr->fst.ig[IG3] == refgd->fst.ig[IG3] && gr->fst.ig[IG4] == refgd->fst.ig[IG4])
      {
      if (gr->grtyp[0] == 'G')
          {
           found = 1;
           index_found = refgd->index;
           break;
           }
      if (gr->fst.igref[IG1] == refgd->fst.igref[IG1] && 
         gr->fst.igref[IG2] == refgd->fst.igref[IG2] &&
         gr->fst.igref[IG3] == refgd->fst.igref[IG3] && 
         gr->fst.igref[IG4] == refgd->fst.igref[IG4])
      {
      if (refgd->ax != NULL && refgd->ay != NULL && gr->ax != NULL && gr->ay != NULL)
         {
         resax = 0;
         resay = 0;
         if (refgd->grtyp[0] == 'Z')
            {
            for (i=0; i < gr->ni; i++)
               {
               if (refgd->ax[i] != gr->ax[i])
                 {
                 resax=1;
                 break;
                 }
               }
            if (resax == 0)
               {
               for (j=0; j < gr->nj; j++)
                  {
                  if (refgd->ay[j] != gr->ay[j])
                    {
                    resay=1;
                    break;
                    }
                  }
               }
/*            resax = memcmp(refgd->ax, gr->ax, (size_t)(gr->ni*sizeof(ftnfloat)));
            resay = memcmp(refgd->ay, gr->ay, (size_t)(gr->nj*sizeof(ftnfloat)));*/
            }
         else
            {
            for (i=0; i < gr->ni*gr->nj; i++)
               {
               if (refgd->ax[i] != gr->ax[i])
                 {
                 resax=1;
                 break;
                 }
               }
            if (resax == 0)
               {
               for (j=0; j < gr->ni*gr->nj; j++)
                  {
                  if (refgd->ay[j] != gr->ay[j])
                    {
                    resay=1;
                    break;
                    }
                  }
               }

/*            resax = memcmp(refgd->ax, gr->ax, (size_t)(gr->ni*gr->nj*sizeof(ftnfloat)));
            resay = memcmp(refgd->ay, gr->ay, (size_t)(gr->ni*gr->nj*sizeof(ftnfloat)));*/
            }

         if (resax == 0 && resay == 0)
            {
            found = 1;
            index_found = refgd->index;
            }
         else
            {
            if (refgd->next_gd == -1)
               {
               end_reached = 1;
               }
            else
               {
               c_gdkey2rowcol(refgd->next_gd, &gdrow, &gdcol);
               refgd = &Grille[gdrow][gdcol];
               }
            }
         }
      else
         {
         found = 1;
         index_found = refgd->index;
         }
      }
      }
    else
      {
      if (refgd->next_gd == -1)
        {
        end_reached = 1;
        }
      else
        {
        c_gdkey2rowcol(refgd->next_gd, &gdrow, &gdcol);
        refgd = &Grille[gdrow][gdcol];
        }
      }
    }  

// JP
  refgd->access_count++;
#ifdef MUTEX
// JP
   pthread_mutex_unlock(&EZ_MTX);
#endif
  
  if (found == -1)
    {
    return -1;
    }
  else  
    {
    return index_found;
    }
  }

void dump_gr_list()
  {
  int i, gd_row, gd_col;
  _Grille *gr;
  
  for (i=0; i < chunks[cur_log_chunk]; i++)
    {
    if (gr_list[i] != NULL)
      {
      gr = gr_list[i];
      printf("%d %d -> ", i, gr->index);
      while (gr->next_gd != -1)
        {
        printf("%d ->", gr->next_gd);
        c_gdkey2rowcol(gr->next_gd, &gd_row, &gd_col);
        gr = &Grille[gd_row][gd_col];
        }
      printf("\n");
      }
    }
  }
  
