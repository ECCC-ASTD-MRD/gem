#include <rpnmacros.h>
#include "ezscint.h"
#include "ez_funcdef.h"

int c_ez_addgrid(int grid_index, _Grille *newgr)
  {
  int i, gdrow, gdcol, gdindex, next_index, nxt_row, nxt_col, cur_gdid;
  _Grille *cur_gr;
  
  cur_gr = gr_list[grid_index];
  if (cur_gr == NULL)
    {
    gdindex = nGrilles;
    c_gdkey2rowcol(gdindex, &gdrow, &gdcol);
    gr_list[grid_index] = &Grille[gdrow][gdcol];
    }
  else
    {
    cur_gdid = cur_gr->index;
    c_gdkey2rowcol(cur_gdid, &gdrow, &gdcol);
    next_index = Grille[gdrow][gdcol].next_gd;
    nxt_row = gdrow;
    nxt_col = gdcol;
    while (next_index != -1)
      {
      c_gdkey2rowcol(next_index, &nxt_row, &nxt_col);
      next_index = Grille[nxt_row][nxt_col].next_gd;      
      }
    Grille[nxt_row][nxt_col].next_gd = nGrilles; 
    }
    
  c_gdkey2rowcol(nGrilles, &gdrow, &gdcol);
  memcpy(&(Grille[gdrow][gdcol]), newgr, sizeof(_Grille));
  Grille[gdrow][gdcol].index = nGrilles;
  Grille[gdrow][gdcol].next_gd = -1;
  
  nGrilles++;
  if (nGrilles >= chunks_sq[cur_log_chunk])
    {
    fprintf(stderr, "<c_ez_addgrid> : Message from the EZSCINT package\n");
    fprintf(stderr, "<c_ez_addgrid> : Maximum number of definable grids attained : %d\n", nGrilles);
    fprintf(stderr, "               : Please contact RPN support to increase the maximum number\n");
    exit(13);
    }
  
  if (0 == (nGrilles % chunks[cur_log_chunk]))
    {
    c_gdkey2rowcol(nGrilles, &gdrow, &gdcol);
    Grille[gdrow] = (_Grille *) calloc(chunks[cur_log_chunk], sizeof(_Grille));
    for (i=0; i < chunks[cur_log_chunk]; i++)
      {
      Grille[gdrow][i].index = -1;
      }
    }
  return (nGrilles-1);    
  }
