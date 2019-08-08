  
  int k,*ind,ij,ijk;
  double lvl;
  char message[128];

  strcpy(message,"(Cvgd) ERROR in ");
  strcat(message,proc_name);
  strcat(message,", cannot allocate ind of int of size\n");
  
  if( my_alloc_int(&ind, nk, message) == VGD_ERROR )
    return(VGD_ERROR);
  
  // Find ip1 indexes
  for( k = 0; k < nk; ++k ){
    if( ( ind[k] = VGD_FindIp1Idx(ip1_list[k],self->ip1_m,self->nl_m)) == -1 ) {
      printf("(Cvgd) ERROR in %s, cannot find ip1 %d in vgrid descriptor.\n",proc_name, ip1_list[k]);
      free(ind);
      return(VGD_ERROR);
    }
  }
  
  // Compute heights
  for( k = 0, ijk=0; k < nk; ++k ){
    for( ij = 0; ij < ni*nj; ++ij, ++ijk ){
      lvl = self->a_m_8[ind[k]];
#if defined(REAL_8)
      levels[ijk] = lvl;
#else
      levels[ijk] = (float) lvl;
#endif
    }
  }
  free(ind);
  return(VGD_OK);
