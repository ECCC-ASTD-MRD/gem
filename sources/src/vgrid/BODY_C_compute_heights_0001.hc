
  double *aa_8;
  int ij, k, ijk, ind;

  aa_8 = malloc(nk*sizeof(double));
  if(! aa_8 ) {
    printf("(Cvgd) ERROR in %s, cannot allocate aa_8 of bouble of size %d\n", proc_name, nk);
    return(VGD_ERROR);
  }  

  for(k=0; k < nk; k++) {
    if( (ind = VGD_FindIp1Idx( ip1_list[k], self->ip1_m, self->nl_m) ) != -1 ) {
      aa_8[k] = self->a_m_8[ind];
    } else {
      if( (ind = VGD_FindIp1Idx( ip1_list[k], self->ip1_w, self->nl_w) ) != -1 ) {
	aa_8[k] = self->a_w_8[ind];
      } else {
	printf("(Cvgd) ERROR in %s, cannot find ip1 %d in vgrid descriptor.\n", proc_name,ip1_list[k]);
	free(aa_8);
	return(VGD_ERROR);	
      }
    }
  }

  for(k=0, ijk=0; k < nk; k++) {
    for(ij=0; ij < ni*nj; ij++, ijk++) {
#if defined(REAL_8)
      levels[ijk] = aa_8[k];
#else
      levels[ijk] = (float) aa_8[k];
#endif
    }
  }

  free(aa_8);

  return(VGD_OK);
