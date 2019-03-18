
  double *aa_8, *bb_8, *cc_8;
  int ij, k, ijk, ind, kind;
  float hyb;

  aa_8 = malloc(nk*sizeof(double));
  if(! aa_8 ) {
    printf("(Cvgd) ERROR in %s, cannot allocate aa_8 of bouble of size %d\n", proc_name, nk);
    return(VGD_ERROR);
  }  
  bb_8 = malloc(nk*sizeof(double));
  if(! bb_8 ) {
    printf("(Cvgd) ERROR in %s, cannot allocate bb_8 of bouble of size %d\n", proc_name, nk);
    free(aa_8);
    return(VGD_ERROR);
  }
  cc_8 = malloc(nk*sizeof(double));
  if(! cc_8 ) {
    printf("(Cvgd) ERROR in %s, cannot allocate cc_8 of bouble of size %d\n", proc_name, nk);
    free(aa_8);
    free(bb_8);
    return(VGD_ERROR);
  }

  for(k=0; k < nk; k++) {
    if( (ind = VGD_FindIp1Idx( ip1_list[k], self->ip1_m, self->nl_m) ) != -1 ) {
      aa_8[k] = self->a_m_8[ind];
      bb_8[k] = self->b_m_8[ind];
      cc_8[k] = self->c_m_8[ind];
    } else {
      if( (ind = VGD_FindIp1Idx( ip1_list[k], self->ip1_t, self->nl_t) ) != -1 ) {
	aa_8[k] = self->a_t_8[ind];
	bb_8[k] = self->b_t_8[ind];
	cc_8[k] = self->c_t_8[ind];
      } else {
	if( (ind = VGD_FindIp1Idx( ip1_list[k], self->ip1_w, self->nl_w) ) != -1 ) {
	  aa_8[k] = self->a_w_8[ind];
	  bb_8[k] = self->b_w_8[ind];
	  cc_8[k] = self->c_w_8[ind];
	} else {
	  printf("(Cvgd) ERROR in %s, cannot find ip1 %d in vgrid descriptor.\n", proc_name,ip1_list[k]);
	  free(aa_8);
	  free(bb_8);
	  free(cc_8);
	  return(VGD_ERROR);	
	}
      }
    }
  }

  if( strcmp(self->ref_namel, VGD_NO_REF_NOMVAR) == 0 ){
    my_sfc_field_ls = sfc_field;
  } else {
    my_sfc_field_ls = sfc_field_ls;
  }

  for(k=0, ijk=0; k < nk; k++) {
    for(ij=0; ij < ni*nj; ij++, ijk++) {
#if defined(REAL_8)
      levels[ijk] = aa_8[k] + bb_8[k]*sfc_field[ij] + cc_8[k]*my_sfc_field_ls[ij];
#else
      levels[ijk] = (float) ( aa_8[k] + bb_8[k]*sfc_field[ij] + cc_8[k]*my_sfc_field_ls[ij] );
#endif
    }
  }
  //Force surface heights to be equal to sfc_field
  //Needed by assimilation section.  
  for(k=0; k < nk; k++) {
    hyb = c_convip_IP2Level(ip1_list[k],&kind);
    if(fabs(hyb) < .000001 && kind == 21) {
      ijk=k*ni*nj;
      for(ij=0; ij < ni*nj; ij++, ijk++) {
        levels[ijk] = sfc_field[ij];
      }
    }
  }

  free(aa_8);
  free(bb_8);
  free(cc_8);

  return(VGD_OK);
