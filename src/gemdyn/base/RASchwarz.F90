!**s/r  RAS -   Restricted Additive Schwarz preconditioner with
!               local solver = "separable" approximate elliptic problem
!

      Subroutine RASchwarz(F_SOL,F_RHS,imin,imax,jmin,jmax, &
                                minx, maxx, miny, maxy,NK)
      use dyn_fisl_options
      use sol         , only: sol_pil_s, sol_pil_n, sol_pil_w, sol_pil_e
      use gem_options , only: G_halox, G_haloy

      use prec
      use glb_ld
      use ptopo

      use, intrinsic :: iso_fortran_env
      implicit none
      integer, intent(in) :: minx, maxx, miny, maxy, NK
      integer, intent(in) :: imin, imax, jmin, jmax
      real(kind=REAL64), dimension(minx:maxx, miny:maxy, NK), intent(in) :: F_RHS
      real(kind=REAL64), dimension(minx:maxx, miny:maxy, NK), intent(out) :: F_SOL
      real(kind=REAL64), dimension(imin:imax, jmin:jmax, NK) :: Rwork_space, Swork_space

      integer :: i0,j0,in,jn
      integer :: ii0,iin,jj0,jjn,ni_val,nj_val,halox,haloy



      i0 = 1    + sol_pil_w
      in = l_ni - sol_pil_e
      j0 = 1    + sol_pil_s
      jn = l_nj - sol_pil_n

      halox=ovlpx
      haloy=ovlpy

      ii0=imin
      iin=imax
      jj0=jmin
      jjn=jmax

      Rwork_space=0.d0
      Rwork_space(i0:in,j0:jn,:)=F_RHS(i0:in,j0:jn,:)
      call rpn_comm_xch_halo_8 (Rwork_space,ii0, iin, jj0, jjn, l_ni, l_nj, NK, &
                                 halox , haloy, G_periodx, G_periody, l_ni,0 )

      if (Ptopo_mycol==1)  ii0  = 0
      if (Ptopo_mycol.eq.Ptopo_npex-2)  iin = l_ni+1
      if (Ptopo_myrow==1)  jj0  = 0
      if (Ptopo_myrow.eq.Ptopo_npey-2)  jjn = l_nj+1

      if (l_east)  iin = l_ni-sol_pil_e
      if (l_west)  ii0 = 1+sol_pil_w
      if (l_north) jjn = l_nj-sol_pil_n
      if (l_south) jj0 = 1+sol_pil_s
      ni_val = iin-ii0+1
      nj_val = jjn-jj0+1

      Swork_space=0.d0
      call pre_jacobi3D2 ( Swork_space(ii0:iin,jj0:jjn,:),Rwork_space(ii0:iin,jj0:jjn,:), &
                                     Prec_xevec_8, ii0, iin, jj0, jjn, ni_val, nj_val, NK,&
                                     Prec_ai_8, Prec_bi_8, Prec_ci_8 )

      F_SOL(i0:in,j0:jn,:)=Swork_space(i0:in,j0:jn,:)

      return
      end

