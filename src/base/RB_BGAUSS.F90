!**s/r  RB_BGAUSS -   RED_BLACK_BlocGAUSS_3D  preconditioner with
!               local solver = "separable" approximate elliptic problem
!

      Subroutine RB_BGAUSS(F_SOL,F_RHS,minx, maxx, miny, maxy,NK)
      use dyn_fisl_options
      use sol         , only: sol_pil_s, sol_pil_n, sol_pil_w, sol_pil_e
      use gem_options , only: G_halox, G_haloy

      use prec
      use glb_ld
      use ptopo

      use, intrinsic :: iso_fortran_env
      implicit none
      integer, intent(in) :: minx, maxx, miny, maxy, NK
      real(kind=REAL64), dimension(minx:maxx, miny:maxy,NK), intent(in) :: F_RHS
      real(kind=REAL64), dimension(minx:maxx, miny:maxy,NK), intent(out) :: F_SOL

      integer :: icol,P_mycol,i0,j0,in,jn,niloc,njloc
      real(kind=REAL64), dimension(l_minx:l_maxx, l_miny:l_maxy,NK) :: Rwork_space
      real(kind=REAL64), dimension(l_minx:l_maxx, l_miny:l_maxy,NK) :: Swork_space
!author
!       Abdessamad Qaddouri -  2018
!


      j0=1+sol_pil_s
      jn=l_nj-sol_pil_n
      i0=1+sol_pil_w
      in=l_ni-sol_pil_e
      niloc= in-i0+1
      njloc= jn-j0+1
      P_mycol=mod(Ptopo_mycol+Ptopo_myrow,2)+1


      call pre_jacobi3D ( Swork_space(i0:in,j0:jn,:),F_RHS(i0:in,j0:jn,:), &
                                      Prec_xevec_8, niloc, njloc, NK,&
                                        Prec_ai_8, Prec_bi_8, Prec_ci_8 )

      do icol =1,Gauss_Niter
              call rpn_comm_xch_halo_8 (Swork_space,l_minx, l_maxx, l_miny, l_maxy, l_ni, l_nj, NK, &
                       G_halox, G_haloy, G_periodx, G_periody, l_ni,0 )
              if (icol == P_mycol) then
                  Rwork_space(i0:in,j0:jn,:)=F_RHS(i0:in,j0:jn,:)
                  call bord_cor( Rwork_space,Swork_space ,l_minx,l_maxx,l_miny,l_maxy,l_ni,l_nj,NK )
                  call pre_jacobi3D ( Swork_space(i0:in,j0:jn,:), &
                                   Rwork_space(i0:in,j0:jn,:), &
                                   Prec_xevec_8, niloc, njloc, NK,&
                                   Prec_ai_8, Prec_bi_8, Prec_ci_8 )
              end if
      end do

      F_SOL(i0:in,j0:jn,:)=Swork_space(i0:in,j0:jn,:)
      return
      end

