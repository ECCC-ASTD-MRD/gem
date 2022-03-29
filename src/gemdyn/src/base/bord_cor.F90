!
!**s/r  bord_cor -  rhs correction in preconditionner
!
      subroutine  bord_cor (Rhs,lhs, Minx, Maxx, Miny, Maxy, nil, njl, Nk)
      use glb_ld
      use opr
      use sol
      use, intrinsic :: iso_fortran_env
      implicit none
!
      integer, intent(in) :: Minx, Maxx, Miny, Maxy, nil, njl, Nk
      real(kind=REAL64), dimension(Minx:Maxx,Miny:Maxy,Nk), intent(in) :: lhs
      real(kind=REAL64), dimension(Minx:Maxx,Miny:Maxy,Nk), intent(inout) :: Rhs
!
!author
!       Abdessamad Qaddouri - December  2006
!
      integer i,j,k,ii,jj
      real(kind=REAL64) stencil2,stencil3,stencil4,stencil5,di_8
!
       do k = 1,Nk
          i=1+sol_pil_w
          ii=i+l_i0-1
          do j=1+sol_pil_s, njl-sol_pil_n
             jj=j+l_j0-1
             di_8= Opr_opsyp0_8(G_nj+jj) / cos( G_yg_8 (jj) )**2
             stencil2= Opr_opsxp2_8(ii) * Opr_opszp0_8(G_nk+k) * Opr_opsyp0_8(G_nj+jj) &
                           / (cos( G_yg_8 (jj) )**2) / (Opr_opsxp0_8(G_ni+ii)*Opr_opsyp0_8(G_nj+jj))
             Rhs(i,j,k) =Rhs(i,j,k)-stencil2*lhs(i-1,j,k)
          end do
!
          i=nil-sol_pil_e
          ii=i+l_i0-1
          do j=1+sol_pil_s, njl-sol_pil_n
             jj=j+l_j0-1
             di_8= Opr_opsyp0_8(G_nj+jj) / cos( G_yg_8 (jj) )**2
             stencil3=Opr_opsxp2_8(2*G_ni+ii) * Opr_opszp0_8(G_nk+k) * Opr_opsyp0_8(G_nj+jj) &
                          / (cos( G_yg_8 (jj) )**2) / (Opr_opsxp0_8(G_ni+ii) * Opr_opsyp0_8(G_nj+jj))
             Rhs(i,j,k) =Rhs(i,j,k)-stencil3*lhs(i+1,j,k)
          end do
!
          j=1+sol_pil_s
          jj=j+l_j0-1
          do i=1+sol_pil_w, nil-sol_pil_e
             ii=i+l_i0-1
             stencil4=Opr_opsxp0_8(G_ni+ii) * Opr_opsyp2_8(jj) * Opr_opszp0_8(G_nk+k) &
                          / (Opr_opsxp0_8(G_ni+ii) * Opr_opsyp0_8(G_nj+jj))
            Rhs(i,j,k) =Rhs(i,j,k)-stencil4*lhs(i,j-1,k)
          end do
!
          j=njl-sol_pil_n
          jj=j+l_j0-1
          do i=1+sol_pil_w, nil-sol_pil_e
             ii=i+l_i0-1
             stencil5=Opr_opsxp0_8(G_ni+ii) * Opr_opsyp2_8(2*G_nj+jj) * Opr_opszp0_8(G_nk+k) &
                           / (Opr_opsxp0_8(G_ni+ii) * Opr_opsyp0_8(G_nj+jj))

             Rhs(i,j,k) =Rhs(i,j,k)-stencil5*lhs(i,j+1,k)
          end do
!
      end do
!
!     ---------------------------------------------------------------
!
      return
      end

