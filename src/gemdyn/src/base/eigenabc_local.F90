!---------------------------------- LICENCE BEGIN -------------------------------
! GEM - Library of kernel routines for the GEM numerical atmospheric model
! Copyright (C) 1990-2010 - Division de Recherche en Prevision Numerique
!                       Environnement Canada
! This library is free software; you can redistribute it and/or modify it
! under the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, version 2.1 of the License. This library is
! distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
! without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
! PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.
! You should have received a copy of the GNU Lesser General Public License
! along with this library; if not, write to the Free Software Foundation, Inc.,
! 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
!---------------------------------- LICENCE END ---------------------------------

!**s/r  Eigenabc_local - local eigenvalues and eigenvectors and inverse
!                      tridiagonal coefficients For the preconditioning
!
      subroutine eigenabc_local ( eval_local, evec_local, ai_local, &
                       bi_local,ci_local,nil,njl,Ni,Nj,Nk,i0,j0,wk)
      use HORgrid_options
      use gem_options
      use glb_ld
      use sol_mem
      use opr
      use, intrinsic :: iso_fortran_env
      implicit none

      integer nil,njl,Ni,Nj,Nk,i0,j0
      real(kind=REAL64)  eval_local(Ni),evec_local(Ni,Ni),ai_local(Ni,Nj,Nk), &
              bi_local(Ni,Nj,Nk),ci_local(Ni,Nj,Nk),wk(*)
!author
!     Abdessamad Qaddouri  - initial version - December 2006
!
      integer i,j,k,iloc,jloc,ii,jj,pnn,info
      real(kind=REAL64) a_8(Ni,Ni), b_8(Ni,Ni), d_8(3*Ni-1),&
                        r_8(Ni),di_8,cst,faz_8
      real(kind=REAL64), parameter :: zero=0.d0, one=1.d0
      real(kind=REAL64) a_81(Ni,Nj),b_81(Ni,Nj),c_81(Ni,Nj)
!
!     ---------------------------------------------------------------
!
      a_8=zero ; b_8=zero

      iloc=0
      do i = 1+sol_pil_w,nil-sol_pil_e-1
         ii=i+i0-1
         iloc=iloc+1
         a_8(iloc,iloc+1) = Opr_opsxp2_8(2*G_ni+ii)
         a_8(iloc,iloc  ) = Opr_opsxp2_8(G_ni+ii)
         a_8(iloc+1,iloc) = a_8(iloc,iloc+1)
         b_8(iloc,iloc+1) = Opr_opsxp0_8(2*G_ni+ii)
         b_8(iloc,iloc  ) = Opr_opsxp0_8(G_ni+ii)
         b_8(iloc+1,iloc) = b_8(iloc,iloc+1)
      end do
      a_8(Ni,Ni) = Opr_opsxp2_8(G_ni+(nil-sol_pil_e+i0-1))
      b_8(Ni,Ni) = Opr_opsxp0_8(G_ni+(nil-sol_pil_e+i0-1))

      pnn=ni
      call DSYGV ( 1, 'V', 'U', pnn, a_8, ni, b_8, ni, r_8,&
                                        d_8, 3*pnn-1, info )
      do j= 1, ni
         faz_8 = sign( one, a_8(1,j) )
         do i= 1, ni
            a_8(i,j) = faz_8 * a_8(i,j)
         end do
      end do

      evec_local = a_8
      eval_local = r_8

! inverse trid

      do k=1,Nk

         cst=wk(k)

         iloc=0
         do i = 1+sol_pil_w,nil-sol_pil_e
            iloc=iloc+1
            jloc=0
            do j=1+sol_pil_s,njl-sol_pil_n-1
               jj=j+j0-1
               jloc=jloc+1
               di_8= Opr_opsyp0_8(G_nj+jj) / (cos( G_yg_8 (jj) )**2)
               b_81(iloc,jloc)=eval_local(iloc) * di_8 + &
                     Opr_opsyp2_8(G_nj+jj)+cst*Opr_opsyp0_8(G_nj+jj)
               c_81(iloc,jloc)=Opr_opsyp2_8(2*G_nj+jj)
               a_81(iloc,jloc+1)=c_81(iloc,jloc)
            end do

            jj=njl-sol_pil_n+j0-1
            di_8= Opr_opsyp0_8(G_nj+jj) / (cos( G_yg_8 (jj)) **2)
            b_81(iloc,Nj)=eval_local(iloc)*di_8+Opr_opsyp2_8(G_nj+jj)+ &
                          cst*Opr_opsyp0_8(G_nj+jj)
            a_81(iloc,1  )= zero
            c_81(iloc,Nj) = zero
            if(Nj > 1) a_81(iloc,Nj) =c_81(iloc, Nj-1)
         end do

         ai_local(:,:,k)=zero;bi_local(:,:,k)=zero;ci_local(:,:,k)=zero

         do i=1,Ni
            bi_local(i,1,k)=b_81(i,1)
            ci_local(i,1,k)=c_81(i,1)
         end do

         do i=1,Ni
         do j=2,Nj
            ci_local(i,j,k)=c_81(i,j)
            ai_local(i,j,k)=a_81(i,j)/bi_local(i,j-1,k)
            bi_local(i,j,k)=b_81(i,j)-ai_local(i,j,k)*c_81(i,j-1)
         end do
         end do
      end do
!
!     ---------------------------------------------------------------
!
      return
      end
