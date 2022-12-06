!-------------------------------------- LICENCE BEGIN ------------------------
!Environment Canada - Atmospheric Science and Technology License/Disclaimer,
!                     version 3; Last Modified: May 7, 2008.
!This is free but copyrighted software; you can use/redistribute/modify it under the terms
!of the Environment Canada - Atmospheric Science and Technology License/Disclaimer
!version 3 or (at your option) any later version that should be found at:
!http://collaboration.cmc.ec.gc.ca/science/rpn.comm/license.html
!
!This software is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
!without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!See the above mentioned License/Disclaimer for more details.
!You should have received a copy of the License/Disclaimer along with this software;
!if not, you can write to: EC-RPN COMM Group, 2121 TransCanada, suite 500, Dorval (Quebec),
!CANADA, H9P 1J3; or send e-mail to service.rpn@ec.gc.ca
!-------------------------------------- LICENCE END --------------------------
!/@*
subroutine litblrad (F_file_S, F_myproc)
   implicit none
!!!#include <arch_specific.hf>

      character(len=*) F_file_S
      integer F_myproc

!Author
!          B. Bilodeau (april 1994) - from lirg123
!
!Revision
!
! 001      B. Dugas (Aug 97) - Redefine IUNIT for FORTRAN file
! 002      M. Desgagne (Oct 98) - call back to rdradf_d (from dynamics)
! 003      B. Bilodeau (Jan 01) - remove call to ozpak
! 004      B. Bilodeau (May 03) - IBM conversion
!              - invert dimension of some radiation tables in order to
!                reduce the cache flooding in radir6
!
!Object
!          to read the radiation table from file (either unformatted
!          fortran binary file or RPN standard file) for infra-red
!          radiation calculation
!
!Arguments
!
!          - input -
! F_file_S  name of the radiation table file
!*@/

#include "radiation.cdk"
#include "ozopnt.cdk"

      external rd_radtab
      integer i,j,k,ij
      real, dimension(ntt,nco2)     :: bcninv, dbcninv
      real, dimension(mxx,nco2)     :: th2oinv
      real, dimension(mxx,ncx,nco2) :: yg3inv
!
!-----------------------------------------------------------------
!
!     calcul des pointeurs qui decoupent le champ g

      G1=5+1
      G2=G1+MXX*NTT
      G3=G2+MXX*NTT
      TH2O=G3+MXX*NTT
      TRO3=TH2O+MXX*NCO2
      YG3=TRO3+MXX
      BCN=YG3+NCO2*MXX*NCX
      DBCN=BCN+NCO2*NTT
      BO3=DBCN+NTT*NCO2
      DBO3=BO3+NTT
      TO3=DBO3+NTT
      UU=TO3+NO3
      TT=UU+MXX
!
      IF(TT+NTT-1 .NE. NTOTAL) THEN
         call physeterror('litblrad', 'Erreur dans les pointeurs.')
         return
      ENDIF

      call phyrdfile (F_file_S, rd_radtab, 'IRTAB', F_myproc)

!     inverser les dimensions de bcn, dbcn, th2o et yg3 pour
!     optimiser l'utilisation de la cache

      do j=1,ntt
         do i=1,nco2
            ij = (j-1)*nco2 + i - 1
            bcninv (j,i) = g(bcn +ij)
            dbcninv(j,i) = g(dbcn+ij)
         end do
      end do

      do j=1,mxx
         do i=1,nco2
            ij = (j-1)*nco2 + i - 1
            th2oinv(j,i) = g(th2o+ij)
         end do
      end do

      do k=1,ncx
         do j=1,mxx
            do i=1,nco2
               ij = (k-1)*mxx*nco2 + (j-1)*nco2 + i - 1
               yg3inv(j,k,i) = g(yg3+ij)
            end do
         end do
      end do

      do j=1,nco2
         do i=1,ntt
            ij=(j-1)*ntt + i - 1
            g(bcn +ij) = bcninv (i,j)
            g(dbcn+ij) = dbcninv(i,j)
         end do
      end do

      do j=1,nco2
         do i=1,mxx
            ij=(j-1)*mxx + i - 1
            g(th2o+ij) = th2oinv(i,j)
         end do
      end do

      do k=1,nco2
         do j=1,mxx
            do i=1,ncx
               ij=(k-1)*mxx*ncx + (j-1)*mxx + i - 1
               g(yg3+ij) = yg3inv(i,j,k)
            end do
         end do
      end do
!
!-----------------------------------------------------------------
!
      return
      end
