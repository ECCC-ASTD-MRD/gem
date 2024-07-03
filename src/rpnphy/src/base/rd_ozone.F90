!-------------------------------------- LICENCE BEGIN ------------------------------------
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
!-------------------------------------- LICENCE END --------------------------------------
!**s/r rd_ozone  -- Perform the actual reading and organization of
!                   the data in the ozone file
!

module rd_ozone
   implicit none
   private
   public :: rd_ozone1
      
contains

   subroutine rd_ozone1(file, rbuf, dim, status)
      use rmn_fst24
      implicit none
!!!#include <arch_specific.hf>
!
      type(fst_file), intent(inout) :: file
      real, intent(inout), target :: rbuf(*)
      integer, intent(inout) :: dim(:)
      integer, intent(inout) :: status
!
!Author
!          M. Desgagne (Spring 2008)
!
!Arguments
!          - Input -
! file     RPS STD file object
!          - Input/Output -
! rbuf     read buffer
!          - Output -
! status   exit status for the routine (=0 if OK,  =-1 otherwise)
!
#include "radiation.cdk"
#include "ozopnt.cdk"
!
      type(fst_record) :: record
      type(fst_query)  :: query
      logical          :: success
      integer :: i,m,NLP,code
!
!-----------------------------------------------------------------
!
      code   = status
      status = -1

      if (code.eq.200) then
         query = file%new_query(nomvar='O3  ')
         success = query%find_next(record)
         if (.not. query%find_next(record)) return

         NLACL=record%ni
         NPCL=record%nj
         dim(1) = 3
         dim(2) = NLACL*NPCL*13 + NLACL + NPCL
         dim(3) = NLACL
         dim(4) = NPCL
         status = 0
         return
      endif

      NLP = NLACL*NPCL

      if ( (code.gt.200) .and. (code.le.300) )then
         if (.not. file%read(record,data=c_loc(rbuf),nomvar='ZLAT')) return
         if (.not. file%read(record,data=c_loc(rbuf(NLACL+1)),nomvar='PREF')) return
         do m=1,12
            success = file%read(record,data=c_loc(rbuf(NLACL+NPCL+(m-1)*NLP+1)),ip3=m,nomvar='O3  ')
         end do
         status = 0
      endif

      if (code.ge.300) then
         NLACL = dim(3)
         NPCL  = dim(4)
         NLP   = NLACL*NPCL

         allocate (gozon12(NLP,12), goz(NLP+NLACL+NPCL))

         DO i=1,NLACL
            goz(NLP      +i)=  rbuf(i)
         ENDDO

         DO i=1,NPCL
            goz(NLP+NLACL+i)= rbuf(NLACL+i)
         ENDDO
         do m=1,12
         do i=1,NLP
            gozon12(i,m) = rbuf(NLACL+NPCL+(m-1)*NLP+i)
         end do
         end do
         status = 0
      endif
!
!-----------------------------------------------------------------
!
      return
   end subroutine rd_ozone1

end module rd_ozone
