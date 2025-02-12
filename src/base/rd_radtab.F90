!**s/r rd_radtab  -- Perform the actual reading and organization of
!                    the data in the radiation file
module rd_radtab
   implicit none
   private
   public :: rd_radtab1
   
contains

   subroutine rd_radtab1(file, rbuf, dim, status)
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
!Object
!          Reads and organizes the data in cible according to
!          flag ozotable:   .true.  ===> for the ozone data
!                           .false. ===> for the radiation table
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
      integer :: code
!
!-----------------------------------------------------------------
!
      code   = status
      status = -1

      if (code.eq.200) then
         query = file%new_query(typvar="C ",nomvar="G1  ")
         success = query%find_next(record)
         if (.not.(( record%ni.eq.mxx .and. record%nj.eq.ntt) .and. success)) return
         dim(1) = 1
         dim(2) = ntotal
         status = 0
         return
      endif

      if ( (code.gt.200) .and. (code.le.300) )then
         success=file%read(record,data=c_loc(rbuf(g1)),nomvar='G1  ')
         success=file%read(record,data=c_loc(rbuf(g2)),nomvar='G2  ')
         success=file%read(record,data=c_loc(rbuf(g3)),nomvar='G3  ')
         success=file%read(record,data=c_loc(rbuf(th2o)),nomvar='2O  ')
         success=file%read(record,data=c_loc(rbuf(tro3)),nomvar='T3  ')
         success=file%read(record,data=c_loc(rbuf(yg3)),nomvar='Y3  ')
         success=file%read(record,data=c_loc(rbuf(bcn)),nomvar='BN  ')
         success=file%read(record,data=c_loc(rbuf(dbcn)),nomvar='DN  ')
         success=file%read(record,data=c_loc(rbuf(bo3)),nomvar='B3  ')
         success=file%read(record,data=c_loc(rbuf(dbo3)),nomvar='D3  ')
         success=file%read(record,data=c_loc(rbuf(to3)),nomvar='3O  ')
         success=file%read(record,data=c_loc(rbuf(uu)),nomvar='2U  ')
         success=file%read(record,data=c_loc(rbuf(tt)),nomvar='2T  ')
         status = 0
      endif

      if (code.ge.300) then
         g      = rbuf(1:ntotal)
         status = 0
      endif
!
!-----------------------------------------------------------------
!
      return
   end subroutine rd_radtab1

end module rd_radtab
