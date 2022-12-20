!===================================================================
! Copyright: MSC-RPN COMM Group Licence/Disclaimer version 2
! http://www.cmc.ec.gc.ca/rpn/modcom/licence.html
!-------------------------------------------------------------------
! Modifications: [Date,Who,What]
! 2004-04, Stephane Chamberland
!    Original Code
! 2004-09, Stephane Chamberland
!    Make 'order' optional
!    New Param: SORT_UP, SORT_DOWN, SORT_FROMVALUES, SORT_FROMINDEX
!    Add: sort_string, sort_*_index, sort_*_multi
!    Remove: sort_int2tables, sort_real2tables, ASSENDING, DESCENDING
! 2005-04, Stephane Chamberland [debug (minor) w/ f90 on Pollux]
!    Solve embiguity for interfaces w/ *_multi s/r
!    Replace == by .eqv. for logical comparison
! 2012-01, Stephane Chamberland
!    Add real8 sort and unique option
!-------------------------------------------------------------------
! Dependencies:
! [Pending]
!-------------------------------------------------------------------
! Description
! [Pending]
!-------------------------------------------------------------------
! PRIVATE FN: (see below)
! PUBLIC FN: (see below)
! TODO: sort_string case independent option
! TODO: Add sorting for types int*8
! TODO: May modify sort_*_index to use indexes in table instead of sorting values [faster?]
!===================================================================
module sort_mod
   implicit none
   private
   public :: SORT_UP,SORT_DOWN,SORT_FROMVALUES,SORT_FROMINDEX,SORT_UNIQUE, &
        & sort,reverse,sort_get_indices
#include <rmn/msg.h>

   !---- Parameters -------------------------------------------------
   logical,parameter :: SORT_UP   = .true.
   logical,parameter :: SORT_DOWN = .false.
   logical,parameter :: SORT_FROMVALUES = .true.
   logical,parameter :: SORT_FROMINDEX  = .false.
   logical,parameter :: SORT_UNIQUE  = .true.
   integer,parameter :: MAX_STR_LEN = 2048

   !---- Interfaces -------------------------------------------------
   interface sort
      module procedure sort_int
      module procedure sort_int_2a
      module procedure sort_real
      module procedure sort_real8
      module procedure sort_string
!!$      module procedure sort_int_index
!!$      module procedure sort_real_index
!!$      module procedure sort_string_index
!!$      module procedure sort_int_multi
!!$      module procedure sort_real_multi
!!$      module procedure sort_string_multi
!!$      module procedure sort_int2tables
!!$      module procedure sort_real2tables
   end interface

   interface sort_get_indices
      module procedure sort_int_get_indices
   end interface

   interface reverse
      module procedure reverse_int
      module procedure reverse_real
      module procedure reverse_real8
      module procedure reverse_string
   end interface

   interface priv_unique
      module procedure priv_unique_int
      module procedure priv_unique_int_2a
      module procedure priv_unique_real
      module procedure priv_unique_real8
      module procedure priv_unique_string
   end interface

   integer,parameter :: MAXITEMS = 4096

contains


   !/@*
   function sort_int(table1,order_L,unique_L) result(dim)
      implicit none
      logical,intent(IN),optional :: order_L,unique_L
      integer,intent(INOUT) :: table1(:)
      integer :: dim
      !*@/
      integer :: jj,ii,tmp1
      !---------------------------------------------------------------
      dim = size(table1)
      do jj=2,dim
         tmp1=table1(jj)
         do ii=jj-1,1,-1
            if (table1(ii)<tmp1) goto 10
            table1(ii+1)=table1(ii)
         enddo
         ii=0
10       table1(ii+1)=tmp1
      enddo

      if (present(order_L)) then
         if (order_L.eqv.SORT_DOWN) call reverse(table1)
      endif

      if (present(unique_L)) then
         if (unique_L) call priv_unique(table1,dim)
      endif
      !---------------------------------------------------------------
      return
   end function sort_int


   !/@*
   function sort_int_2a(table1,table2,order_L,unique_L) result(dim)
      implicit none
      logical,intent(IN),optional :: order_L,unique_L
      integer,intent(INOUT) :: table1(:),table2(:)
      integer :: dim
      !*@/
      integer :: jj,ii,tmp1,tmp2
      !---------------------------------------------------------------
      dim = size(table1)
      do jj=2,dim
         tmp1=table1(jj)
         tmp2=table2(jj)
         do ii=jj-1,1,-1
            if (table1(ii)<tmp1) goto 10
            table1(ii+1)=table1(ii)
            table2(ii+1)=table2(ii)
         enddo
         ii=0
10       table1(ii+1)=tmp1
         table2(ii+1)=tmp2
      enddo

      if (present(order_L)) then
         if (order_L.eqv.SORT_DOWN) call reverse(table1)
         if (order_L.eqv.SORT_DOWN) call reverse(table2)
      endif

      if (present(unique_L)) then
         if (unique_L) call priv_unique(table1,table2,dim)
      endif
      !---------------------------------------------------------------
      return
   end function sort_int_2a


   !/@*
   function sort_int_get_indices(indexes,table1,order_L,unique_L) result(dim)
      implicit none
      logical,intent(IN),optional :: order_L,unique_L
      integer,intent(IN)  :: table1(:)
      integer,intent(OUT) :: indexes(:)
      integer :: dim
      !*@/
      integer :: ii,table1copy(MAXITEMS)
      logical :: order1_L,unique1_L
      !---------------------------------------------------------------
      dim = size(table1)
      if (size(indexes) < dim .or. dim > MAXITEMS) then
         call msg(MSG_ERROR,'(sort) Dataset too big')
         indexes = 1
         dim = -1
         return
      endif
      order1_L = SORT_UP
      unique1_L = .false.
      if (present(order_L)) order1_L = order_L
      if (present(unique_L)) unique1_L = unique_L

      do ii=1,dim
         indexes = ii
      enddo
      table1copy(1:dim) = table1
      dim = sort_int_2a(table1copy(1:dim),indexes(1:dim),order1_L,unique1_L)
      !---------------------------------------------------------------
      return
   end function sort_int_get_indices


   !/@*
   function sort_real(table1,order_L,unique_L) result(dim)
      implicit none
      logical,intent(IN),optional :: order_L,unique_L
      real,intent(INOUT) :: table1(:)
      integer :: dim
      !*@/
      integer :: jj,ii
      real    :: tmp1
      !---------------------------------------------------------------
      dim = size(table1)
      do jj=2,dim
         tmp1=table1(jj)
         do ii=jj-1,1,-1
            if (table1(ii)<tmp1) goto 10
            table1(ii+1)=table1(ii)
         enddo
         ii=0
10       table1(ii+1)=tmp1
      enddo

      if (present(order_L)) then
         if (order_L.eqv.SORT_DOWN) call reverse(table1)
      endif

      if (present(unique_L)) then
         if (unique_L) call priv_unique(table1,dim)
      endif
      !---------------------------------------------------------------
      return
   end function sort_real


   !/@*
   function sort_real8(table1,order_L,unique_L) result(dim)
      implicit none
      logical,intent(IN),optional :: order_L,unique_L
      real(8),intent(INOUT) :: table1(:)
      integer :: dim
      !*@/
      integer :: jj,ii
      real(8) :: tmp1
      !---------------------------------------------------------------
      dim = size(table1)
      do jj=2,dim
         tmp1=table1(jj)
         do ii=jj-1,1,-1
            if (table1(ii)<tmp1) goto 10
            table1(ii+1)=table1(ii)
         enddo
         ii=0
10       table1(ii+1)=tmp1
      enddo

      if (present(order_L)) then
         if (order_L.eqv.SORT_DOWN) call reverse(table1)
      endif

      if (present(unique_L)) then
         if (unique_L) call priv_unique(table1,dim)
      endif
      !---------------------------------------------------------------
     return
   end function sort_real8


   !/@*
   function sort_string(table1,str_len0,order_L,unique_L) result(dim)
      implicit none
      logical,intent(IN),optional :: order_L,unique_L
      integer,intent(IN),optional :: str_len0 !sort using only first str_len char
      character(len=*),intent(INOUT) :: table1(:)
      integer :: dim
      !*@/
      integer jj,ii,str_len
      character(len=MAX_STR_LEN) :: tmp1,tmp2
      !---------------------------------------------------------------
      str_len = min(len(table1(1)),len(tmp1))
      if (present(str_len0)) str_len = min(str_len0,str_len)
      dim = size(table1)
      do jj=2,dim
         tmp1=table1(jj)
         do ii=jj-1,1,-1
            tmp2=table1(ii)
            if (llt(tmp2(1:str_len),tmp1(1:str_len))) goto 10
            table1(ii+1)=table1(ii)
         enddo
         ii=0
10       table1(ii+1)=tmp1
      enddo

      if (present(order_L)) then
         if (order_L.eqv.SORT_DOWN) call reverse(table1)
      endif

      if (present(unique_L)) then
         if (unique_L) call priv_unique(table1,dim,str_len)
      endif
      !---------------------------------------------------------------
      return
   end function sort_string


   !/@*
   subroutine reverse_int(table1)
      implicit none
      integer,intent(INOUT) :: table1(:)
      !*@/
      integer :: ii,tmp,dim
      !---------------------------------------------------------------
      dim = size(table1)
      do ii=1,dim/2
         tmp=table1(ii)
         table1(ii)=table1(dim-ii+1)
         table1(dim-ii+1)=tmp
      enddo
      !---------------------------------------------------------------
      return
   end subroutine reverse_int

   !/@*
   subroutine reverse_real(table1)
      implicit none
      real,intent(INOUT) :: table1(:)
      !*@/
      integer :: ii,dim
      real    :: tmp
      !---------------------------------------------------------------
      dim = size(table1)
      do ii=1,dim/2
         tmp=table1(ii)
         table1(ii)=table1(dim-ii+1)
         table1(dim-ii+1)=tmp
      enddo
      !---------------------------------------------------------------
      return
   end subroutine reverse_real


   !/@*
   subroutine reverse_real8(table1)
      implicit none
      real(8),intent(INOUT) :: table1(:)
      !*@/
      integer :: ii,dim
      real(8) :: tmp
      !---------------------------------------------------------------
      dim = size(table1)
      do ii=1,dim/2
         tmp=table1(ii)
         table1(ii)=table1(dim-ii+1)
         table1(dim-ii+1)=tmp
      enddo
      !---------------------------------------------------------------
      return
   end subroutine reverse_real8


   !/@*
   subroutine reverse_string(table1)
      implicit none
      character(len=*),intent(INOUT) :: table1(:)
      !*@/
      integer :: ii,dim
      character(len=MAX_STR_LEN) :: tmp
      !---------------------------------------------------------------
      dim = size(table1)
      do ii=1,dim/2
         tmp=table1(ii)
         table1(ii)=table1(dim-ii+1)
         table1(dim-ii+1)=tmp
      enddo
      !---------------------------------------------------------------
      return
   end subroutine reverse_string

   !TODO: review/test the following sort s/r

!!$   !=================================================================
!!$   subroutine sort_int_index(itable,vtable,dim,stype,order)
!!$      implicit none
!!$      logical,intent(IN),optional :: order
!!$      logical,intent(IN) :: stype
!!$      integer,intent(IN) :: dim
!!$      integer,dimension(dim),intent(INOUT) :: itable,vtable
!!$
!!$      !----
!!$      integer :: jj,ii,tmpi,tmpv
!!$      integer,dimension(dim) :: itable2
!!$      !---------------------------------------------------------------
!!$
!!$      if (stype.eqv.SORT_FROMVALUES) then
!!$         !---- Sort according to values, create indexes table
!!$         do ii=1,dim
!!$            itable(ii)=ii
!!$         enddo
!!$
!!$         do jj=2,dim
!!$            tmpi=itable(jj)
!!$            tmpv=vtable(jj)
!!$            do ii=jj-1,1,-1
!!$               if (vtable(ii)<tmpv) goto 10
!!$               itable(ii+1)=itable(ii)
!!$               vtable(ii+1)=vtable(ii)
!!$            enddo
!!$            ii=0
!!$10          itable(ii+1)=tmpi
!!$            vtable(ii+1)=tmpv
!!$         enddo
!!$
!!$         if (present(order)) then
!!$            if (order.eqv.SORT_DOWN) then
!!$               call reverse(itable,dim)
!!$               call reverse(vtable,dim)
!!$            endif
!!$         endif
!!$      else
!!$         !---- Sort values according to indexes table
!!$         !- Make a copy of index table to preserve it
!!$         do ii=1,dim
!!$            itable2(ii)=itable(ii)
!!$         enddo
!!$
!!$         do jj=2,dim
!!$            tmpi=itable2(jj)
!!$            tmpv=vtable(jj)
!!$            do ii=jj-1,1,-1
!!$               if (itable2(ii)<tmpi) goto 20
!!$               itable2(ii+1)=itable2(ii)
!!$               vtable(ii+1)=vtable(ii)
!!$            enddo
!!$            ii=0
!!$20          itable2(ii+1)=tmpi
!!$            vtable(ii+1)=tmpv
!!$         enddo
!!$      endif
!!$
!!$   end subroutine sort_int_index
!!$
!!$
!!$   !=================================================================
!!$   subroutine sort_real_index(itable,vtable,dim,stype,order)
!!$      implicit none
!!$      logical,intent(IN),optional :: order
!!$      logical,intent(IN) :: stype
!!$      integer,intent(IN) :: dim
!!$      integer,dimension(dim),intent(INOUT) :: itable
!!$      real,dimension(dim),intent(INOUT) :: vtable
!!$
!!$      !----
!!$      integer :: jj,ii,tmpi
!!$      real :: tmpv
!!$      integer,dimension(dim) :: itable2
!!$      !---------------------------------------------------------------
!!$
!!$      if (stype.eqv.SORT_FROMVALUES) then
!!$         !---- Sort according to values, create indexes table
!!$         do ii=1,dim
!!$            itable(ii)=ii
!!$         enddo
!!$
!!$         do jj=2,dim
!!$            tmpi=itable(jj)
!!$            tmpv=vtable(jj)
!!$            do ii=jj-1,1,-1
!!$               if (vtable(ii)<tmpv) goto 10
!!$               itable(ii+1)=itable(ii)
!!$               vtable(ii+1)=vtable(ii)
!!$            enddo
!!$            ii=0
!!$10          itable(ii+1)=tmpi
!!$            vtable(ii+1)=tmpv
!!$         enddo
!!$
!!$         if (present(order)) then
!!$            if (order.eqv.SORT_DOWN) then
!!$               call reverse(itable,dim)
!!$               call reverse(vtable,dim)
!!$            endif
!!$         endif
!!$      else
!!$         !---- Sort values according to indexes table
!!$         !- Make a copy of index table to preserve it
!!$         do ii=1,dim
!!$            itable2(ii)=itable(ii)
!!$         enddo
!!$
!!$         do jj=2,dim
!!$            tmpi=itable2(jj)
!!$            tmpv=vtable(jj)
!!$            do ii=jj-1,1,-1
!!$               if (itable2(ii)<tmpi) goto 20
!!$               itable2(ii+1)=itable2(ii)
!!$               vtable(ii+1)=vtable(ii)
!!$            enddo
!!$            ii=0
!!$20          itable2(ii+1)=tmpi
!!$            vtable(ii+1)=tmpv
!!$         enddo
!!$      endif
!!$
!!$   end subroutine sort_real_index
!!$
!!$   !=================================================================
!!$   subroutine sort_string_index(itable,vtable,dim,str_len,stype,order)
!!$      implicit none
!!$      logical,intent(IN),optional :: order
!!$      logical,intent(IN) :: stype
!!$      integer,intent(IN) :: dim,str_len
!!$      integer,dimension(dim),intent(INOUT) :: itable
!!$      character(len=*),dimension(dim),intent(INOUT) :: vtable
!!$
!!$      !----
!!$      integer :: jj,ii,tmpi
!!$      character(len=str_len) :: tmpv
!!$      integer,dimension(dim) :: itable2
!!$      !---------------------------------------------------------------
!!$
!!$      if (stype.eqv.SORT_FROMVALUES) then
!!$         !---- Sort according to values, create indexes table
!!$         do ii=1,dim
!!$            itable(ii)=ii
!!$         enddo
!!$
!!$         do jj=2,dim
!!$            tmpi=itable(jj)
!!$            tmpv=vtable(jj)
!!$            do ii=jj-1,1,-1
!!$               if (vtable(ii)<tmpv) goto 10
!!$               itable(ii+1)=itable(ii)
!!$               vtable(ii+1)=vtable(ii)
!!$            enddo
!!$            ii=0
!!$10          itable(ii+1)=tmpi
!!$            vtable(ii+1)=tmpv
!!$         enddo
!!$
!!$         if (present(order)) then
!!$            if (order.eqv.SORT_DOWN) then
!!$               call reverse(itable,dim)
!!$               call reverse(vtable,dim,str_len)
!!$            endif
!!$         endif
!!$      else
!!$         !---- Sort values according to indexes table
!!$         !- Make a copy of index table to preserve it
!!$         do ii=1,dim
!!$            itable2(ii)=itable(ii)
!!$         enddo
!!$
!!$         do jj=2,dim
!!$            tmpi=itable2(jj)
!!$            tmpv=vtable(jj)
!!$            do ii=jj-1,1,-1
!!$               if (itable2(ii)<tmpi) goto 20
!!$               itable2(ii+1)=itable2(ii)
!!$               vtable(ii+1)=vtable(ii)
!!$            enddo
!!$            ii=0
!!$20          itable2(ii+1)=tmpi
!!$            vtable(ii+1)=tmpv
!!$         enddo
!!$      endif
!!$
!!$   end subroutine sort_string_index
!!$
!!$   !=================================================================
!!$   subroutine sort_int_multi(table1,dim,order,&
!!$        & table2,table3,table4,table5,table6,table7,table8,table9)
!!$      implicit none
!!$      logical,intent(IN) :: order
!!$      integer,intent(IN) :: dim
!!$      integer,dimension(dim),intent(INOUT) :: table1,table2
!!$      integer,dimension(dim),intent(INOUT),optional :: &
!!$           & table3,table4,table5,table6,table7, &
!!$           & table8,table9
!!$      !----
!!$      integer,dimension(dim) :: itable
!!$      !---------------------------------------------------------------
!!$      call sort_int_index(itable,table1,dim,SORT_FROMVALUES,order)
!!$      call sort_int_index(itable,table2,dim,SORT_FROMINDEX,SORT_UP)
!!$      if (present(table3)) then
!!$         call sort_int_index(itable,table3,dim,SORT_FROMINDEX,SORT_UP)
!!$      endif
!!$      if (present(table4)) then
!!$         call sort_int_index(itable,table4,dim,SORT_FROMINDEX,SORT_UP)
!!$      endif
!!$      if (present(table5)) then
!!$         call sort_int_index(itable,table5,dim,SORT_FROMINDEX,SORT_UP)
!!$      endif
!!$      if (present(table6)) then
!!$         call sort_int_index(itable,table6,dim,SORT_FROMINDEX,SORT_UP)
!!$      endif
!!$      if (present(table7)) then
!!$         call sort_int_index(itable,table7,dim,SORT_FROMINDEX,SORT_UP)
!!$      endif
!!$      if (present(table8)) then
!!$         call sort_int_index(itable,table8,dim,SORT_FROMINDEX,SORT_UP)
!!$      endif
!!$      if (present(table9)) then
!!$         call sort_int_index(itable,table9,dim,SORT_FROMINDEX,SORT_UP)
!!$      endif
!!$   end subroutine sort_int_multi
!!$
!!$   !=================================================================
!!$   subroutine sort_real_multi(table1,dim,order,&
!!$        & table2,table3,table4,table5,table6,table7,table8,table9)
!!$      implicit none
!!$      logical,intent(IN) :: order
!!$      integer,intent(IN) :: dim
!!$      real,dimension(dim),intent(INOUT) :: table1,table2
!!$      real,dimension(dim),intent(INOUT),optional :: &
!!$           & table3,table4,table5,table6,table7, &
!!$           & table8,table9
!!$      !----
!!$      integer,dimension(dim) :: itable
!!$      !---------------------------------------------------------------
!!$      call sort_real_index(itable,table1,dim,SORT_FROMVALUES,order)
!!$      call sort_real_index(itable,table2,dim,SORT_FROMINDEX,SORT_UP)
!!$      if (present(table3)) then
!!$         call sort_real_index(itable,table3,dim,SORT_FROMINDEX,SORT_UP)
!!$      endif
!!$      if (present(table4)) then
!!$         call sort_real_index(itable,table4,dim,SORT_FROMINDEX,SORT_UP)
!!$      endif
!!$      if (present(table5)) then
!!$         call sort_real_index(itable,table5,dim,SORT_FROMINDEX,SORT_UP)
!!$      endif
!!$      if (present(table6)) then
!!$         call sort_real_index(itable,table6,dim,SORT_FROMINDEX,SORT_UP)
!!$      endif
!!$      if (present(table7)) then
!!$         call sort_real_index(itable,table7,dim,SORT_FROMINDEX,SORT_UP)
!!$      endif
!!$      if (present(table8)) then
!!$         call sort_real_index(itable,table8,dim,SORT_FROMINDEX,SORT_UP)
!!$      endif
!!$      if (present(table9)) then
!!$         call sort_real_index(itable,table9,dim,SORT_FROMINDEX,SORT_UP)
!!$      endif
!!$   end subroutine sort_real_multi
!!$
!!$   !=================================================================
!!$   subroutine sort_string_multi(table1,dim,str_len,order,&
!!$        & table2,table3,table4,table5,table6,table7,table8,table9)
!!$      implicit none
!!$      logical,intent(IN) :: order
!!$      integer,intent(IN) :: dim,str_len
!!$      character(len=*),dimension(dim),intent(INOUT) :: table1,table2
!!$      character(len=*),dimension(dim),intent(INOUT),optional :: &
!!$           & table3,table4,table5,table6,table7, &
!!$           & table8,table9
!!$      !----
!!$      integer,dimension(dim) :: itable
!!$      !---------------------------------------------------------------
!!$      call sort_string_index(itable,table1,dim,str_len,SORT_FROMVALUES,order)
!!$      call sort_string_index(itable,table2,dim,str_len,SORT_FROMINDEX,SORT_UP)
!!$      if (present(table3)) then
!!$         call sort_string_index(itable,table3,dim,str_len,SORT_FROMINDEX,SORT_UP)
!!$      endif
!!$      if (present(table4)) then
!!$         call sort_string_index(itable,table4,dim,str_len,SORT_FROMINDEX,SORT_UP)
!!$      endif
!!$      if (present(table5)) then
!!$         call sort_string_index(itable,table5,dim,str_len,SORT_FROMINDEX,SORT_UP)
!!$      endif
!!$      if (present(table6)) then
!!$         call sort_string_index(itable,table6,dim,str_len,SORT_FROMINDEX,SORT_UP)
!!$      endif
!!$      if (present(table7)) then
!!$         call sort_string_index(itable,table7,dim,str_len,SORT_FROMINDEX,SORT_UP)
!!$      endif
!!$      if (present(table8)) then
!!$         call sort_string_index(itable,table8,dim,str_len,SORT_FROMINDEX,SORT_UP)
!!$      endif
!!$      if (present(table9)) then
!!$         call sort_string_index(itable,table9,dim,str_len,SORT_FROMINDEX,SORT_UP)
!!$      endif
!!$   end subroutine sort_string_multi


!!$  !=================================================================
!!$  !
!!$  !=================================================================
!!$  SUBROUTINE sort_int2tables(table1,table2,dim,order)
!!$    IMPLICIT NONE
!!$    LOGICAL,INTENT(IN) :: order
!!$    INTEGER,INTENT(IN) :: dim
!!$    INTEGER,DIMENSION(dim),INTENT(INOUT) :: table1,table2
!!$
!!$    !----
!!$    INTEGER jj,ii,tmp1,tmp2
!!$    !---------------------------------------------------------------
!!$    DO jj=2,dim
!!$       tmp1=table1(jj)
!!$       tmp2=table2(jj)
!!$       DO ii=jj-1,1,-1
!!$          IF (table1(ii)<tmp1) GOTO 10
!!$          table1(ii+1)=table1(ii)
!!$          table2(ii+1)=table2(ii)
!!$       ENDDO
!!$       ii=0
!!$10     table1(ii+1)=tmp1
!!$       table2(ii+1)=tmp2
!!$    ENDDO
!!$
!!$    IF (.NOT.order) THEN
!!$       CALL reverse(table1,dim)
!!$       CALL reverse(table2,dim)
!!$    ENDIF
!!$
!!$  END SUBROUTINE sort_int2tables
!!$
!!$  !=================================================================
!!$  !
!!$  !=================================================================
!!$  SUBROUTINE sort_real2tables(table1,table2,dim,order)
!!$    IMPLICIT NONE
!!$    LOGICAL,INTENT(IN) :: order
!!$    INTEGER,INTENT(IN) :: dim
!!$    REAL,DIMENSION(dim),INTENT(INOUT) :: table1,table2
!!$
!!$    !----
!!$    INTEGER jj,ii
!!$    REAL    tmp1,tmp2
!!$    !---------------------------------------------------------------
!!$    DO jj=2,dim
!!$       tmp1=table1(jj)
!!$       tmp2=table2(jj)
!!$       DO ii=jj-1,1,-1
!!$          IF (table1(ii)<tmp1) GOTO 10
!!$          table1(ii+1)=table1(ii)
!!$          table2(ii+1)=table2(ii)
!!$       ENDDO
!!$       ii=0
!!$10     table1(ii+1)=tmp1
!!$       table2(ii+1)=tmp2
!!$    ENDDO
!!$
!!$    IF (.NOT.order) THEN
!!$       CALL reverse(table1,dim)
!!$       CALL reverse(table2,dim)
!!$    ENDIF
!!$
!!$  END SUBROUTINE sort_real2tables


   !=================================================================
   !
   !=================================================================
   subroutine priv_unique_int(table1,dim)
      implicit none
      integer,intent(INOUT) :: dim
      integer,dimension(dim),intent(INOUT) :: table1
      !----
      integer :: ii,ii2
      !---------------------------------------------------------------
      ii2=1
      do ii=2,dim
         if (table1(ii) /= table1(ii-1)) ii2 = ii2+1 
         table1(ii2) = table1(ii)
      enddo
      dim = ii2
      return
   end subroutine priv_unique_int

   subroutine priv_unique_int_2a(table1,table2,dim)
      implicit none
      integer,intent(INOUT) :: dim
      integer,dimension(dim),intent(INOUT) :: table1,table2(:)
      !----
      integer :: ii,ii2
      !---------------------------------------------------------------
      ii2=1
      do ii=2,dim
         if (table1(ii) /= table1(ii-1)) ii2 = ii2+1 
         table1(ii2) = table1(ii)
         table2(ii2) = table2(ii)
      enddo
      dim = ii2
      return
   end subroutine priv_unique_int_2a

   !=================================================================
   subroutine priv_unique_real(table1,dim)
      implicit none
      integer,intent(INOUT) :: dim
      real,dimension(dim),intent(INOUT) :: table1
      !----
      integer :: ii,ii2
      !---------------------------------------------------------------
      ii2=1
      do ii=2,dim
         if (table1(ii) /= table1(ii-1)) ii2 = ii2+1 
         table1(ii2) = table1(ii)
      enddo
      dim = ii2
      return
   end subroutine priv_unique_real

   !=================================================================
   subroutine priv_unique_real8(table1,dim)
      implicit none
      integer,intent(INOUT) :: dim
      real(8),dimension(dim),intent(INOUT) :: table1
      !----
      integer :: ii,ii2
      !---------------------------------------------------------------
      ii2=1
      do ii=2,dim
         if (table1(ii) /= table1(ii-1)) ii2 = ii2+1 
         table1(ii2) = table1(ii)
      enddo
      dim = ii2
      return
   end subroutine priv_unique_real8

   !=================================================================
   subroutine priv_unique_string(table1,dim,str_len)
      implicit none
      integer,intent(IN) :: str_len
      integer,intent(INOUT) :: dim
      character(len=*),dimension(dim),intent(INOUT) :: table1
      !----
      integer :: ii,ii2
      !---------------------------------------------------------------
      ii2=1
      do ii=2,dim
         if (table1(ii)(1:str_len) /= table1(ii-1)(1:str_len)) ii2 = ii2+1 
         table1(ii2) = table1(ii)
      enddo
      dim = ii2
      return
   end subroutine priv_unique_string


end module sort_mod
