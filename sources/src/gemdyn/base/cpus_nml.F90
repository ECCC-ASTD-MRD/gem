!**s/r cpus_nml - Read namelist &cpus

      integer function cpus_nml (F_unf)
      use lun
      use cpus_options
      implicit none
#include <arch_specific.hf>
!
      integer, intent(in) :: F_unf

      character(len=64) :: nml_S
      logical nml_must
!
!-------------------------------------------------------------------
!
! boiler plate - start
      if ( F_unf < 0 ) then
         cpus_nml= 0
         if ( Lun_out >= 0) write (Lun_out,nml=cpus)
         return
      end if

      cpus_nml= -1 ; nml_must= .false. ; nml_S= 'cpus'

      rewind(F_unf)
      read (F_unf, nml=cpus, end= 1001, err=1003)
      cpus_nml= 0 ; goto 1000
 1001 if (Lun_out >= 0) write (Lun_out, 6005) trim(nml_S)
      if (.not.nml_must) then
         cpus_nml= 1
         if (Lun_out >= 0) write (Lun_out, 6002) trim(nml_S)
      end if
      goto 1000
 1003 if (Lun_out >= 0) write (Lun_out, 6007) trim(nml_S)

 1000 if (cpus_nml < 0 ) return
      if ((Lun_out>=0).and.(cpus_nml==0)) write (Lun_out, 6004) trim(nml_S)
      cpus_nml= 1

 6002 format (' Skipping reading of namelist ',A)
 6004 format (' Reading of namelist ',A,' is successful')
 6005 format (' Namelist ',A,' NOT AVAILABLE')
 6007 format (/,' NAMELIST ',A,' IS INVALID'/)
! boiler plate - end
!
!-------------------------------------------------------------------
!
      return
      end function cpus_nml
