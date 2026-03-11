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

!**s/r INs_gz3d - Obtain 3D GZ directly from input file or
!                 from VGD descriptor when input kind=21

      logical function INs_gz3d ()
      use INs_base
      implicit none

      logical, external :: INs_gz_from_vgd
      
      character(len=1) :: grid_S(3)
      character(len=6) :: vname,dumc
      integer :: version,lislon,deb,err
      integer :: ip1_list(2*(INs_nka+2))
      real :: level!, wk(INs_nid,INs_njd,200)
      real(kind=REAL64) :: pref_a_8
!!$ 0   KIND_ABOVE_SEA  height (m) above mean sea level         (-20,000 -> 100,000)
!!$ 1   KIND_SIGMA      sigma coordinates       (0.0 -> 1.0)
!!$ 2   KIND_PRESSURE   pressure (mb)   (0 -> 1100)
!!$ 3   KIND_ARBITRARY  arbitrary number, no units      (-4.8e8 -> 1.0e10)
!!$ 4   KIND_ABOVE_GND  height (m) above ground         (-20,000 -> 100,000)
!!$ 5   KIND_HYBRID     hybrid coordinates      (0.0 -> 1.0)
!!$ 6   KIND_THETA      theta coordinates       (1 -> 200,000)
!!$10   KIND_HOURS      time (hours)    (0.0 -> 1.0e10)
!!$15   KIND_SAMPLES    reserved (integer value)        (0 -> 1 999 999)
!!$17   KIND_MTX_IND    conversion matrix x subscript
!!$                     (shared with kind=1) (1.0 -> 1.0e10)
!!$21   KIND_M_PRES     pressure-meters (shared with kind=5)
!!$                     (0 -> 1,000,000) fact=1E+4 
!
!-----------------------------------------------------------------------
!
      if (INs_1o1_L) err = fstopc('MSGLVL','INFORM',RMN_OPT_SET)
      INs_gz3d= .false.
!  GZcaract(1:5)     
!!$ 1 - nk gz found
!!$ 2 - Inp_kind de GZ
!!$ 3 - me_L     
!!$ 4 - mels_L     
!!$ 5 - pref_a
      GZcaract(:)= 0

      grid_S=['Q','U','V']
      
      call INs_getvgd ()

      err = INs_read_mt ('GEOPOTENTIAL',grid_S,vname,GZ,3,GZIP1,&
                      lislon,Rot_ig1, Rot_ig2, Rot_ig3, Rot_ig4,&
                      F_quiet_L=.true.)

      if ((err == 0) .and. (lislon>1)) then ! Found GZ in input file
                       
         call convip (GZIP1(lislon),level,Inp_kind,-1,dumc,.false.)
         GZcaract(1)= lislon
         GZcaract(2)= Inp_kind
         GZcaract(3)= 1
         if (Lun_out > 0) write(lun_out,9100) lislon,&
                            'GEOPOTENTIAL',GZcaract(2)
         deb = 3*lislon+1
         err = INs_read_mt ('MELS','Q', vname,GZ(1:,1:,deb),1,ip1_list,&
                             lislon,Rot_ig1, Rot_ig2, Rot_ig3, Rot_ig4,&
                             F_quiet_L=.true.)
         if ((err == 0) .and. (lislon>0)) then
            GZcaract(4)= 1
         endif
         INs_gz3d= .true.
         if (Inp_kind == 21) return
         err= vgd_get ( Inp_vgd_src, key='VERS',value=version )
         if ( (Inp_kind == 5) .or. &
             ((Inp_kind == 1).and.(version == 3)) ) then
             err= vgd_get ( Inp_vgd_src, key='PREF',value=pref_a_8 )
             GZcaract(5)= pref_a_8
             if (err/=0) INs_gz3d= .false.
         end if
      else
         if (Inp_vgdkind==21) INs_gz3d= INs_gz_from_vgd ()
      endif
                             
      if (INs_1o1_L) err = fstopc ('MSGLVL','WARNIN',RMN_OPT_SET)

 9100 format(1x,'INs_data: ', i3,' levels of ',a,&
             ' found in input file with kind=',i5)
!
!-----------------------------------------------------------------------
!
      return
      end function INs_gz3d
