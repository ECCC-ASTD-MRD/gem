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

module INs_base
      use iso_c_binding
      use, intrinsic :: iso_fortran_env
      use IOs
      use INs
      use tdpack
      use rmn_fst24
      implicit none
#include <rmnlib_basics.hf>
      
contains

!**   s/r INs_read_mt - Read variable F_var_S and perform horizontal
!                       interpolation to F_nd Arakawa grid destinations

      integer function INs_read_mt ( F_var_S, F_hgrid_S, F_nomvar_S,&
             F_dest, F_nd, F_ip1, F_nka, F_ig1, F_ig2, F_ig3, F_ig4,&
             F_hint_S, F_quiet_L, F_type_S )
      implicit none
      character(len=*)          ,intent(in)  :: F_var_S
      character(len=*), dimension(*),intent(in) :: F_hgrid_S
      character(len=6), intent(OUT)          :: F_nomvar_S
      character(len=*), optional,intent(in)  :: F_hint_S
      logical         , optional,intent(in)  :: F_quiet_L
      character(len=1), optional,intent(in)  :: F_type_S
      integer                   ,intent(in ) :: F_nd,F_ig1,F_ig2,F_ig3,F_ig4
      integer                   ,intent(out) :: F_nka
      integer, intent(OUT) :: F_ip1(*)
      real, intent(OUT) :: F_dest(G_ni+2*G_halox,G_nj+2*G_haloy,*)

      integer, external :: samegrid_gid, samegrid_rot
      character(len=1) typ,grd
      character(len=4) nomvar,var,dumc
      character(len=12) lab,interp_S
      logical :: quiet_L
      integer, parameter :: nlis = 1024
      integer i, k, idst, err, nz, n1,n2,n3, nrec, liste(nlis),&
              liste_sorted(nlis),lislon,maxdim_wk2
      integer subid,nicore,njcore,datev,ni_dest,nj_dest
      integer mpx,local_nk,irest,kstart, dstf_gid,src_gid, ip1
      integer dte, det, ipas, p1, p2, p3, g1, g2, g3, g4, bit, &
              dty, swa, lng, dlf, ubc, ex1, ex2, ex3
      real :: surface_level
      real, dimension(:  ), allocatable, target :: wk1
      real, dimension(:  ), pointer     :: posx,posy
      real(kind=REAL64) add, mult

      type(fst_query)  :: query
      type(fst_record) :: recs(nlis) 
      logical          :: success
!
!---------------------------------------------------------------------
!
      INs_read_mt= -1
      F_nka= -1 ; local_nk= 0
      add= 0.d0 ; mult= 1.d0
      quiet_L=.false.
      if (present(F_quiet_L)) quiet_L= F_quiet_L
      typ= ' '
      if (present(F_type_S)) typ= F_type_S
      
      nomvar = F_var_S ; ip1= -1
      select case (F_var_S)
         case ('OROGRAPHY')
            if (Inp_kind == 2  ) then
               nomvar= '@NUL'
!!$               if (Inp_src_PX_L) then
!!$                  nomvar= 'GZ' ; surface_level= 1. ; p1=5
!!$                  call convip ( ip1, surface_level,p1,1,dumc,.false. )
!!$               endif
            endif
!!$            if (Inp_kind == 1 ) then
!!$               nomvar= 'GZ' ; ip1= 12000
!!$               if (Inp_src_PX_L) then
!!$                  nomvar= 'GZ' ; surface_level= 1. ; p1=5
!!$                  call convip ( ip1, surface_level,p1,1,dumc,.false. )
!!$               endif
!!$            endif
            if (Inp_kind == 5 ) then
               nomvar= 'GZ' ; surface_level= 1.
               call convip ( ip1, surface_level,Inp_kind,1,dumc,.false. )
            endif
            if ( Inp_src_hauteur_L ) then
              ! nomvar= 'GZ'
               nomvar= 'ME'
               ip1=0
              ! if (Inp_kind==21) surface_level= 0.
              ! if (Inp_kind==5 ) surface_level= 1.
              ! call convip ( ip1, surface_level,Inp_kind,1,dumc,.false. )
            endif
            if ( nomvar == 'GZ' ) mult= 10.d0 * grav_8
         case ('SFCPRES')
            nomvar= 'P0'
            if (Inp_kind == 2  ) nomvar= '@NUL'
            !if (Inp_kind == 1  ) nomvar= 'P0'
            !if (Inp_kind == 5  ) nomvar= 'P0'
            !if (Inp_src_hauteur_L ) nomvar= 'P0'
            if ( nomvar == 'P0' ) mult= 100.d0
         case ('TEMPERATURE')
            nomvar= 'TT'
            !if (Inp_kind == 2  ) nomvar= 'TT'
            !if (Inp_kind == 1  ) nomvar= 'TT'
            !if (Inp_kind == 5  ) nomvar= 'TT'
            !if (Inp_src_hauteur_L ) nomvar= 'TT'
            if ( nomvar == 'TT' ) add= tcdk_8
         case ('GEOPOTENTIAL')
            nomvar= 'GZ' ; mult= 10.d0
         case ('PX')
            mult= 100.d0
         case ('URT1')
            mult= knams_8
         case ('VRT1')
            mult= knams_8
      end select

      datev= Inp_cmcdate
      if ( F_var_S(1:min(3,len_trim(F_var_S))) == 'TR/' ) then
         nomvar= F_var_S(4:)
         if (Tr3d_anydate_L) datev= -1
      end if

      F_nomvar_S= trim(nomvar)
      if (typ /= " ") F_nomvar_S= trim(F_nomvar_S)//':'//typ

      if ( nomvar == '@NUL' ) return

      query = Inp_file%new_query(datev=datev,nomvar=nomvar,&
                                 ip1=ip1,typvar=typ)
      lislon = query%find_all(recs)

      if (lislon == 0) goto 999

      src_gid= ezqkdef (recs(1)%ni,recs(1)%nj,recs(1)%grtyp,&
                        recs(1)%ig1,recs(1)%ig2,recs(1)%ig3,&
                        recs(1)%ig4,Inp_file%get_unit())

      if ((trim(nomvar) == 'URT1').or.(trim(nomvar) == 'VRT1').or.&
          (trim(nomvar) == 'UT1' ).or.(trim(nomvar) == 'VT1' )) then
         err= samegrid_rot ( src_gid, &
         Rot_ig1, Rot_ig2, Rot_ig3, Rot_ig4)
         if (err < 0) then
            lislon= 0
            goto 999
         end if
      end if

      call record_sort_ip1 (recs,liste_sorted,lislon)

      if (lislon > 1) then
         F_ip1(1:lislon) = liste_sorted(1:lislon)
      else
         F_ip1(1) = recs(1)%ip1
      end if

      F_nka= lislon
      if (Grd_yinyang_L) then
         call splitW ( YYG_myproc, YYG_numproc,lislon,1,&
                       local_nk, kstart, i)
      else
         call splitW ( myproc_IOS, numproc_IOS,lislon,1,&
                       local_nk, kstart, i)
      endif
      if (kstart<0) goto 999

      !allocate (wk1(n1*n2))
      allocate (wk1(recs(1)%ni*recs(1)%nj))
      interp_S= 'CUBIC'
      if (present(F_hint_S)) interp_S= F_hint_S
      
      do idst= 1, F_nd
          
          if (local_nk > 0) then
             if (F_hgrid_S(idst) == 'Q') then
                posx => geomh_longs
                posy => geomh_latgs
             end if
             if (F_hgrid_S(idst) == 'U') then
                posx => geomh_longu
                posy => geomh_latgs
             end if
             if (F_hgrid_S(idst) == 'V') then
                posx => geomh_longs
                posy => geomh_latgv
             end if
             if (F_hgrid_S(idst) == 'F') then
                posx => geomh_longu
                posy => geomh_latgv
             end if
             ni_dest= G_ni+2*G_halox
             nj_dest= G_nj+2*G_haloy
             dstf_gid = ezgdef_fmem (ni_dest, nj_dest, 'Z', 'E', &
                           F_ig1, F_ig2, F_ig3, F_ig4, posx, posy)

            if ( recs(1)%grtyp == 'U' ) then
                nicore = G_ni-Glb_pil_w-Glb_pil_e
                njcore = G_nj-Glb_pil_s-Glb_pil_n
                if (recs(1)%ni >= nicore .and. recs(1)%nj/2 >= njcore) then
                   subid= samegrid_gid ( src_gid, F_ig1,F_ig2,F_ig3,F_ig4,&
                      posx(1+Glb_pil_w), posy(1+Glb_pil_s), nicore,njcore )
                else
                   subid=-1
                end if
                if (subid >= 0) then
                   interp_S = 'NEAREST'
                   err = ezsetopt ('USE_1SUBGRID', 'YES')
                   err = ezsetival('SUBGRIDID', subid)
                else
                   err = ezsetopt ('USE_1SUBGRID', 'NO')
                end if
             end if

             err = ezdefset ( dstf_gid , src_gid )
             err = ezsetopt ('INTERP_DEGREE', interp_S)
             if (lun_out>0) write(lun_out,1001) &
                'Interpolating: ',trim(F_var_S),trim(nomvar),', nka= ',&
                lislon,',valid: ',Inp_datev,' on ',F_hgrid_S(idst),&
                'grid, levels:',kstart,kstart+local_nk-1
          end if

          do i=1,local_nk
             success=recs(kstart+i-1)%read(data=c_loc(wk1))
             k= (idst-1)*lislon+kstart+i-1
             err = ezsint (F_dest(1,1,k),wk1)
             F_dest(:,:,k)= F_dest(:,:,k)*mult + add
          end do
          if (lun_out>0) then
             if (err==2) &
             write(lun_out,1001) &
             'EXTRApolating: ',trim(F_var_S),trim(nomvar),', nka= ',&
             lislon,',valid: ',Inp_datev,' on ',F_hgrid_S(idst),' grid'
          endif
          err= ezsetopt ( 'USE_1SUBGRID', 'NO' )
      end do
      deallocate (wk1)

 999  if (lislon > 0) then
         INs_read_mt= 0
      else
         if ((.not.quiet_L).and.(lun_out>0)) write(lun_out,'(7a)') &
              ' FIELD: ',trim(F_var_S),':',trim(nomvar),' valid: ',&
              Inp_datev, 'NOT FOUND'
      end if
      call MPI_barrier (MY_WORLD_COMM,err)

 1001 format (2a,':',2a,i3,5a,2i3)
!
!---------------------------------------------------------------------
!
      return
      end function INs_read_mt
    
!**s/r INs_read_uv - Read UU and VV and perform horizontal
!                    interpolation to U and V points respectively

      integer function INs_read_uv (F_dest, F_ip1, F_nka, &
                               F_ig1, F_ig2, F_ig3, F_ig4 )
      implicit none

      integer, intent(IN ) :: F_ig1, F_ig2, F_ig3, F_ig4
      integer, intent(OUT) :: F_nka
      integer, intent(OUT) :: F_ip1(*)
      real, intent(OUT) :: F_dest(G_ni+2*G_halox,G_nj+2*G_haloy,*)

      character(len=1) typ,grd
      character(len=4) var,dumc
      character(len=12) lab
      integer, parameter :: nlis = 1024
      integer :: i, k, idst, err, nz, n1,n2,n3, nrec
      integer :: liste_u(nlis),liste_v(nlis),liste_sorted(nlis)
      integer :: datev,local_nk,src_gid,nku,nkv
      integer :: dstu_gid,dstv_gid,kstart,erru,errv
      integer :: dte, det, ipas, p1, p2, p3, g1, g2, g3, g4, bit, &
                 dty, swa, lng, dlf, ubc, ex1, ex2, ex3
      real, dimension(:), pointer     :: posxu,posyu,posxv,posyv
      real, dimension(:), allocatable, target :: uv,u,v

      type(fst_record) :: recs_u(nlis),recs_v(nlis)
      type(fst_query)  :: query_u,query_v
      logical          :: success
!
!---------------------------------------------------------------------
!
      INs_read_uv= -1
      F_nka= -1 ; local_nk= 0

      datev= Inp_cmcdate

      query_u = Inp_file%new_query(datev=Inp_cmcdate,nomvar='UU  ')
      nku = query_u%find_all(recs_u)

      query_v = Inp_file%new_query(datev=Inp_cmcdate,nomvar='VV  ')
      nkv = query_v%find_all(recs_v)
      if ((nku/=nkv).or.(nku<3)) goto 999

      src_gid = ezqkdef (recs_u(1)%ni,recs_u(1)%nj,recs_u(1)%grtyp,&
                         recs_u(1)%ig1,recs_u(1)%ig2,recs_u(1)%ig3,&
                         recs_u(1)%ig4,Inp_file%get_unit())

      call record_sort_ip1 (recs_u,liste_sorted,nku)
      call record_sort_ip1 (recs_v,liste_sorted,nkv)

      F_ip1(1:nku) = liste_sorted(1:nku)
      F_ip1(nku+1:2*nku) = F_ip1(1:nku)
      
      F_nka= 2*nku
      if (Grd_yinyang_L) then
         call splitW ( YYG_myproc, YYG_numproc,nku,1,&
                       local_nk, kstart, i)
      else
         call splitW ( myproc_IOS, numproc_IOS,nku,1,&
                       local_nk, kstart, i)
      endif
      if (kstart<0) goto 999

      allocate (u(recs_u(1)%ni*recs_u(1)%nj), &
                v(recs_u(1)%ni*recs_u(1)%nj))
      allocate (uv(INs_nid*INs_njd))

      err = ezsetopt ('INTERP_DEGREE', 'CUBIC')
      posxu => geomh_longu
      posyu => geomh_latgs
      posxv => geomh_longs
      posyv => geomh_latgv

      if (local_nk > 0) then
         if (lun_out>0) write(Lun_out,1001) 'Interpolating: UU, nka= ',&
                    nku,', valid: ',Inp_datev,' on U grid'
         dstu_gid = ezgdef_fmem ( INs_nid, INs_njd, 'Z', 'E', &
                                  F_ig1, F_ig2, F_ig3, F_ig4, &
                                                 posxu, posyu )
         if (lun_out>0) write(Lun_out,1001) 'Interpolating: VV, nka= ',&
                    nkv,', valid: ',Inp_datev,' on V grid'
         dstv_gid = ezgdef_fmem ( INs_nid, INs_njd, 'Z', 'E', &
                                  F_ig1, F_ig2, F_ig3, F_ig4, &
                                                 posxv, posyv )
         do i=1,local_nk
            k= kstart+i-1
            success = recs_u(k)%read(data=c_loc(u))
            success = recs_v(k)%read(data=c_loc(v))
            err = ezdefset ( dstu_gid , src_gid )
            erru= ezuvint  ( F_dest(1,1,k),uv, u,v )
            F_dest(:,:,k)= F_dest(:,:,k) * knams_8
            err = ezdefset ( dstv_gid , src_gid )
            errv= ezuvint  ( uv,F_dest(1,1,nku+k), u,v )
            F_dest(:,:,nku+k)= F_dest(:,:,nku+k) * knams_8
         end do
         if ((erru==2).and.(Lun_out>0)) &
            write(Lun_out,1002) 'EXTRApolating: UU, nka= ',&
                      nku,', valid: ',Inp_datev,' on U grid'
         if ((errv==2).and.(Lun_out>0)) &
            write(Lun_out,1002) 'EXTRApolating: VV, nka= ',&
                      nkv,', valid: ',Inp_datev,' on V grid'
      endif
      deallocate (u,v,uv)

 999  if (nku > 0) then
         INs_read_uv= 0
      else
         if (Lun_out>0) write(Lun_out,'(3a)') &
         'Variable: UU,VV valid: ',Inp_datev, 'NOT FOUND'
      end if
      call MPI_barrier (MY_WORLD_COMM,err)

 1001 format (a,i3,3a)
 1002 format (a,i3,3a)
!
!---------------------------------------------------------------------
!
      return
      end function INs_read_uv

!**s/r INs_hwnd - Read and interpolate horizontal winds UU,VV
      
      integer function INs_hwnd (F_vname, F_nk, F_deb, F_data, &
                                 F_ig1, F_ig2, F_ig3, F_ig4, F_dim )
      implicit none

      character(len=*), intent (OUT) :: F_vname
      integer, intent (IN ) :: F_ig1, F_ig2, F_ig3, F_ig4, F_dim
      integer, intent (OUT) :: F_nk, F_deb
      real   , intent (OUT) :: F_data(F_dim)
      
      character(len=6) vname
      integer nkau,nkav,err,deb1,deb2,initial
!
!---------------------------------------------------------------------
!
      INs_hwnd= -1 ; initial= INs_nplans
      F_vname= '' ; F_nk=0
      deb1= INs_nplans+1
      deb2= INs_nplans*INs_hord+1
      F_deb= deb1
      err= INs_read_mt ('URT1','U',vname,F_data(deb2:),1,DIP1(deb1),nkau,F_ig1, F_ig2, F_ig3, F_ig4)
      if (nkau>0) then
         INs_nplans= INs_nplans+nkau
         deb1= INs_nplans+1
         deb2= INs_nplans*INs_hord+1
         err= INs_read_mt ('VRT1','V',vname,F_data(deb2:),1,DIP1(deb1),nkav,F_ig1, F_ig2, F_ig3, F_ig4)
         if (nkav/=nkau) then
            print*, 'URT1 / VRT1 are unusable'
            INs_nplans= initial
            F_deb=-1
         else
            INs_nplans= INs_nplans+nkav
            INs_hwnd= 0
            F_vname= 'UVRT1' ; F_nk=2*nkau
         endif
      endif

      if (INs_hwnd/=0) then
         deb1= INs_nplans+1
         deb2= INs_nplans*INs_hord+1
         F_deb= deb1
         err= INs_read_uv (F_data(deb2:),DIP1(deb1),nkau,&
                           F_ig1, F_ig2, F_ig3, F_ig4)
         if (nkau>0) then
            INs_nplans= INs_nplans+nkau
            INs_hwnd=0
            F_vname= 'UV' ; F_nk=nkau
         endif
      end if
!     
!---------------------------------------------------------------------
!
      return
      end function INs_hwnd

end module INs_base
