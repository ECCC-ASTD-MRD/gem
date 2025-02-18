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

module inp_base
   use iso_c_binding
   use vertical_interpolation
   use vGrid_Descriptors
   use gem_options
   use inp_options
   use geomh
   use glb_ld
   use glb_pil
   use inp_mod
   use hgc
   use lun
   use ver
   use tdpack
   use, intrinsic :: iso_fortran_env
   implicit none
#include <arch_specific.hf>
#include <rmnlib_basics.hf>

   public

contains

!**s/r inp_get - Read variable F_var_S valid at F_datev, perform hor.
!                interpolation to F_hgrid_S hor. grid and perform
!                vertical interpolation to F_vgrid_S vertical grid

      integer function inp_get ( F_var_S, F_hgrid_S, F_ver_ip1         ,&
                         F_sfc_src, F_sfcLS_src, F_sfc_dst, F_sfcLS_dst,&
                         F_gz, F_GZ_ip1, F_dest , Minx,Maxx,Miny,Maxy  ,&
                         F_nk, F_inttype_S, F_quiet_L )

      implicit none

      character(len=*)          , intent(in) :: F_var_S,F_hgrid_S
      character(len=*), optional, intent(in) :: F_inttype_S
      logical         , optional, intent(in) :: F_quiet_L
      integer                   , intent(in) :: Minx,Maxx,Miny,Maxy, F_nk
      integer, dimension(:)  , pointer, intent(in) :: F_ver_ip1,F_gz_ip1
      real, dimension (:,:), pointer, intent(in) :: &
                              F_sfc_src,F_sfcLS_src,F_sfc_dst,F_sfcLS_dst
      real, dimension(Minx:Maxx,Miny:Maxy,F_nk), intent(out):: F_dest
      real, dimension(:,:,:), pointer, intent(in ) :: F_gz

!     local variables
      character(len=12) :: inttype
      logical quiet_L
      integer nka
      integer, dimension (:    ), pointer :: ip1_list
      real   , dimension (:,:,:), pointer :: wrkr,srclev,dstlev
!
!---------------------------------------------------------------------
!
      inp_get= -1
      if ( any (Inp_blacklist_S(1:MAX_blacklist) == trim(F_var_S)) ) &
           return

      nullify (ip1_list, wrkr)
      quiet_L=.false.
      if (present(F_quiet_L)) quiet_L= F_quiet_L

      inp_get= inp_read_mt ( F_var_S, F_hgrid_S, wrkr, 1, &
                             ip1_list, nka, F_quiet_L=quiet_L )

      if (inp_get < 0) then
         if (associated(ip1_list)) deallocate (ip1_list)
         if (associated(wrkr    )) deallocate (wrkr    )
         return
      end if

      allocate ( srclev(l_minx:l_maxx,l_miny:l_maxy,nka) ,&
                 dstlev(l_minx:l_maxx,l_miny:l_maxy,F_nk) )

      call inp_src_levels (srclev, nka, ip1_list, Inp_vgd_src, &
         F_sfc_src, F_sfcLS_src, F_gz, F_GZ_ip1,Minx,Maxx,Miny,Maxy)

      call inp_dst_levels ( dstlev, Ver_vgdobj, F_ver_ip1,&
                            F_sfc_dst, F_sfcLS_dst, F_nk )

      inttype= 'cubic'
      if (present(F_inttype_S)) inttype= F_inttype_S
      call vertint2 ( F_dest,dstlev,F_nk, wrkr,srclev,nka           ,&
                      l_minx,l_maxx,l_miny,l_maxy                   ,&
                      1-G_halox,l_ni+G_halox, 1-G_haloy,l_nj+G_haloy,&
                      varname=F_var_S, inttype= inttype,&
                      levtype=Inp_levtype_S )

      deallocate (ip1_list,wrkr,srclev,dstlev)
!
!---------------------------------------------------------------------
!
      return
      end function inp_get
!
!**s/r inp_read_mt - Parallel read of variable F_var_S and horizontal
!                    interpolation to F_nd Arakawa grid destinations

      integer function inp_read_mt ( F_var_S, F_hgrid_S, F_dest, &
                         F_nd, F_ip1, F_nka, F_hint_S, F_quiet_L )

      use rmn_fst24

      implicit none

      character(len=*)          ,intent(in)  :: F_var_S
      character(len=*), dimension(*),intent(in)  :: F_hgrid_S
      character(len=*), optional,intent(in)  :: F_hint_S
      logical         , optional,intent(in)  :: F_quiet_L
      integer                   ,intent(in ) :: F_nd
      integer                   ,intent(out) :: F_nka
      integer, dimension(:    ), pointer,intent(inout) :: F_ip1
      real   , dimension(:,:,:), pointer,intent(inout) :: F_dest

      integer, external :: samegrid_gid, samegrid_rot, inp_is_real_wind
      character(len=1) typ
      character(len=4) nomvar,var,dumc
      character(len=12) lab,interp_S
      logical :: quiet_L
      integer, parameter :: nlis = 1024

      type(fst_query)  :: query
      type(fst_record) :: recs(nlis) 
      logical          :: success
      
       integer i, k, idst, err, nz, &
              liste_sorted(nlis),lislon,maxdim_wk2
      integer ni_dest,nj_dest
      integer subid,nicore,njcore,datev
      integer mpx,local_nk,irest,kstart, src_gid, vcode, ip1, p1

      integer, dimension(:  ), allocatable :: zlist
      real :: surface_level
      real   , dimension(:  ), allocatable, target :: wk1
      real   , dimension(:,:), allocatable :: wk4
      real   , dimension(:  ), pointer     :: posx,posy
      real(kind=REAL64) add, mult
      common /bcast_i / lislon,nz
!
!---------------------------------------------------------------------
!
      inp_read_mt= -1
      F_nka= -1 ; local_nk= 0
      add= 0.d0 ; mult= 1.d0
      if (associated(F_ip1 )) deallocate (F_ip1 )
      if (associated(F_dest)) deallocate (F_dest)
      nullify (F_ip1,F_dest)
      quiet_L=.false.
      if (present(F_quiet_L)) quiet_L= F_quiet_L

      nomvar = F_var_S ; ip1= -1
      select case (F_var_S)
         case ('OROGRAPHY')
            if (Inp_kind == 2  ) then
               nomvar= '@NUL'
               if (Inp_src_PX_L) then
                  nomvar= 'GZ' ; surface_level= 1. ; p1=5
                  call convip ( ip1, surface_level,p1,1,dumc,.false. )
               endif
            endif
            if (Inp_kind == 1 ) then
               nomvar= 'GZ' ; ip1= 12000
               if (Inp_src_PX_L) then
                  nomvar= 'GZ' ; surface_level= 1. ; p1=5
                  call convip ( ip1, surface_level,p1,1,dumc,.false. )
               endif
            endif
            if (Inp_kind == 5 ) then
               nomvar= 'GZ' ; surface_level= 1.
               call convip ( ip1, surface_level,Inp_kind,1,dumc,.false. )
            endif
            if ( Inp_src_hauteur_L ) then
               nomvar= 'GZ'
               if (Inp_kind==21) surface_level= 0.
               if (Inp_kind==5 ) surface_level= 1.
               call convip ( ip1, surface_level,Inp_kind,1,dumc,.false. )
            endif
            if ( nomvar == 'GZ' ) mult= 10.d0 * grav_8
         case ('SFCPRES')
            if (Inp_kind == 2  ) nomvar= '@NUL'
            if (Inp_kind == 1  ) nomvar= 'P0'
            if (Inp_kind == 5  ) nomvar= 'P0'
            if (Inp_src_hauteur_L ) nomvar= 'P0'
            if ( nomvar == 'P0' ) mult= 100.d0
         case ('TEMPERATURE')
            if (Inp_kind == 2  ) nomvar= 'TT'
            if (Inp_kind == 1  ) nomvar= 'TT'
            if (Inp_kind == 5  ) nomvar= 'TT'
            if (Inp_src_hauteur_L ) nomvar= 'TT'
            if ( nomvar == 'TT' ) add= tcdk_8
         case ('GEOPOTENTIAL')
            nomvar= 'GZ' ; mult= 10.d0
         case ('PX')
            mult= 100.d0
         case ('UU')
            mult= knams_8
         case ('VV')
            mult= knams_8
      end select

      datev= Inp_cmcdate
      if ( F_var_S(1:min(3,len_trim(F_var_S))) == 'TR/' ) then
         nomvar= F_var_S(4:)
         if (Tr3d_anydate_L) datev= -1
      end if

      if ( nomvar == '@NUL' ) return

      nz= -1
      maxdim_wk2 = 1
      if (Inp_iome >= 0) then
         vcode= -1
         query = Inp_file%new_query(datev=datev,nomvar=nomvar,ip1=ip1)
         lislon = query%find_all(recs)

         if (lislon == 0) goto 999

         src_gid= ezqkdef (recs(1)%ni,recs(1)%nj,recs(1)%grtyp,recs(1)%ig1,recs(1)%ig2,recs(1)%ig3,recs(1)%ig4,Inp_file%get_unit())

         if ((trim(nomvar) == 'URT1').or.(trim(nomvar) == 'VRT1').or.&
             (trim(nomvar) == 'UT1' ).or.(trim(nomvar) == 'VT1' )) then
             err= samegrid_rot (src_gid, Hgc_ig1ro, Hgc_ig2ro, Hgc_ig3ro, Hgc_ig4ro)
             if (err < 0) then
                lislon= 0
                goto 999
             end if
         end if

         call record_sort_ip1 (recs,liste_sorted,lislon)

         allocate (F_ip1(max(1,lislon)))
         if (lislon > 1) then
            F_ip1(1:lislon) = liste_sorted(1:lislon)
         else
            F_ip1(1) = recs(1)%ip1
         end if

         nz= (lislon + Inp_npes - 1) / Inp_npes
!        lislon should be >0 here so, nz  >=1

         maxdim_wk2=nz*F_nd
         mpx      = mod( Inp_iome, Inp_npes )
         local_nk = lislon / Inp_npes
         irest  = lislon  - local_nk * Inp_npes
         kstart = mpx * local_nk + 1
         if ( mpx < irest ) then
            local_nk   = local_nk + 1
            kstart = kstart + mpx
         else
            kstart = kstart + irest
         end if

         ni_dest= G_ni+2*G_halox
         nj_dest= G_nj+2*G_haloy
         allocate (wk4(ni_dest*nj_dest,maxdim_wk2))
         allocate (wk1(recs(1)%ni*recs(1)%nj))

         interp_S= 'CUBIC'
         if (present(F_hint_S)) interp_S= F_hint_S

         do idst= 1, F_nd !IDST loop
         if (local_nk > 0) then

            if (F_hgrid_S(idst) == 'Q') then
               posx => geomh_lonQ
               posy => geomh_latQ
            end if
            if (F_hgrid_S(idst) == 'U') then
               posx => geomh_lonF
               posy => geomh_latQ
            end if
            if (F_hgrid_S(idst) == 'V') then
               posx => geomh_lonQ
               posy => geomh_latF
            end if
            if (F_hgrid_S(idst) == 'F') then
               posx => geomh_lonF
               posy => geomh_latF
            end if

            dstf_gid = ezgdef_fmem (ni_dest, nj_dest, 'Z', 'E', &
                    Hgc_ig1ro, Hgc_ig2ro, Hgc_ig3ro, Hgc_ig4ro, &
                    posx, posy)

            if ( recs(1)%grtyp == 'U' ) then
              nicore = G_ni-Glb_pil_w-Glb_pil_e
              njcore = G_nj-Glb_pil_s-Glb_pil_n
              if (recs(1)%ni >= nicore .and. recs(1)%nj/2 >= njcore) then
                 subid= samegrid_gid ( &
                    src_gid, Hgc_ig1ro,Hgc_ig2ro,Hgc_ig3ro,Hgc_ig4ro,&
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
            write(output_unit,1001) 'Interpolating: ',trim(F_var_S),trim(nomvar),', nka= ',&
               lislon,',valid: ',Inp_datev,' on ',F_hgrid_S(idst),' grid'
         end if

         err = -1
         do i=1,local_nk
            success=recs(kstart+i-1)%read(data=c_loc(wk1))
            err = ezsint(wk4(1,(idst-1)*nz+i), wk1)
         end do
         if (err == 2) then
            write(output_unit,1001) 'EXTRApolating: ',trim(F_var_S),trim(nomvar),', nka= ',&
               lislon,',valid: ',Inp_datev,' on ',F_hgrid_S(idst),' grid'
         end if
         err = ezsetopt ( 'USE_1SUBGRID', 'NO' )

         END DO !IDST loop

         deallocate (wk1)
      else !Inp_iome >= 0
         allocate (wk4(1,1))
         maxdim_wk2 = 1
      end if !Inp_iome >= 0

 999  call rpn_comm_bcast ( lislon, 2, "MPI_INTEGER", Inp_iobcast, &
                            "grid", err ) !NOTE: bcast lislon AND nz
      F_nka= lislon

      if (F_nka > 0) then

         inp_read_mt= 0
         if (F_nka >= 1) then
            if (Inp_iome < 0) allocate ( F_ip1(F_nka) )
            call rpn_comm_bcast ( F_ip1, F_nka, "MPI_INTEGER", &
                                  Inp_iobcast, "grid", err )
         end if

         allocate (zlist(nz)) ; zlist= -1
         do i=1, local_nk
            zlist(i)= i + kstart - 1
         end do

         allocate ( F_dest(l_minx:l_maxx,l_miny:l_maxy,lislon*F_nd) );F_dest=0.

         do idst=1, F_nd
            k=min((idst-1)*nz+1,maxdim_wk2)
            call glbdist_os (wk4(1,k),F_dest(l_minx,l_miny,(idst-1)*lislon+1),&
                             l_minx,l_maxx,l_miny,l_maxy,F_nka,&
                             G_ni+G_halox,G_nj+G_haloy,zlist,nz,mult,add)
         end do
         deallocate (wk4,zlist)

      else

         inp_read_mt= -1
         if ((Inp_iome >= 0).and.(.not.quiet_L)) write(output_unit,'(7a)') &
            ' FIELD: ',trim(F_var_S),':',trim(nomvar),' valid: ',&
            Inp_datev, 'NOT FOUND'

      end if

 1001 format (2a,':',2a,i3,5a)
!
!---------------------------------------------------------------------
!
      return
      end function inp_read_mt

!**s/r inp_oro - Read orography from input dataset valid at F_datev

      subroutine inp_oro ( F_topo, F_topo_LS, F_meqr, F_datev,&
                           Minx, Maxx, Miny, Maxy )
      use gmm_geof
      use glb_ld
      use rmn_gmm
      use dyn_fisl_options
      use HORgrid_options
      use lam_options
      use mem_nest
      use step_options
      implicit none

      character(len=*), intent(in) :: F_datev
      integer         , intent(in) :: Minx, Maxx, Miny, Maxy
      real, dimension(Minx:Maxx,Miny:Maxy), intent(out) :: F_topo,F_topo_LS
      real, dimension (:,:,:), pointer    , intent(inout) :: F_meqr

      integer i,j,err,err_ls,nka
      integer, dimension (:), pointer :: ip1_list
      real, dimension (:,:,:), pointer :: wrk
      real, dimension (:,:), pointer :: ls
      real step_current
      real(kind=REAL64) diffd
!
!---------------------------------------------------------------------
!
      if (associated(F_meqr)) deallocate (F_meqr)
      nullify (F_meqr, wrk, ip1_list)
      err_ls= -1

      err = inp_read_mt ( 'OROGRAPHY', 'Q', wrk, 1, ip1_list, nka )

      if ( associated(ip1_list) ) then
         deallocate (ip1_list) ; nullify (ip1_list)
      end if
      if ( associated(wrk) ) then
         allocate (F_meqr(l_minx:l_maxx,l_miny:l_maxy,2))
         F_meqr(:,:,1) = wrk(:,:,1)
         deallocate (wrk) ; nullify (wrk)
         err_ls = inp_read_mt ( 'MELS', 'Q', wrk, 1, ip1_list, nka, F_quiet_L=.true.)
         if ( associated(wrk) ) then
            F_meqr(:,:,2) = wrk(:,:,1)
            deallocate (wrk,ip1_list) ; nullify (wrk,ip1_list)
         endif
      endif
      if ( trim(F_datev) == trim(Step_runstrt_S) ) then
         if ( associated(F_meqr) ) then
            topo_low(:,:,1) = F_meqr(:,:,1)
            if (err_ls == 0) then
               topo_low(:,:,2) = F_meqr(:,:,2)
            else
               call mc2_topols (topo_low(:,:,2),F_meqr,&
                       l_minx,l_maxx,l_miny,l_maxy,Schm_orols_np)
            endif
        else
            Vtopo_L= .false.
         end if
      end if

      call difdatsd (diffd,Step_runstrt_S,F_datev)
      step_current = diffd*86400.d0 / Step_dt + Step_initial
      call var_topo ( F_topo, F_topo_LS, step_current, &
                      l_minx,l_maxx,l_miny,l_maxy )

      if ( associated(F_meqr) .and. .not. Grd_yinyang_L) then
      if ( Lam_blendoro_L ) then
         allocate (ls(l_minx:l_maxx,l_miny:l_maxy))
         if (err_ls == 0) then
            ls(:,:) = F_meqr(:,:,2)
         else
            call mc2_topols (ls,F_meqr,&
                   l_minx,l_maxx,l_miny,l_maxy,Schm_orols_np)
         endif
         do j= 1-G_haloy, l_nj+G_haloy
            do i= 1-G_halox, l_ni+G_halox
               F_topo(i,j)= F_topo(i,j)*(1.-nest_weightm(i,j,G_nk+1)) +&
                            F_meqr(i,j,1)*nest_weightm(i,j,G_nk+1)
               F_topo_LS(i,j)= F_topo_LS(i,j)*(1.-nest_weightm(i,j,G_nk+1)) +&
                               ls(i,j)*nest_weightm(i,j,G_nk+1)
            enddo
         enddo
         deallocate (ls) ; nullify (ls)
      end if
      end if
!
!---------------------------------------------------------------------
!
      return
      end subroutine inp_oro

!**s/r inp_tv - Read TT and HU from input dataset and compute TV

      subroutine inp_tv ( F_tv, F_tt, F_hu, F_TT_ip1, &
                          F_HU_ip1, F_nka_tt, F_nka_hu )
      implicit none

      integer, intent(  out) :: F_nka_tt, F_nka_hu
      integer, dimension (:), pointer, intent(inout) :: F_TT_ip1,F_HU_ip1
      real   , dimension (:,:,:), pointer, intent(inout) :: &
                                             F_tv, F_tt, F_hu
      integer err(2),kt,kh
!
!---------------------------------------------------------------------
!
      if (associated(F_tv)) deallocate(F_tv)
      if (associated(F_tt)) deallocate(F_tt)
      if (associated(F_hu)) deallocate(F_hu)
      if (associated(F_TT_ip1)) deallocate(F_TT_ip1)
      if (associated(F_HU_ip1)) deallocate(F_HU_ip1)
      nullify (F_tv, F_tt, F_hu, F_TT_ip1, F_HU_ip1)

      err(1) = inp_read_mt ( 'TR/HU'      ,'Q',F_hu,1,F_HU_ip1,F_nka_hu )
      err(2) = inp_read_mt ( 'TEMPERATURE','Q',F_tt,1,F_TT_ip1,F_nka_tt )

      call gem_error ( minval(err(:)),'inp_tv','MISSING DATA')

      allocate (F_tv(l_minx:l_maxx,l_miny:l_maxy,F_nka_tt))
      F_tv=F_tt

      do kt=1, F_nka_tt
         inner:      do kh=1, F_nka_hu
         if (F_TT_ip1(kt) == F_HU_ip1(kh)) then
            call mfottv2 (F_tv(l_minx,l_miny,kt),F_tv(l_minx,l_miny,kt),&
                  F_hu(l_minx,l_miny,kh),l_minx,l_maxx,l_miny,l_maxy,1 ,&
                  1-G_halox,l_ni+G_halox,1-G_haloy,l_nj+G_haloy, .true. )
            exit inner
         end if
      end do inner
      end do
!     
!---------------------------------------------------------------------
!
      return
      end subroutine inp_tv

!**s/r inp_src_surface - Read surface information from input dataset

      subroutine inp_src_surface ( F_sq,F_su,F_sv,F_LSsq,F_LSsu,F_LSsv,F_topo,F_nka )
      use dynkernel_options
      implicit none

      integer, intent(in) :: F_nka
      real   , dimension (:,:), pointer, intent(inout) :: F_sq,F_su,F_sv,F_LSsq,F_LSsu,F_LSsv
      real   , dimension (*  ),          intent(in   ) :: F_topo

      integer err,err2,i,j,k,nk,kind
      integer, dimension (:    ), pointer     :: ip1_list
      real   , dimension (:,:,:), pointer     :: wrk
      real   , dimension (:    ), allocatable :: rna
      character(len=1) :: grid_S(3)
!
!---------------------------------------------------------------------
!
      if (associated(F_sq)) deallocate (F_sq)
      if (associated(F_su)) deallocate (F_su)
      if (associated(F_sv)) deallocate (F_sv)
      if (associated(F_LSsq)) deallocate (F_LSsq)
      if (associated(F_LSsu)) deallocate (F_LSsu)
      if (associated(F_LSsv)) deallocate (F_LSsv)
      nullify (F_sq,F_su,F_sv,F_LSsq,F_LSsu,F_LSsv,ip1_list,wrk)

      if (Inp_dst_hauteur_L.and..not.Schm_autobar_L) then

         err = -1
         nullify (wrk,ip1_list)

         err = inp_read_mt ( 'SFCPRES', 'Q', wrk, 1,&
                                ip1_list, nk )

         if (associated(wrk)) then
            allocate ( F_sq(l_minx:l_maxx,l_miny:l_maxy) )
            F_sq(:,:) = wrk(:,:,1)
            deallocate (wrk,ip1_list) ; nullify (wrk,ip1_list)
            err = 0
         endif

         call gem_error ( err,'inp_src_gz', 'MISSING SURFACE P0 DATA')

         return
      endif

      grid_S=['Q','U','V']

      err = inp_read_mt ( 'SFCPRES', grid_S, wrk, 3,&
                           ip1_list, nk)
      if (associated(wrk)) then
         allocate ( F_sq(l_minx:l_maxx,l_miny:l_maxy),&
                    F_su(l_minx:l_maxx,l_miny:l_maxy),&
                    F_sv(l_minx:l_maxx,l_miny:l_maxy) )
         F_sq(:,:) = wrk(:,:,1)
         F_su(:,:) = wrk(:,:,2)
         F_sv(:,:) = wrk(:,:,3)
         deallocate (wrk,ip1_list) ; nullify (wrk,ip1_list)
         err2 = inp_read_mt ( 'P0LS', grid_S, wrk, 3,&
                              ip1_list, nk,F_quiet_L=.true.)
         if (associated(wrk)) then
            allocate ( F_LSsq(l_minx:l_maxx,l_miny:l_maxy),&
                       F_LSsu(l_minx:l_maxx,l_miny:l_maxy),&
                       F_LSsv(l_minx:l_maxx,l_miny:l_maxy) )
            F_LSsq(:,:) = wrk(:,:,1)*100.
            F_LSsu(:,:) = wrk(:,:,2)*100.
            F_LSsv(:,:) = wrk(:,:,3)*100.
            deallocate (wrk,ip1_list) ; nullify (wrk,ip1_list)
         endif
      else
         if (Inp_kind == 2.or.Schm_autobar_L) then
            if (Inp_src_PX_L) then
               allocate ( F_sq(l_minx:l_maxx,l_miny:l_maxy),&
                          F_su(l_minx:l_maxx,l_miny:l_maxy),&
                          F_sv(l_minx:l_maxx,l_miny:l_maxy) )
               F_sq(:,:) = PX3d%valq(:,:,PX3d%nk)
               F_su(:,:) = PX3d%valu(:,:,PX3d%nk)
               F_sv(:,:) = PX3d%valv(:,:,PX3d%nk)
               err=0
            else

            err = inp_read_mt ( 'GEOPOTENTIAL', 'Q', wrk, 1,&
                                 ip1_list, nk )
            if (nk == F_nka) then
               allocate (rna(nk))
               do k=1,nk
                  call convip(ip1_list(k),rna(k),kind,-1,' ',.false.)
               end do
               allocate ( F_sq(l_minx:l_maxx,l_miny:l_maxy), &
                          F_su(l_minx:l_maxx,l_miny:l_maxy), &
                          F_sv(l_minx:l_maxx,l_miny:l_maxy) )
               call gz2p0 ( F_sq, wrk, F_topo, rna     ,&
                            l_minx,l_maxx,l_miny,l_maxy,&
                      nk,1-G_halox,l_ni+G_halox,1-G_haloy,l_nj+G_haloy )
               do j=1-G_haloy,l_nj+G_haloy
                  do i=1-G_halox,l_ni+G_halox-1
                     F_su(i,j)= (F_sq(i,j) + F_sq(i+1,j)) * 0.5d0
                  end do
               end do
               if (l_east) F_su(l_ni+G_halox,:)=F_sq(l_ni+G_halox,:)
               do j=1-G_haloy,l_nj+G_haloy-1
                  do i=1-G_halox,l_ni+G_halox
                     F_sv(i,j)= (F_sq(i,j) + F_sq(i,j+1)) * 0.5d0
                  end do
               end do
               if (l_north) F_sv (:,l_nj+G_haloy)= F_sq(:,l_nj+G_haloy)
               deallocate (rna,wrk,ip1_list) ; nullify (wrk,ip1_list)
               call gem_xch_halo (F_su, l_minx, l_maxx, l_miny, l_maxy, 1)
               call gem_xch_halo (F_sv, l_minx, l_maxx, l_miny, l_maxy, 1)
            else
               err = -1
            end if
            end if
         end if
      endif

      call gem_error ( err,'inp_src_surface','MISSING SURFACE DATA' )
!
!---------------------------------------------------------------------
!
      return
      end subroutine inp_src_surface

!**s/r inp_dst_surface - compute destination surface information

      subroutine inp_dst_surface (F_p0,F_q,F_u,F_v,F_lsq,F_lsu,F_lsv   ,&
                                  F_sq0, F_sLSq0, &
                                  F_ip1list, F_me, F_tv,F_topo  ,&
                                  F_orols,F_sleve_L,Minx,Maxx,Miny,Maxy,&
                                  F_nka, F_i0, F_in, F_j0, F_jn)
      use cstv
      use gmm_pw
      use gmm_geof
      use rmn_gmm
      implicit none

      logical, intent(in) :: F_sleve_L
      integer, intent(in) :: Minx, Maxx, Miny, Maxy,F_nka,&
                             F_i0,F_in,F_j0,F_jn
      integer, dimension (:) , pointer, intent(inout) :: F_ip1list
      real, dimension(*), intent(in) :: F_tv
      real, dimension(:,:), pointer, intent(in) :: F_sq0,F_sLSq0
      real, dimension(Minx:Maxx,Miny:Maxy), intent(out) :: F_p0
      real, dimension(Minx:Maxx,Miny:Maxy), intent(in) :: F_topo,F_orols
      real   , dimension (:,:), pointer, intent(inout) :: F_q,F_u,F_v,&
                                                    F_lsq,F_lsu,F_lsv
      real   , dimension (:,:,:), pointer , intent(inout) :: F_me

      character(len=4) vname
      integer i,j,typ
      real lev
      real, dimension (:,:,:), pointer :: srclev
      real(kind=REAL64) :: oneoRT
      logical sfcTT_L
!
!---------------------------------------------------------------------
!
      sfcTT_L= .false.
      if ( associated(F_ip1list) ) then
         call convip (F_ip1list(F_nka), lev, typ,-1, vname, .false.)
         sfcTT_L = (typ == 4) .or. &
                 ( (typ /= 2) .and. (abs(lev-1.) <= 1.e-5) )
      endif

      allocate (F_q(l_minx:l_maxx,l_miny:l_maxy), &
                F_u(l_minx:l_maxx,l_miny:l_maxy), &
                F_v(l_minx:l_maxx,l_miny:l_maxy) )

      if ( associated(F_me) .and. sfcTT_L &
                            .and. (Inp_kind /= 21) ) then
         if (lun_out > 0) then
            write(lun_out,'(" PERFORMING surface pressure adjustment")')
         end if
         allocate (srclev(l_minx:l_maxx,l_miny:l_maxy,F_nka))
         F_q= F_sq0
         call inp_3dpres ( Inp_vgd_src,F_ip1list,F_sq0,F_sLSq0,&
                           srclev,1,F_nka )
         srclev(:,:,F_nka)= F_sq0(:,:)
         call adj_ss2topo ( F_p0, F_topo, srclev, F_me, F_tv,&
                           Minx,Maxx,Miny,Maxy,F_nka,F_i0,F_in,F_j0,F_jn)
         deallocate (F_me,srclev) ; nullify (F_me)
      else
         if (lun_out > 0) then
            write(lun_out,'(" NO surface pressure adjustment")')
         end if
         F_p0 = F_sq0
      end if

      if (Inp_dst_hauteur_L) then
            F_q(1-G_halox:l_ni+G_halox,1-G_haloy:l_nj+G_haloy)= &
         F_topo(1-G_halox:l_ni+G_halox,1-G_haloy:l_nj+G_haloy) / grav_8
      else
         F_q = F_p0
      endif

      do j=1-G_haloy,l_nj+G_haloy
         do i=1-G_halox,l_ni+G_halox-1
            F_u(i,j)= (F_q(i,j)+F_q(i+1,j))*.5d0
         end do
      end do
      if (l_east)  then
         F_u(l_ni+G_halox,:)= F_q(l_ni+G_halox,:)
      endif
      do j=1-G_haloy,l_nj+G_haloy-1
         do i=1-G_halox,l_ni+G_halox
            F_v(i,j)= (F_q(i,j)+F_q(i,j+1))*.5d0
         end do
      end do
      if (l_north) F_v(:,l_nj+G_haloy)= F_q(:,l_nj+G_haloy)

      call gem_xch_halo (F_u, l_minx, l_maxx, l_miny, l_maxy, 1)
      call gem_xch_halo (F_v, l_minx, l_maxx, l_miny, l_maxy, 1)

      if ( F_sleve_L ) then
         allocate (F_lsq(l_minx:l_maxx,l_miny:l_maxy), &
                   F_lsu(l_minx:l_maxx,l_miny:l_maxy), &
                   F_lsv(l_minx:l_maxx,l_miny:l_maxy) )
         if (Inp_dst_hauteur_L) then
            F_lsq(1-G_halox:l_ni+G_halox,1-G_haloy:l_nj+G_haloy) = F_orols(1-G_halox:l_ni+G_halox,1-G_haloy:l_nj+G_haloy) / grav_8
         else
            oneoRT= 1.d0 / (rgasd_8 * Tcdk_8)
            F_lsq(1-G_halox:l_ni+G_halox,1-G_haloy:l_nj+G_haloy) = Cstv_pref_8 * exp(-F_orols(1-G_halox:l_ni+G_halox,1-G_haloy:l_nj+G_haloy) * oneoRT)
         endif
         do j=1-G_haloy,l_nj+G_haloy
            do i=1-G_halox,l_ni+G_halox-1
               F_lsu(i,j)= (F_lsq(i,j)+F_lsq(i+1,j))*.5d0
            end do
         end do
         if (l_east)  then
            F_lsu(l_ni+G_halox,:)= F_lsq(l_ni+G_halox,:)
         endif
         do j=1-G_haloy,l_nj+G_haloy-1
            do i=1-G_halox,l_ni+G_halox
               F_lsv(i,j)= (F_lsq(i,j)+F_lsq(i,j+1))*.5d0
            end do
         end do
         if (l_north) F_lsv(:,l_nj+G_haloy)= F_lsq(:,l_nj+G_haloy)

         call gem_xch_halo (F_lsu, l_minx, l_maxx, l_miny, l_maxy, 1)
         call gem_xch_halo (F_lsv, l_minx, l_maxx, l_miny, l_maxy, 1)
      end if
!
!---------------------------------------------------------------------
!
      return
      end subroutine inp_dst_surface

!**s/r inp_read_uv - Read horizontal winds UU,VV F valid at F_datev
!                    and perform vectorial horizontal interpolation
!                    on proper Arakawa grid u and v respectively

      subroutine inp_read_uv ( F_u, F_v, F_target_S, F_ip1, F_nka )
      use rmn_fst24

      implicit none

      character(len=*)                  , intent(in)  :: F_target_S
      integer                           , intent(out) :: F_nka
      integer, dimension(:    ), pointer, intent(inout) :: F_ip1
      real   , dimension(:,:,:), pointer, intent(inout) :: F_u, F_v

!     local variables
      character(len=1) typ
      character(len=4) var
      character(len=12) lab
      integer, parameter :: nlis = 1024
      integer liste_sorted(nlis),lislon
      integer i,err, nz, same_rot, ni_dest, nj_dest
      integer mpx,local_nk,irest,kstart, src_gid, dst_gid, vcode
      integer dstu_gid,dstv_gid,erru,errv
      integer, dimension(:  ), allocatable :: zlist
      real(kind = real32), dimension(:,:), allocatable :: uhr,vhr
      real(kind = real32), dimension(:  ), allocatable, target :: uv,u,v
      real(kind = real32), dimension(:  ), pointer     :: posxu,posyu,posxv,posyv

      type(fst_record) :: recs_u(nlis),recs_v(nlis)
      type(fst_query)  :: query_u,query_v
      logical          :: success

      common /bcast_i_uv / lislon,nz,same_rot
!
!---------------------------------------------------------------------
!
      local_nk= 0 ; F_nka= -1
      if (associated(F_ip1)) deallocate (F_ip1)
      if (associated(F_u  )) deallocate (F_u  )
      if (associated(F_v  )) deallocate (F_v  )
      nullify (F_ip1, F_u, F_v)

      if (Inp_iome >= 0) then
         vcode= -1 ; nz= -1 ; same_rot= -1

         query_u = Inp_file%new_query(datev=Inp_cmcdate,nomvar='UU  ')
         lislon = query_u%find_all(recs_u)

         query_v = Inp_file%new_query(datev=Inp_cmcdate,nomvar='VV  ')
         lislon = query_v%find_all(recs_v)

         if (lislon == 0) goto 999

         src_gid = ezqkdef (recs_u(1)%ni,recs_u(1)%nj,recs_u(1)%grtyp,recs_u(1)%ig1,recs_u(1)%ig2,recs_u(1)%ig3,recs_u(1)%ig4,Inp_file%get_unit())

         i= lislon
         call record_sort_ip1 (recs_u,liste_sorted,i)
         call record_sort_ip1 (recs_v,liste_sorted,lislon)

         allocate (F_ip1(max(1,lislon)))
         if (lislon > 1) then
            F_ip1(1:lislon) = liste_sorted(1:lislon)
         else
            F_ip1(1) = recs_u(1)%ip1
         end if

         nz= (lislon + Inp_npes - 1) / Inp_npes

         mpx      = mod( Inp_iome, Inp_npes )
         local_nk = lislon / Inp_npes
         irest  = lislon  - local_nk * Inp_npes
         kstart = mpx * local_nk + 1
         if ( mpx < irest ) then
            local_nk   = local_nk + 1
            kstart = kstart + mpx
         else
            kstart = kstart + irest
         end if

         ni_dest= G_ni+2*G_halox
         nj_dest= G_nj+2*G_haloy
         allocate (uhr(ni_dest*nj_dest,nz),vhr(ni_dest*nj_dest,nz))
         allocate (u(recs_u(1)%ni*recs_u(1)%nj), v(recs_u(1)%ni*recs_u(1)%nj))
         allocate (uv(ni_dest*nj_dest))

         if (local_nk > 0) then

            err = ezsetopt ('INTERP_DEGREE', 'CUBIC')


            if (trim(F_target_S) == 'UV') then

               posxu => geomh_lonF
               posyu => geomh_latQ
               posxv => geomh_lonQ
               posyv => geomh_latF

               write(output_unit,1001) 'Interpolating: UU, nka= ',&
                             lislon,', valid: ',Inp_datev,' on U grid'
               dstu_gid = ezgdef_fmem ( ni_dest, nj_dest, 'Z', 'E', &
                       Hgc_ig1ro, Hgc_ig2ro, Hgc_ig3ro, Hgc_ig4ro, &
                                                      posxu, posyu )
               write(output_unit,1001) 'Interpolating: VV, nka= ',&
                             lislon,', valid: ',Inp_datev,' on V grid'
               dstv_gid = ezgdef_fmem ( ni_dest, nj_dest, 'Z', 'E', &
                       Hgc_ig1ro, Hgc_ig2ro, Hgc_ig3ro, Hgc_ig4ro, &
                                                      posxv, posyv )
               do i=1,local_nk
                  success = recs_u(kstart+i-1)%read(data=c_loc(u))
                  success = recs_v(kstart+i-1)%read(data=c_loc(v))
 
                  err = ezdefset ( dstu_gid , src_gid )
                  erru = ezuvint  ( uhr(1,i),uv, u,v )
                  err = ezdefset ( dstv_gid , src_gid )
                  errv = ezuvint  ( uv,vhr(1,i), u,v )
               end do
               if (erru == 2) &
                      write(output_unit,1002) 'EXTRApolating: UU, nka= ',&
                             lislon,', valid: ',Inp_datev,' on U grid'
               if (errv == 2) &
                      write(output_unit,1002) 'EXTRApolating: VV, nka= ',&
                             lislon,', valid: ',Inp_datev,' on V grid'

            else !Q

               posxu => geomh_lonQ
               posyu => geomh_latQ

               write(output_unit,1001) 'Interpolating: UV, nka= ',&
                              lislon,', valid: ',Inp_datev,' on Q grid'
               dst_gid = ezgdef_fmem ( ni_dest, nj_dest, 'Z', 'E', &
                       Hgc_ig1ro, Hgc_ig2ro, Hgc_ig3ro, Hgc_ig4ro, &
                                                      posxu, posyu )
               err = ezdefset ( dst_gid , src_gid )
               do i=1,local_nk
                  success = recs_u(kstart+i-1)%read(data=c_loc(u))
                  success = recs_v(kstart+i-1)%read(data=c_loc(v))

                  err = ezuvint  ( uhr(1,i),vhr(1,i), u,v )
               end do
               if (err == 2) &
                 write(output_unit,1002) 'EXTRApolating: UV, nka= ',&
                             lislon,', valid: ',Inp_datev,' on Q grid'

            end if


         end if
         deallocate (uv,u,v)

      else
         allocate (uhr(1,1), vhr(1,1))
      end if

 999  call rpn_comm_bcast (lislon, 3, "MPI_INTEGER", Inp_iobcast, &
                        "grid", err) !NOTE: bcast lislon AND nz, samerot
      F_nka= lislon

      if (F_nka > 0) then

         if (F_nka >= 1) then
            if (Inp_iome < 0) allocate ( F_ip1(F_nka) )
            call rpn_comm_bcast ( F_ip1, F_nka, "MPI_INTEGER", &
                                  Inp_iobcast, "grid", err )
         end if

         allocate (zlist(nz)) ; zlist= -1
         do i=1, local_nk
            zlist(i)= i + kstart - 1
         end do

         allocate ( F_u(l_minx:l_maxx,l_miny:l_maxy,lislon),&
                    F_v(l_minx:l_maxx,l_miny:l_maxy,lislon) )

         call glbdist_os (uhr,F_u, l_minx,l_maxx,l_miny,l_maxy,F_nka,&
                      G_ni+G_halox,G_nj+G_haloy,zlist,nz,knams_8,0.d0)

         call glbdist_os (vhr,F_v, l_minx,l_maxx,l_miny,l_maxy,F_nka,&
                      G_ni+G_halox,G_nj+G_haloy,zlist,nz,knams_8,0.d0)
         deallocate (uhr,vhr,zlist)

      else

         if (Inp_iome >= 0) write(output_unit,'(3a)') &
                  'Variable: UU,VV valid: ',Inp_datev, 'NOT FOUND'
         call gem_error ( -1, 'inp_read_uv', &
                  'Missing input data: horizontal winds')

      end if

 1001 format (a,i3,3a)
 1002 format (a,i3,3a)
!
!---------------------------------------------------------------------
!
      return
      end subroutine inp_read_uv

      subroutine inp_3dpres ( F_vgd, F_ip1, F_sfc, F_sfcL, &
                              F_dest, k0, kn, F_inlog_S)

      character(len=*), optional, intent(in) :: F_inlog_S
      type(vgrid_descriptor) , intent(in) :: F_vgd
      integer                , intent(in) :: k0,kn
      integer, dimension(:)    , pointer, intent(in ) :: F_ip1
      real   , dimension(:,:  ), pointer, intent(in ) :: F_sfc,F_sfcL
      real   , dimension(:,:,:), pointer, intent(inout) :: F_dest

!     local variables
      integer istat, nka, kind,i,j,k
      logical inlog, vgd
      real, dimension (:,:  ), pointer :: pres,presl
      real, dimension (:,:,:), pointer :: ptr3d
!
!---------------------------------------------------------------------
!
      inlog= .false.
      if ( present(F_inlog_S) ) then
         if (F_inlog_S == 'in_log') inlog= .true.
      end if

      pres  => F_sfc (1-G_halox:l_ni+G_halox,1-G_haloy:l_nj+G_haloy)
      ptr3d => F_dest(1-G_halox:l_ni+G_halox,1-G_haloy:l_nj+G_haloy,k0:kn)
      vgd = vgd_get ( F_vgd, key='KIND',value=kind,quiet=.true. )==VGD_OK

      if (vgd) then
         if ( associated (F_sfcL) ) then
            presl=> F_sfcL (1-G_halox:l_ni+G_halox,1-G_haloy:l_nj+G_haloy)
            istat= vgd_levels (F_vgd,F_ip1(k0:kn),ptr3d,sfc_field=pres,&
                               sfc_field_ls=presl, in_log=inlog)
         else
           istat= vgd_levels (F_vgd,F_ip1(k0:kn),ptr3d, pres, &
                               in_log=inlog)
         end if

      else

         nka= kn-k0+1
         istat= inp_match (F_dest, PX3d%valq, F_ip1, PX3d%ip1, &
                       l_minx,l_maxx,l_miny,l_maxy,nka, PX3d%nk)
        if (nka /= kn-k0+1) then
            istat=-1
         else
            if (inlog) then
            do k=1,nka
               do j= 1-G_haloy, l_nj+G_haloy
                  do i= 1-G_halox, l_ni+G_halox
                     F_dest(i,j,k) = log(F_dest(i,j,k))
                  end do
               end do
            end do
            endif
         endif

      endif

      call handle_error_l(istat==VGD_OK,'inp_3dpres',&
                          'Error computing pressure')
!
!---------------------------------------------------------------------
!
      return
      end subroutine inp_3dpres

      subroutine inp_src_levels ( F_dest, F_nk, F_ip1, F_vgd, F_sfc,&
                         F_sfcL, F_gz, F_gz_ip1, Minx,Maxx,Miny,Maxy)

      use dynkernel_options

      implicit none

      type(vgrid_descriptor) , intent(in) :: F_vgd
      integer                , intent(inout) :: F_nk
      integer                , intent(in)  :: Minx,Maxx,Miny,Maxy
      integer, dimension(:)    , pointer, intent(in) :: F_ip1,F_gz_ip1
      real   , dimension(:,:  ), pointer, intent(in) :: F_sfc,F_sfcL
      real   , dimension(:,:,:), pointer, intent(in) :: F_gz
      real   , dimension(:,:,:), pointer, intent(inout) :: F_dest

!     local variables
      integer istat
!
!---------------------------------------------------------------------
!
      F_nk= ubound(F_ip1,1)
      if (Inp_src_hauteur_L.and..not.Schm_autobar_L) then
         istat= inp_match (F_dest, F_gz, F_ip1, F_gz_ip1, &
                       Minx,Maxx,Miny,Maxy,F_nk, ubound(F_gz_ip1,1))
      else
         call inp_3dpres ( F_vgd,  F_ip1,  F_sfc,  F_sfcL, F_dest, &
                           1, F_nk, F_inlog_S='in_log' )
      endif
!
!---------------------------------------------------------------------
!
      return
      end subroutine inp_src_levels

      subroutine inp_dst_levels ( F_dest, F_vgd, F_ip1, &
                                  F_sfc, F_sfcL, F_nk )
      implicit none

      type(vgrid_descriptor) , intent(in) :: F_vgd
      integer, intent(in ) :: F_nk
      integer, dimension(:)    , pointer, intent(in ) :: F_ip1
      real   , dimension(:,:  ), pointer, intent(in ) :: F_sfc,F_sfcL
      real   , dimension(:,:,:), pointer, intent(inout) :: F_dest
!
!---------------------------------------------------------------------
!
      if (Inp_dst_hauteur_L) then
         call inp_3dhgts ( F_vgd, F_ip1, F_sfc, F_sfcL, F_dest, &
                           1, F_nk)
      else
         call inp_3dpres ( F_vgd, F_ip1, F_sfc, F_sfcL, F_dest, &
                           1, F_nk, F_inlog_S='in_log' )
      endif
!
!---------------------------------------------------------------------
!
      return
      end subroutine inp_dst_levels
!
!---------------------------------------------------------------------
!
!**s/r inp_match - To match heights
!
      integer function inp_match ( F_ho, F_hi, F_ip1_list_o, &
            F_ip1_list_i, Minx,Maxx,Miny,Maxy, nko, nki) result(status)

      integer, intent(inout) :: nko
      integer, intent(in)    :: Minx,Maxx,Miny,Maxy,nki
      integer, dimension (:), pointer, intent(in) :: &
                             F_ip1_list_i,F_ip1_list_o
      real, dimension (:,:,:), pointer, intent(inout) :: F_ho
      real, dimension (Minx:Maxx,Miny:Maxy,nki), intent(in) :: F_hi
      integer, parameter :: INP_OK = 0, INP_ERROR = -1

!     Local variables
      integer :: index_diag_AGL, ko, ki, kind, i, j, k
      logical :: found_L
      real :: lev
      character(len=128) :: message_S
!
!---------------------------------------------------------------------
!
      status = INP_ERROR

      index_diag_AGL = -1

      do ko=1, nko
         found_L = .false.
         call convip (F_ip1_list_o(ko), lev, kind, -1, message_S,.false.)
         if( kind == 4 )then
            found_L = .true.
            ! Pour le moment sort_ip1 ne garde que le dernier kind = 4
            ! de la liste. Donc cette erreur ne se produira pas.
            if( index_diag_AGL /= -1 )&
                 call gem_error (-1,'inp_match',&
                        'more than one diagnostic level, review code')
            index_diag_AGL = ko
            F_ho(1-G_halox:l_ni+G_halox,1-G_haloy:l_nj+G_haloy,ko ) = &
            F_hi(1-G_halox:l_ni+G_halox,1-G_haloy:l_nj+G_haloy,nki) + lev
         else
            do ki=1, nki
               if (F_ip1_list_o(ko) == F_ip1_list_i(ki)) then
                   found_L = .true.
                   F_ho(:,:,ko) = F_hi(:,:,ki)
                   exit
               end if
            end do
         end if
         if(.not. found_L)then
            write(message_S,*)'Missing field: GZ for ip1 = ',&
                 F_ip1_list_o(ko)
            call gem_error (-1,'inp_match',message_S)
         end if
      end do

      if(index_diag_AGL /= -1)then
         ! Check monotonicity
         found_L = .false.
         do j = 1, l_nj
            do i = 1, l_ni
               if( F_ho(i,j, index_diag_AGL-1) <= &
                   F_ho(i,j, index_diag_AGL) ) then
                  found_L = .true.
               end if
            end do
         end do
         if(index_diag_AGL < nko)then
            ! Check level below
            do j = 1, l_nj
               do i = 1, l_ni
                  if( F_ho(i,j, index_diag_AGL+1) >= &
                      F_ho(i,j, index_diag_AGL) ) then
                     found_L = .true.
                  end if
               end do
            end do
         end if
         if(found_L)then
            ! Remove diag level
            do k = index_diag_AGL, nko-1
               F_ho(:,:,k) = F_ho(:,:,k+1)
            end do
            nko = nko - 1
         end if
      end if

      status = INP_OK
!
!---------------------------------------------------------------------
!
      return
      end function inp_match

      subroutine inp_3dhgts ( F_vgd, F_ip1, F_sfc, F_sfcL, F_dest,k0,kn )

      type(vgrid_descriptor) , intent(IN) :: F_vgd
      integer                , intent(IN) :: k0,kn
      integer, dimension(:)    , pointer, intent(in ) :: F_ip1
      real   , dimension(:,:  ), pointer, intent(in ) :: F_sfc,F_sfcL
      real   , dimension(:,:,:), pointer, intent(inout) :: F_dest

!     local variables
      integer :: istat
      real, dimension (:,:  ), pointer :: hgt,hgtls
      real, dimension (:,:,:), pointer :: ptr3d
!
!---------------------------------------------------------------------
!
      hgt   => F_sfc (1-G_halox:l_ni+G_halox,1-G_haloy:l_nj+G_haloy)
      ptr3d => F_dest(1-G_halox:l_ni+G_halox,1-G_haloy:l_nj+G_haloy,k0:kn)
      if (associated(F_sfcL)) then
         hgtls => F_sfcL(1-G_halox:l_ni+G_halox,1-G_haloy:l_nj+G_haloy)
      else
         nullify (hgtls)
      endif
      istat= vgd_levels ( F_vgd, F_ip1(k0:kn), ptr3d, &
                          sfc_field=hgt, sfc_field_ls=hgtls )
!
!---------------------------------------------------------------------
!
      return
      end subroutine inp_3dhgts

!**s/r inp_hwnd - Read horizontal winds UU,VV F valid at F_datev
!                 and perform vectorial horizontal interpolation
!                 and  vertical interpolation to momentum levels

      subroutine inp_hwnd ( F_u,F_v, F_stag_L    ,&
                      F_ssqr,F_ssur,F_ssvr, F_ssqrLS,F_ssurLS,F_ssvrLS,&
                      F_ssq0,F_ssu0,F_ssv0, F_ssq0LS,F_ssu0LS,F_ssv0LS,&
                            F_gz_q, F_gz_u, F_gz_v, F_GZ_ip1          ,&
                            Minx,Maxx,Miny,Maxy, F_nk )
      implicit none

      logical                , intent(in)  :: F_stag_L
      integer                , intent(in)  :: Minx,Maxx,Miny,Maxy, F_nk
      integer, dimension(:)  , pointer, intent(in)  :: F_GZ_ip1
      real, dimension (:,:), pointer, intent(in) :: &
                             F_ssqr,F_ssur,F_ssvr, F_ssqrLS,F_ssurLS,F_ssvrLS,&
                             F_ssq0,F_ssu0,F_ssv0, F_ssq0LS,F_ssu0LS,F_ssv0LS
      real, dimension (:,:,:), pointer, intent(in) :: F_gz_q, F_gz_u, &
                                                      F_gz_v
      real, dimension(Minx:Maxx,Miny:Maxy,F_nk), intent(out) :: F_u, F_v

!     local variables
      integer nka
      integer, contiguous, dimension (:    ), pointer :: ip1_list, ip1_target
      real   , contiguous, dimension (:,:,:), pointer :: srclev,dstlev
      real   , contiguous, dimension (:,:,:), pointer :: ur,vr
!
!---------------------------------------------------------------------
!
      nullify (ip1_list, ur, vr)
      if (F_stag_L) then
         call inp_read_uv ( ur, vr, 'UV' , ip1_list, nka )
      else
         call inp_read_uv ( ur, vr, 'Q ' , ip1_list, nka )
      end if

      allocate ( srclev(l_minx:l_maxx,l_miny:l_maxy,nka) ,&
                 dstlev(l_minx:l_maxx,l_miny:l_maxy,G_nk) )
      allocate (ip1_target(1:G_nk))
      ip1_target(1:G_nk)= Ver_ip1%m(1:G_nk)
      if (F_stag_L) then

         call inp_src_levels (srclev, nka, ip1_list, Inp_vgd_src, F_ssur,&
                            F_ssurLS, F_gz_u, F_GZ_ip1,Minx,Maxx,Miny,Maxy)


         call inp_dst_levels (dstlev, Ver_vgdobj, ip1_target,&
                              F_ssu0, F_ssu0Ls, G_nk)

         call vertint2 ( F_u,dstlev,G_nk, ur,srclev,nka,&
                         l_minx,l_maxx,l_miny,l_maxy   ,&
                         1-G_halox,l_ni+G_halox, 1-G_haloy,l_nj+G_haloy,&
                         levtype=Inp_levtype_S )

         call inp_src_levels (srclev, nka, ip1_list, Inp_vgd_src, F_ssvr,&
                            F_ssvrLS, F_gz_v, F_GZ_ip1,Minx,Maxx,Miny,Maxy)

         call inp_dst_levels (dstlev, Ver_vgdobj, ip1_target,&
                              F_ssv0, F_ssv0Ls, G_nk)

         call vertint2 ( F_v,dstlev,G_nk, vr,srclev,nka,&
                         l_minx,l_maxx,l_miny,l_maxy   ,&
                         1-G_halox,l_ni+G_halox, 1-G_haloy,l_nj+G_haloy,&
                         levtype=Inp_levtype_S )

      else

         call inp_src_levels (srclev, nka, ip1_list, Inp_vgd_src, F_ssqr,&
                            F_ssqrLS, F_gz_q, F_GZ_ip1,Minx,Maxx,Miny,Maxy)

         call inp_dst_levels (dstlev, Ver_vgdobj, ip1_target,&
                              F_ssq0, F_ssq0Ls, G_nk)

         call vertint2 ( F_u,dstlev,G_nk, ur,srclev,nka,&
                         l_minx,l_maxx,l_miny,l_maxy   ,&
                         1-G_halox,l_ni+G_halox, 1-G_haloy,l_nj+G_haloy,&
                         levtype=Inp_levtype_S )
         call vertint2 ( F_v,dstlev,G_nk, vr,srclev,nka,&
                         l_minx,l_maxx,l_miny,l_maxy   ,&
                         1-G_halox,l_ni+G_halox, 1-G_haloy,l_nj+G_haloy,&
                         levtype=Inp_levtype_S)

      end if

      deallocate (ip1_list,ip1_target,ur,vr,srclev,dstlev)
      nullify    (ip1_list,ip1_target,ur,vr,srclev,dstlev)
!
!---------------------------------------------------------------------
!
      return
      end subroutine inp_hwnd

      logical function inp_src_vert3d (F_var_S, F_3dv)
      implicit none

      character(len=*), intent(IN) :: F_var_S
      character(len=1) :: grid_S(3)
      type(Inp_vrt3d), intent(inout) :: F_3dv

      character(len=4  ) :: dumc
      integer, parameter :: nlis = 1024
      integer err_mt, lislon
      real level
      real, dimension(:,:,:), pointer :: val
!!$ 0 	KIND_ABOVE_SEA 	height (m) above mean sea level 	(-20,000 -> 100,000)
!!$ 1 	KIND_SIGMA 	sigma coordinates 	(0.0 -> 1.0)
!!$ 2 	KIND_PRESSURE 	pressure (mb) 	(0 -> 1100)
!!$ 3 	KIND_ARBITRARY 	arbitrary number, no units 	(-4.8e8 -> 1.0e10)
!!$ 4 	KIND_ABOVE_GND 	height (m) above ground 	(-20,000 -> 100,000)
!!$ 5 	KIND_HYBRID 	hybrid coordinates 	(0.0 -> 1.0)
!!$ 6 	KIND_THETA 	theta coordinates 	(1 -> 200,000)
!!$10 	KIND_HOURS 	time (hours) 	(0.0 -> 1.0e10)
!!$15 	KIND_SAMPLES 	reserved (integer value) 	(0 -> 1 999 999)
!!$17 	KIND_MTX_IND 	conversion matrix x subscript
!!$                     (shared with kind=1) (1.0 -> 1.0e10)
!!$21 	KIND_M_PRES 	pressure-meters (shared with kind=5)
!!$                   	(0 -> 1,000,000) fact=1E+4 
!
!---------------------------------------------------------------------
!
      F_3dv%nk= -1 ; F_3dv%kind= -1
      if (associated(F_3dv%valq)) deallocate (F_3dv%valq)
      if (associated(F_3dv%valu)) deallocate (F_3dv%valu)
      if (associated(F_3dv%valv)) deallocate (F_3dv%valv)
      nullify (F_3dv%valq,F_3dv%valu,F_3dv%valv,val)
      inp_src_vert3d = .false.
      grid_S=['Q','U','V']

      err_mt = inp_read_mt (F_var_S,grid_S, val,3,&
                         F_3dv%ip1, lislon,F_quiet_L=.true.)
      inp_src_vert3d= ((err_mt == 0) .and. (lislon>1))

      if (inp_src_vert3d) then
!!$do k=1,lislon
!!$         call convip (F_3dv%ip1(k), level, &
!!$                      F_3dv%kind, -1,dumc,.false. )
!!$print*, k,F_3dv%ip1(k), level,F_3dv%kind
!!$                   end do
!!$call gem_error(-1,'','')
         F_3dv%nk= lislon
         call convip (F_3dv%ip1(F_3dv%nk), level, &
                      F_3dv%kind, -1,dumc,.false. )
         if (Lun_out > 0) write(lun_out,9000) F_3dv%nk,F_var_S,F_3dv%kind
 9000 format(x,i3,' levels of ',a,' provided in input file with kind=',i5)

!!$         limit_near_sfc= abs(level-1.)
!!$         if (F_3dv%kind==21) limit_near_sfc=level  
!!$         if (F_3dv%kind==2 ) limit_near_sfc=0.
!!$         if( limit_near_sfc > epsilon(level)) then
!!$            err= -1
!!$            message_S= 'Missing near surface '//trim(F_var_S)
!!$         else
            allocate (F_3dv%valq(l_minx:l_maxx,l_miny:l_maxy,1:F_3dv%nk))
            allocate (F_3dv%valu(l_minx:l_maxx,l_miny:l_maxy,1:F_3dv%nk))
            allocate (F_3dv%valv(l_minx:l_maxx,l_miny:l_maxy,1:F_3dv%nk))
            F_3dv%valq(l_minx:l_maxx,l_miny:l_maxy,1:F_3dv%nk)= val(l_minx:l_maxx,l_miny:l_maxy,           1:  F_3dv%nk)
            F_3dv%valu(l_minx:l_maxx,l_miny:l_maxy,1:F_3dv%nk)= val(l_minx:l_maxx,l_miny:l_maxy,  F_3dv%nk+1:2*F_3dv%nk)
            F_3dv%valv(l_minx:l_maxx,l_miny:l_maxy,1:F_3dv%nk)= val(l_minx:l_maxx,l_miny:l_maxy,2*F_3dv%nk+1:3*F_3dv%nk)
!         endif
!         deallocate (val)
      endif
      if (associated(val)) deallocate (val)
      nullify (val)
!      call gem_error ( err,'inp_src_vert3d', message_S)
!
!---------------------------------------------------------------------
!
       return
       end function inp_src_vert3d

end module inp_base
      subroutine reshape_wk (src,dst,nix,njx,ni,nj,hx,hy)
      use, intrinsic :: iso_fortran_env
      implicit none
      integer nix,njx,ni,nj,hx,hy
      real src(nix,njx), dst(ni,nj)
      integer i,j
      do j=1,nj
      do i=1,ni
         dst(i,j) = src(i+hx,j+hy)
      end do
      end do
      return
      end
