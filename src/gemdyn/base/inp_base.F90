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

      integer function inp_get ( F_var_S, F_hgrid_S, F_ver_ip1       ,&
                         F_sfc_src, F_sfc_dst, F_sfcLS_dst           ,&
                         F_gz, F_GZ_ip1, F_dest , Minx,Maxx,Miny,Maxy,&
                         F_nk, F_inttype_S, F_quiet_L )
      character(len=*)          , intent(in) :: F_var_S,F_hgrid_S
      character(len=*), optional, intent(in) :: F_inttype_S
      logical         , optional, intent(in) :: F_quiet_L
      integer                   , intent(in) :: Minx,Maxx,Miny,Maxy, F_nk
      integer, dimension(:)  , pointer, intent(in) :: F_ver_ip1,F_gz_ip1
      real, dimension (:,:), pointer, intent(in) :: &
                                      F_sfc_src,F_sfc_dst,F_sfcLS_dst
      real, dimension(Minx:Maxx,Miny:Maxy,F_nk), intent(out):: F_dest
      real, dimension(:,:,:), pointer, intent(in ) :: F_gz

!     local variables
      character(len=12) :: inttype
      logical quiet_L
      integer nka
      integer, dimension (:    ), pointer :: ip1_list
      real   , dimension (:,:  ), pointer :: dummy
      real   , dimension (:,:,:), pointer :: wrkr,srclev,dstlev
!
!---------------------------------------------------------------------
!
      inp_get= -1
      if ( any (Inp_blacklist_S(1:MAX_blacklist) == trim(F_var_S)) ) &
           return

      nullify (ip1_list, wrkr, dummy)
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
                 dstlev(l_minx:l_maxx,l_miny:l_maxy,G_nk) )

      call inp_src_levels (srclev, nka, ip1_list, Inp_vgd_src, &
         F_sfc_src, dummy, F_gz, F_GZ_ip1,Minx,Maxx,Miny,Maxy)

      call inp_dst_levels ( dstlev, Ver_vgdobj, F_ver_ip1,&
                            F_sfc_dst, F_sfcLS_dst )

      inttype= 'cubic'
      if (present(F_inttype_S)) inttype= F_inttype_S
      call vertint2 ( F_dest,dstlev,G_nk, wrkr,srclev,nka       ,&
                      l_minx,l_maxx,l_miny,l_maxy, 1,l_ni,1,l_nj,&
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
      implicit none

      character(len=*)          ,intent(in)  :: F_var_S
      character(len=*), dimension(*),intent(in)  :: F_hgrid_S
      character(len=*), optional,intent(in)  :: F_hint_S
      logical         , optional,intent(in)  :: F_quiet_L
      integer                   ,intent(in ) :: F_nd
      integer                   ,intent(out) :: F_nka
      integer, dimension(:    ), pointer,intent  (out) :: F_ip1
      real   , dimension(:,:,:), pointer,intent(inout) :: F_dest

!     local variables
      integer, external :: RPN_COMM_shuf_ezdist, &
                           samegrid_gid, samegrid_rot, inp_is_real_wind
      character(len=1) typ,grd
      character(len=4) nomvar,var,dumc
      character(len=12) lab,interp_S
      logical, dimension (:), allocatable :: zlist_o
      logical :: quiet_L
      integer, parameter :: nlis = 1024
      integer i, idst, err, nz, n1,n2,n3, nrec, liste(nlis),&
              liste_sorted(nlis),lislon,cnt
      ! Remove the following line by 2021
      integer ni_dest,nj_dest,ut1_is_urt1
      integer subid,nicore,njcore,datev
      integer mpx,local_nk,irest,kstart, src_gid, vcode, ip1
      integer dte, det, ipas, p1, p2, p3, g1, g2, g3, g4, bit, &
              dty, swa, lng, dlf, ubc, ex1, ex2, ex3
      integer, dimension(:  ), allocatable :: zlist
      real :: surface_level
      real   , dimension(:  ), allocatable :: wk3
      real   , dimension(:,:), allocatable :: wk1,wk2
      real   , dimension(:  ), pointer     :: posx,posy
      real(kind=REAL64) add, mult
      common /bcast_i / lislon,nz
!
!---------------------------------------------------------------------
!
      inp_read_mt= -1
      F_nka= -1 ; local_nk= 0
      add= 0.d0 ; mult= 1.d0
      ! Remove the following line by 2021
      ut1_is_urt1 = -1
      if (associated(F_ip1 )) deallocate (F_ip1 )
      if (associated(F_dest)) deallocate (F_dest)
      nullify (F_ip1,F_dest)
      quiet_L=.false.
      if (present(F_quiet_L)) quiet_L= F_quiet_L

      nomvar = F_var_S ; ip1= -1
      select case (F_var_S)
         case ('OROGRAPHY')
            if (Inp_kind == 2  ) nomvar= '@NUL'
            if (Inp_kind == 1 ) then
               nomvar= 'GZ' ; ip1= 12000
            endif
            if (Inp_kind == 5 ) then
               nomvar= 'GZ' ; surface_level= 1.
               call convip ( ip1, surface_level,Inp_kind,1,dumc,.false. )
            endif
            if ( Inp_src_hauteur_L ) then
               nomvar= 'GZ' ; surface_level= 0.
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
            if (Inp_kind == 2  ) nomvar= 'GZ'
            if (Inp_kind == 1  ) nomvar= '@NUL'
            if (Inp_kind == 5  ) nomvar= 'GZ'
            if (Inp_src_hauteur_L ) nomvar= 'GZ'
            if ( nomvar == 'GZ' ) mult= 10.d0
         case ('UU')
            mult= knams_8
         case ('VV')
            mult= knams_8
      end select

      datev= Inp_cmcdate
      if ( F_var_S(1:3) == 'TR/' ) then
         nomvar= F_var_S(4:)
         if (Tr3d_anydate_L) datev= -1
      end if

      if ( nomvar == '@NUL' ) return

      if (Inp_iome >= 0) then
         vcode= -1 ; nz= -1
         nrec= fstinl (Inp_handle, n1,n2,n3, datev,' ', &
                       ip1,-1,-1,' ', nomvar,liste,lislon,nlis)
         if (lislon == 0) goto 999

         err= fstprm (liste(1), DTE, DET, IPAS, n1, n2, n3,&
                  BIT, DTY, P1, P2, P3, TYP, VAR, LAB, GRD,&
                  G1,G2,G3,G4,SWA,LNG,DLF,UBC,EX1,EX2,EX3)

         src_gid= ezqkdef (n1, n2, GRD, g1, g2, g3, g4, Inp_handle)

         if ((trim(nomvar) == 'URT1').or.(trim(nomvar) == 'VRT1').or.&
             (trim(nomvar) == 'UT1' ).or.(trim(nomvar) == 'VT1' )) then
             err= samegrid_rot ( src_gid, &
                        Hgc_ig1ro,Hgc_ig2ro,Hgc_ig3ro,Hgc_ig4ro)
             if (err < 0) then
                lislon= 0
                goto 999
             end if
         end if

         call sort_ip1 (liste,liste_sorted,lislon)

         allocate (F_ip1(max(1,lislon)))
         if (lislon > 1) then
            F_ip1(1:lislon) = liste_sorted(1:lislon)
         else
            F_ip1(1) = p1
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
         allocate (wk3(ni_dest*nj_dest))
         allocate (wk2(G_ni*G_nj,(nz+1)*F_nd))
         allocate (wk1(n1*n2,max(local_nk,1)))

         cnt= 0
         do i= kstart, kstart+local_nk-1
            cnt= cnt+1
            err= fstluk (wk1(1,cnt), liste(i), n1,n2,n3)
            ! Remove the following line by 2021
            ut1_is_urt1 = inp_is_real_wind (wk1(1,cnt),n1*n2,nomvar)
         end do

         interp_S= 'CUBIC'
         if (present(F_hint_S)) interp_S= F_hint_S

         do idst= 1, F_nd
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

            if ( GRD == 'U' ) then
              nicore = G_ni-Glb_pil_w-Glb_pil_e
              njcore = G_nj-Glb_pil_s-Glb_pil_n
              if (n1 >= nicore .and. n2/2 >= njcore) then
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
            write(output_unit,1001) 'Interpolating: ',trim(nomvar),', nka= ',&
               lislon,',valid: ',Inp_datev,' on ',F_hgrid_S(idst),' grid'
         end if

         do i=1,local_nk
            err = ezsint(wk3, wk1(1,i))
            call reshape_wk ( wk3,wk2(1,(idst-1)*(nz+1)+i),&
                 ni_dest,nj_dest,G_ni,G_nj,G_halox,G_haloy )
         end do
         if (err == 2) then
            write(output_unit,1002) 'EXTRApolating: ',trim(nomvar),', nka= ',&
               lislon,',valid: ',Inp_datev,' on ',F_hgrid_S(idst),' grid'
         end if
         err = ezsetopt ( 'USE_1SUBGRID', 'NO' )
         end do
         deallocate (wk1,wk3)
      else
         allocate (wk2(1,1))
      end if

 999  call rpn_comm_bcast ( lislon, 2, "MPI_INTEGER", Inp_iobcast, &
                            "grid", err ) !NOTE: bcas listlon AND nz
      F_nka= lislon
      ! Remove the following line by 2021
      call rpn_comm_allreduce ( ut1_is_urt1, Inp_ut1_is_urt1, 1, &
                                "MPI_INTEGER", "MPI_MAX", "grid", err )

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

         allocate ( F_dest(l_minx:l_maxx,l_miny:l_maxy,lislon*F_nd), &
                    zlist_o(lislon) )

         zlist_o= .FALSE.

         do idst=1, F_nd
         zlist_o= .FALSE.
         err = RPN_COMM_shuf_ezdist ( Inp_comm_setno, Inp_comm_id ,&
                                      wk2(1,(idst-1)*(nz+1)+1), nz,&
                           F_dest(l_minx,l_miny,(idst-1)*lislon+1),&
                                            lislon, zlist, zlist_o )
         end do
         deallocate (wk2,zlist,zlist_o)

         F_dest(1:l_ni,1:l_nj,:) = F_dest(1:l_ni,1:l_nj,:) * mult + add
         if (nomvar == 'ST1') &
         F_dest(1:l_ni,1:l_nj,:)= Inp_pref_a_8 * &
                                  exp(F_dest(1:l_ni,1:l_nj,:))

      else

         inp_read_mt= -1
         if ((Inp_iome >= 0).and.(.not.quiet_L)) write(output_unit,'(7a)') &
            ' FIELD: ',trim(F_var_S),':',trim(nomvar),' valid: ',&
            Inp_datev, 'NOT FOUND'

      end if
 1001 format (3a,i3,5a)
 1002 format (3a,i3,5a)
!
!---------------------------------------------------------------------
!
      return
      end function inp_read_mt

!**s/r inp_oro - Read orography from input dataset valid at F_datev

      subroutine inp_oro ( F_topo,F_meqr,F_datev,Minx, Maxx, Miny, Maxy)
      use gmm_geof
      use gem_options
      use glb_ld
      use gmm_itf_mod
      use HORgrid_options
      use lam_options
      use nest_blending, only: nest_blend
      use step_options
      implicit none

      character(len=*), intent(in) :: F_datev
      integer         , intent(in) :: Minx, Maxx, Miny, Maxy
      real, dimension(Minx:Maxx,Miny:Maxy), intent(out) :: F_topo
      real, dimension (:,:), pointer    , intent(inout) :: F_meqr

      integer err,nka
      integer, dimension (:), pointer :: ip1_list
      real, dimension (:,:,:), pointer :: wrk
      real topo_temp(l_minx:l_maxx,l_miny:l_maxy),step_current
      real(kind=REAL64) diffd
!
!---------------------------------------------------------------------
!
      if (associated(F_meqr)) deallocate (F_meqr)
      nullify (F_meqr, wrk, ip1_list)

      err = inp_read_mt ( 'OROGRAPHY', 'Q', wrk, 1, ip1_list, nka )

      if ( associated(ip1_list) ) then
         deallocate (ip1_list) ; nullify (ip1_list)
      end if
      if ( associated(wrk) ) then
         allocate (F_meqr(l_minx:l_maxx,l_miny:l_maxy))
         F_meqr(:,:) = wrk(:,:,1)
         deallocate (wrk) ; nullify (wrk)
      endif

      if ( trim(F_datev) == trim(Step_runstrt_S) ) then
         if ( associated(F_meqr) ) then
            err= gmm_get(gmmk_topo_low_s , topo_low )
            topo_low(1:l_ni,1:l_nj)= F_meqr(1:l_ni,1:l_nj)
         else
            Vtopo_L= .false.
         end if
      end if

      call difdatsd (diffd,Step_runstrt_S,F_datev)
      step_current = diffd*86400.d0 / Step_dt + Step_initial
      call var_topo2 ( F_topo, step_current, &
                       l_minx,l_maxx,l_miny,l_maxy )

      if ( associated(F_meqr) .and. .not. Grd_yinyang_L) then
      if ( Lam_blendoro_L ) then
         topo_temp(1:l_ni,1:l_nj)= F_meqr(1:l_ni,1:l_nj)
         call nest_blend ( F_topo, topo_temp, l_minx,l_maxx, &
                           l_miny,l_maxy, 'M', level=G_nk+1 )
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

      if (Inp_kind /= 105) then ! Not in self cascade mode
         do kt=1, F_nka_tt
inner:      do kh=1, F_nka_hu
               if (F_TT_ip1(kt) == F_HU_ip1(kh)) then
                  call mfottv2 ( &
                          F_tv(l_minx,l_miny,kt),F_tv(l_minx,l_miny,kt),&
                          F_hu(l_minx,l_miny,kh),l_minx,l_maxx         ,&
                                 l_miny,l_maxy,1, 1,l_ni,1,l_nj, .true. )
                  exit inner
               end if
            end do inner
         end do
      end if
!
!---------------------------------------------------------------------
!
      return
      end subroutine inp_tv

!**s/r inp_src_surface - Read surface information from input dataset

      subroutine inp_src_surface ( F_sq,F_su,F_sv,F_gz_q,F_gz_u,F_gz_v,&
                                   F_topo,F_ip1,F_nka_gz,F_nka )
      use dynkernel_options

      implicit none

      integer, intent(in) :: F_nka
      integer, intent(out):: F_nka_gz
      integer, dimension (:    ), pointer, intent(inout) :: F_ip1
      real   , dimension (:,:,:), pointer, intent(inout) :: &
                                         F_gz_q,F_gz_u,F_gz_v
      real   , dimension (:,:), pointer, intent(inout) :: F_sq,F_su,F_sv
      real   , dimension (*  ),          intent(in   ) :: F_topo

      integer err,k,nk,kind
      integer, dimension (:    ), pointer     :: ip1_list
      real   , dimension (:,:,:), pointer     :: wrk
      real   , dimension (:    ), allocatable :: rna
!
!---------------------------------------------------------------------
!
      if (associated(F_sq)) deallocate (F_sq)
      if (associated(F_su)) deallocate (F_su)
      if (associated(F_sv)) deallocate (F_sv)
      nullify (F_sq,F_su,F_sv,ip1_list,wrk)

      if (Inp_dst_hauteur_L.and..not.Schm_autobar_L) then
         call inp_src_gz ( F_sq,F_gz_q,F_gz_u,F_gz_v,F_ip1,F_nka_gz)
         return
      endif

      err = inp_read_mt ( 'SFCPRES', ['Q','U','V'], wrk, 3,&
                           ip1_list, nk)
      if (associated(wrk)) then
         allocate ( F_sq(l_minx:l_maxx,l_miny:l_maxy),&
                    F_su(l_minx:l_maxx,l_miny:l_maxy),&
                    F_sv(l_minx:l_maxx,l_miny:l_maxy) )
         F_sq(:,:) = wrk(:,:,1)
         F_su(:,:) = wrk(:,:,2)
         F_sv(:,:) = wrk(:,:,3)
         deallocate (wrk,ip1_list) ; nullify (wrk,ip1_list)

      else
         if (Inp_kind == 2.or.Schm_autobar_L) then
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
                            nk,1,l_ni,1,l_nj )
               F_su(:,:) = rna(nk)
               F_sv(:,:) = rna(nk)
               deallocate (rna,wrk,ip1_list) ; nullify (wrk,ip1_list)
            else
               err = -1
            end if
         end if
      endif

      call gem_error ( err,'inp_src_surface','MISSING SURFACE DATA' )
!
!---------------------------------------------------------------------
!
      return
      end subroutine inp_src_surface

!**s/r inp_src_gz - Read surface information from input dataset

      subroutine inp_src_gz (F_sq,F_gz_q,F_gz_u,F_gz_v,F_ip1,F_nka)
      implicit none

      integer, intent(out)  :: F_nka
      integer, dimension (:    ), pointer, intent(inout) :: F_ip1
      real   , dimension (:,:  ), pointer, intent(inout) :: F_sq
      real   , dimension (:,:,:), pointer, intent(inout) :: &
                                         F_gz_q,F_gz_u,F_gz_v

      character(len=128) :: message_S
      integer err(4), nk, kind
      integer, dimension (:    ), pointer     :: ip1_list
      real   , dimension (:,:,:), pointer     :: wrk
      real lev, limit_near_sfc
!
!---------------------------------------------------------------------
!
      err = -1
      if (associated(F_ip1 )) deallocate (F_ip1 )
      if (associated(F_sq  )) deallocate (F_sq  )
      if (associated(F_gz_q)) deallocate (F_gz_q)
      if (associated(F_gz_u)) deallocate (F_gz_u)
      if (associated(F_gz_v)) deallocate (F_gz_v)
      nullify (F_ip1,F_sq,F_gz_q,F_gz_u,F_gz_v,wrk,ip1_list)

      err(1) = inp_read_mt ( 'SFCPRES', 'Q', wrk, 1,&
                             ip1_list, nk )

      if (associated(wrk)) then
         allocate ( F_sq(l_minx:l_maxx,l_miny:l_maxy) )
         F_sq(:,:) = wrk(:,:,1)
         deallocate (wrk,ip1_list) ; nullify (wrk,ip1_list)
         err(2) = 0
      else
         message_S= 'MISSING SURFACE P0 DATA'
      endif
      err(3)= inp_read_mt ( 'GEOPOTENTIAL', ['Q','U','V'], &
                             wrk, 3, F_ip1, F_nka )
      if (err(3) == 0) then
         call convip (F_ip1(F_nka), lev, kind, -1, message_S, .false.)
         limit_near_sfc=abs(lev-1.)
         if (Inp_src_hauteur_L) limit_near_sfc=lev
         if( limit_near_sfc > epsilon(lev)) then
            err(4) = -1
            message_S= 'Missing field: GZ at hyb or eta = 1.0'
         else
            err(4) = 0
            allocate ( F_gz_q(l_minx:l_maxx,l_miny:l_maxy,F_nka),&
                       F_gz_u(l_minx:l_maxx,l_miny:l_maxy,F_nka),&
                       F_gz_v(l_minx:l_maxx,l_miny:l_maxy,F_nka) )
            F_gz_q = wrk(:,:,        1:  F_nka)
            F_gz_u = wrk(:,:,  F_nka+1:2*F_nka)
            F_gz_v = wrk(:,:,2*F_nka+1:3*F_nka)
            deallocate (wrk)
         endif
      else
         message_S= 'MISSING GZ DATA'
      endif

      call gem_error ( minval(err),'inp_src_gz', message_S)
!
!---------------------------------------------------------------------
!
      return
      end subroutine inp_src_gz

!**s/r inp_surface - compute destination surface information

      subroutine inp_dst_surface (F_p0,F_q,F_u,F_v,F_lsq,F_lsu,F_lsv  ,&
                                  F_ip1list, F_sq0, F_me, F_tv        ,&
                                  F_topo,F_sleve_L,Minx,Maxx,Miny,Maxy,&
                                  F_nka, F_i0, F_in, F_j0, F_jn)
      use gmm_pw
      use gmm_geof
      use gmm_itf_mod
      implicit none

      logical, intent(in) :: F_sleve_L
      integer, intent(in) :: Minx, Maxx, Miny, Maxy,F_nka,&
                             F_i0,F_in,F_j0,F_jn
      integer, dimension (:) , pointer, intent(inout) :: F_ip1list
      real, dimension(*), intent(in) :: F_tv
      real, dimension(:,:), pointer, intent(in) :: F_sq0
      real, dimension(Minx:Maxx,Miny:Maxy), intent(out) :: F_p0
      real, dimension(Minx:Maxx,Miny:Maxy), intent(in) :: F_topo
      real   , dimension (:,:), pointer, intent(inout) :: F_q,F_u,F_v,&
                                                    F_lsq,F_lsu,F_lsv
      real   , dimension (:,:), pointer , intent(inout) :: F_me

      character(len=4) vname
      integer i,j,istat,typ
      real lev
      real, dimension (:,:  ), pointer :: dummy
      real, dimension (:,:,:), pointer :: srclev
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
                            .and. .not.Inp_src_hauteur_L) then
         if (lun_out > 0) then
            write(lun_out,'(" PERFORMING surface pressure adjustment")')
         end if
         allocate (srclev(l_minx:l_maxx,l_miny:l_maxy,F_nka))
         F_q(1:l_ni,1:l_nj)= F_sq0(1:l_ni,1:l_nj)
         nullify(dummy)
         call inp_3dpres ( Inp_vgd_src,F_ip1list,F_sq0,dummy,&
                           srclev,1,F_nka )
         srclev(1:l_ni,1:l_nj,F_nka)= F_sq0(1:l_ni,1:l_nj)
         call adj_ss2topo ( F_p0, F_topo, srclev, F_me, F_tv,&
                           Minx,Maxx,Miny,Maxy,F_nka,F_i0,F_in,F_j0,F_jn)
         deallocate (F_me,srclev) ; nullify (F_me)
      else
         if (lun_out > 0) then
            write(lun_out,'(" NO surface pressure adjustment")')
         end if
         F_p0(1:l_ni,1:l_nj)= F_sq0(1:l_ni,1:l_nj)
      end if

      if (Inp_dst_hauteur_L) then
         F_q(1:l_ni,1:l_nj)= F_topo(1:l_ni,1:l_nj) / grav_8
      else
         F_q(1:l_ni,1:l_nj)= F_p0(1:l_ni,1:l_nj)
      endif

      call rpn_comm_xch_halo ( F_q, l_minx,l_maxx,l_miny,l_maxy,&
           l_ni,l_nj,1,G_halox,G_haloy,G_periodx,G_periody,l_ni,0 )

      do j=1,l_nj
      do i=1,l_ni-east
         F_u(i,j)= (F_q(i,j)+F_q(i+1,j))*.5d0
      end do
      end do
      if (l_east) F_u(l_ni,1:l_nj)= F_q(l_ni,1:l_nj)

      do j=1,l_nj-north
      do i=1,l_ni
         F_v(i,j)= (F_q(i,j)+F_q(i,j+1))*.5d0
      end do
      end do
      if (l_north) F_v(1:l_ni,l_nj)= F_q(1:l_ni,l_nj)

      if ( F_sleve_L ) then
         allocate (F_lsu(l_minx:l_maxx,l_miny:l_maxy), &
                   F_lsv(l_minx:l_maxx,l_miny:l_maxy) )
         if (Inp_dst_hauteur_L) then
            istat = gmm_get (gmmk_sls_s     , F_lsq )
         else
            istat = gmm_get (gmmk_pw_p0_ls_s, F_lsq )
         endif

         call rpn_comm_xch_halo ( F_lsq, l_minx,l_maxx,l_miny,l_maxy,&
              l_ni,l_nj,1,G_halox,G_haloy,G_periodx,G_periody,l_ni,0 )
         do j=1,l_nj
         do i=1,l_ni-east
            F_lsu(i,j)= (F_lsq(i,j)+F_lsq(i+1,j))*.5d0
         end do
         end do
         if (l_east) F_lsu(l_ni,1:l_nj)= F_lsq(l_ni,1:l_nj)

         do j=1,l_nj-north
         do i=1,l_ni
            F_lsv(i,j)= (F_lsq(i,j)+F_lsq(i,j+1))*.5d0
         end do
         end do
         if (l_north) F_lsv(1:l_ni,l_nj)= F_lsq(1:l_ni,l_nj)
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
      implicit none

      character(len=*)                  , intent(in)  :: F_target_S
      integer                           , intent(out) :: F_nka
      integer, dimension(:    ), pointer, intent(out) :: F_ip1
      real   , dimension(:,:,:), pointer, intent(out) :: F_u, F_v

!     local variables
      integer, external :: RPN_COMM_shuf_ezdist, samegrid_rot
      character(len=1) typ,grd
      character(len=4) var
      character(len=12) lab
      logical, dimension (:), allocatable :: zlist_o
      integer, parameter :: nlis = 1024
      integer liste_u(nlis),liste_v(nlis),liste_sorted(nlis),lislon
      integer i,err, nz, n1,n2,n3, nrec, cnt, same_rot, ni_dest, nj_dest
      integer mpx,local_nk,irest,kstart, src_gid, dst_gid, vcode
      integer dte, det, ipas, p1, p2, p3, g1, g2, g3, g4, bit, &
              dty, swa, lng, dlf, ubc, ex1, ex2, ex3
      integer, dimension(:  ), allocatable :: zlist
      real   , dimension(:,:), allocatable :: u,v,uhr,vhr,uv,wku,wkv
      real   , dimension(:  ), pointer     :: posxu,posyu,posxv,posyv
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
         nrec= fstinl (Inp_handle,n1,n2,n3,Inp_cmcdate,' ',-1,-1,-1,' ',&
                       'UU',liste_u,lislon,nlis)
         nrec= fstinl (Inp_handle,n1,n2,n3,Inp_cmcdate,' ',-1,-1,-1,' ',&
                       'VV',liste_v,lislon,nlis)
         if (lislon == 0) goto 999

         err= fstprm (liste_u(1), DTE, DET, IPAS, n1, n2, n3,&
                  BIT, DTY, P1, P2, P3, TYP, VAR, LAB, GRD,&
                  G1,G2,G3,G4,SWA,LNG,DLF,UBC,EX1,EX2,EX3)

         src_gid = ezqkdef (n1, n2, GRD, g1, g2, g3, g4, Inp_handle)
!         same_rot= samegrid_rot ( src_gid, &
!                        Hgc_ig1ro,Hgc_ig2ro,Hgc_ig3ro,Hgc_ig4ro)

         call sort_ip1 (liste_u,liste_sorted,lislon)
         call sort_ip1 (liste_v,liste_sorted,lislon)

         allocate (F_ip1(max(1,lislon)))
         if (lislon > 1) then
            F_ip1(1:lislon) = liste_sorted(1:lislon)
         else
            F_ip1(1) = p1
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
         allocate (uhr(G_ni*G_nj,nz),vhr(G_ni*G_nj,nz),&
                   wku(ni_dest*nj_dest,nz),wkv(ni_dest*nj_dest,nz))
         allocate (u(n1*n2,max(local_nk,1)), v(n1*n2,max(local_nk,1)))

         cnt= 0
         do i= kstart, kstart+local_nk-1
            cnt= cnt+1
            err= fstluk ( u(1,cnt), liste_u(i), n1,n2,n3)
            err= fstluk ( v(1,cnt), liste_v(i), n1,n2,n3)
         end do

         if (local_nk > 0) then

            err = ezsetopt ('INTERP_DEGREE', 'CUBIC')

            if (trim(F_target_S) == 'UV') then

               posxu => geomh_lonF
               posyu => geomh_latQ
               posxv => geomh_lonQ
               posyv => geomh_latF

               write(output_unit,1001) 'Interpolating: UU, nka= ',&
                             lislon,', valid: ',Inp_datev,' on U grid'
               dst_gid = ezgdef_fmem ( ni_dest, nj_dest, 'Z', 'E', &
                       Hgc_ig1ro, Hgc_ig2ro, Hgc_ig3ro, Hgc_ig4ro, &
                                                      posxu, posyu )
               err = ezdefset ( dst_gid , src_gid )
               allocate (uv(ni_dest*nj_dest,nz))
               if (same_rot > 0) then
                  do i=1,local_nk
                     err = ezsint (wku(1,i), u(1,i))
                  end do
               else
                  do i=1,local_nk
                     err = ezuvint  ( wku(1,i),uv(1,i), u(1,i),v(1,i))
                  end do
               end if
               if (err == 2) then
                 write(output_unit,1002) 'EXTRApolating: UU, nka= ',&
                             lislon,', valid: ',Inp_datev,' on U grid'
               end if

               write(output_unit,1001) 'Interpolating: VV, nka= ',&
                             lislon,', valid: ',Inp_datev,' on V grid'

               dst_gid = ezgdef_fmem ( ni_dest, nj_dest, 'Z', 'E', &
                       Hgc_ig1ro, Hgc_ig2ro, Hgc_ig3ro, Hgc_ig4ro, &
                                                      posxv, posyv )
               err = ezdefset ( dst_gid , src_gid )

               if (same_rot > 0) then
                  do i=1,local_nk
                     err = ezsint (wkv(1,i), v(1,i))
                  end do
               else
                  do i=1,local_nk
                     err = ezuvint  ( uv(1,i),wkv(1,i), u(1,i),v(1,i))
                  end do
               end if
               if (err == 2) then
                 write(output_unit,1002) 'EXTRApolating: VV, nka= ',&
                             lislon,', valid: ',Inp_datev,' on V grid'
               end if
               deallocate (uv)

            else

               posxu => geomh_lonQ
               posyu => geomh_latQ

               write(output_unit,1001) 'Interpolating: UV, nka= ',&
                              lislon,', valid: ',Inp_datev,' on Q grid'
               dst_gid = ezgdef_fmem ( ni_dest, nj_dest, 'Z', 'E', &
                       Hgc_ig1ro, Hgc_ig2ro, Hgc_ig3ro, Hgc_ig4ro, &
                                                      posxu, posyu )
               err = ezdefset ( dst_gid , src_gid )
               if (same_rot > 0) then
                  do i=1,local_nk
                     err = ezsint (wku(1,i), u(1,i))
                     err = ezsint (wkv(1,i), v(1,i))
                  end do
               else
                  do i=1,local_nk
                     err = ezuvint  ( wku(1,i),wkv(1,i), u(1,i),v(1,i))
                  end do
               end if
               if (err == 2) then
                 write(output_unit,1002) 'EXTRApolating: UV, nka= ',&
                             lislon,', valid: ',Inp_datev,' on Q grid'
               end if

            end if

            do i=1,local_nk
               call reshape_wk ( wku(1,i),uhr(1,i),ni_dest,nj_dest, &
                                 G_ni,G_nj,G_halox,G_haloy )
               call reshape_wk ( wkv(1,i),vhr(1,i),ni_dest,nj_dest, &
                                 G_ni,G_nj,G_halox,G_haloy )
            end do

         end if

         deallocate (u,v,wku,wkv)

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

         allocate ( F_u(l_minx:l_maxx,l_miny:l_maxy,lislon), &
                    F_v(l_minx:l_maxx,l_miny:l_maxy,lislon), &
                    zlist_o(lislon) )

         zlist_o= .FALSE.

         err = RPN_COMM_shuf_ezdist ( Inp_comm_setno, Inp_comm_id, &
                              uhr, nz, F_u, lislon, zlist, zlist_o )
         zlist_o= .FALSE.

         err = RPN_COMM_shuf_ezdist ( Inp_comm_setno, Inp_comm_id, &
                              vhr, nz, F_v, lislon, zlist, zlist_o )

         deallocate (uhr,vhr,zlist,zlist_o)

         F_u(1:l_ni,1:l_nj,:) = F_u(1:l_ni,1:l_nj,:) * knams_8
         F_v(1:l_ni,1:l_nj,:) = F_v(1:l_ni,1:l_nj,:) * knams_8

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
      integer istat
      logical inlog
      real, dimension (:,:  ), pointer :: pres,presl
      real, dimension (:,:,:), pointer :: ptr3d
!
!---------------------------------------------------------------------
!
      inlog= .false.
      if ( present(F_inlog_S) ) then
         if (F_inlog_S == 'in_log') inlog= .true.
      end if

      pres  => F_sfc (1:l_ni,1:l_nj      )
      ptr3d => F_dest(1:l_ni,1:l_nj,k0:kn)

      if ( associated (F_sfcL) ) then
         presl=> F_sfcL (1:l_ni,1:l_nj)
         istat= vgd_levels (F_vgd, F_ip1(k0:kn), ptr3d, sfc_field=pres, &
                            sfc_field_ls=presl, in_log=inlog)
      else
         istat= vgd_levels (F_vgd, F_ip1(k0:kn), ptr3d, pres, &
                            in_log=inlog)
      end if
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

      type(vgrid_descriptor) , intent(in ) :: F_vgd
      integer                , intent(out) :: F_nk
      integer                , intent(in)  :: Minx,Maxx,Miny,Maxy
      integer, dimension(:)    , pointer, intent(in ) :: F_ip1,F_gz_ip1
      real   , dimension(:,:  ), pointer, intent(in ) :: F_sfc,F_sfcL
      real   , dimension(:,:,:), pointer, intent(in ) :: F_gz
      real   , dimension(:,:,:), pointer, intent(inout) :: F_dest

!     local variables
      integer istat
!
!---------------------------------------------------------------------
!
      F_nk= ubound(F_ip1,1)
      if (Inp_dst_hauteur_L.and..not.Schm_autobar_L) then
         istat= inp_match_heights (F_dest, F_gz, F_ip1, F_gz_ip1, &
                       Minx,Maxx,Miny,Maxy,F_nk, ubound(F_gz_ip1,1))
      else
         call inp_3dpres ( F_vgd,  F_ip1,  F_sfc,  F_sfcL, F_dest, &
                           1, F_nk, F_inlog_S='in_log')
      endif
!
!---------------------------------------------------------------------
!
      return
      end subroutine inp_src_levels

      subroutine inp_dst_levels ( F_dest, F_vgd, F_ip1, F_sfc, F_sfcL )
      implicit none

      type(vgrid_descriptor) , intent(in) :: F_vgd
      integer, dimension(:)    , pointer, intent(in ) :: F_ip1
      real   , dimension(:,:  ), pointer, intent(in ) :: F_sfc,F_sfcL
      real   , dimension(:,:,:), pointer, intent(inout) :: F_dest
!
!---------------------------------------------------------------------
!
      if (Inp_dst_hauteur_L) then
         call inp_3dhgts ( F_vgd, F_ip1, F_sfc, F_sfcL, F_dest, 1, G_nk)
      else
         call inp_3dpres ( F_vgd, F_ip1, F_sfc, F_sfcL, F_dest, &
                           1, G_nk, F_inlog_S='in_log')
      endif
!
!---------------------------------------------------------------------
!
      return
      end subroutine inp_dst_levels
!
!---------------------------------------------------------------------
!
!**s/r inp_match_heights - To match heights
!
      integer function inp_match_heights ( F_ho, F_hi, F_ip1_list_o, &
            F_ip1_list_i, Minx,Maxx,Miny,Maxy, nko, nki) result(status)

      integer, intent(inout) :: nko
      integer, intent(in)    :: Minx,Maxx,Miny,Maxy,nki
      integer, dimension (:), pointer, intent(in) :: &
                             F_ip1_list_i,F_ip1_list_o
      real, dimension (:,:,:), pointer, intent(out) :: F_ho
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
                 call gem_error (-1,'inp_match_heights',&
                        'more than one diagnostic level, review code')
            index_diag_AGL = ko
            F_ho(:,:,ko) = F_hi(:,:,nki) + lev
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
            call gem_error (-1,'inp_match_heights',message_S)
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
      end function inp_match_heights

      subroutine inp_3dhgts ( F_vgd, F_ip1, F_sfc, F_sfcL, F_dest,k0,kn )

      type(vgrid_descriptor) , intent(IN) :: F_vgd
      integer                , intent(IN) :: k0,kn
      integer, dimension(:)    , pointer, intent(in ) :: F_ip1
      real   , dimension(:,:  ), pointer, intent(in ) :: F_sfc,F_sfcL
      real   , dimension(:,:,:), pointer, intent(out) :: F_dest

!     local variables
      integer :: istat
      real, dimension (:,:  ), pointer :: hgt,hgtls
      real, dimension (:,:,:), pointer :: ptr3d
!
!---------------------------------------------------------------------
!
      hgt  => F_sfc (1:l_ni,1:l_nj      )
      ptr3d => F_dest(1:l_ni,1:l_nj,k0:kn)

      hgtls=> F_sfcL (1:l_ni,1:l_nj)
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
                            F_ssqr,F_ssur,F_ssvr, F_ssq0,F_ssu0,F_ssv0,&
                            F_ssq0LS,F_ssu0LS,F_ssv0LS                ,&
                            F_gz_q, F_gz_u, F_gz_v, F_GZ_ip1          ,&
                            Minx,Maxx,Miny,Maxy, F_nk )
      implicit none

      logical                , intent(in)  :: F_stag_L
      integer                , intent(in)  :: Minx,Maxx,Miny,Maxy, F_nk
      integer, dimension(:)  , pointer, intent(in)  :: F_GZ_ip1
      real, dimension (:,:), pointer, intent(in) :: &
                             F_ssqr,F_ssur,F_ssvr, F_ssq0,F_ssu0,F_ssv0,&
                             F_ssq0LS,F_ssu0LS,F_ssv0LS
      real, dimension (:,:,:), pointer, intent(in) :: F_gz_q, F_gz_u, &
                                                      F_gz_v
      real, dimension(Minx:Maxx,Miny:Maxy,F_nk), intent(out) :: F_u, F_v

!     local variables
      integer nka
      integer, dimension (:    ), pointer :: ip1_list, ip1_target
      real   , dimension (:,:  ), pointer :: dummy
      real   , dimension (:,:,:), pointer :: srclev,dstlev
      real   , dimension (:,:,:), pointer :: ur,vr
!
!---------------------------------------------------------------------
!
      nullify (ip1_list, ur, vr, dummy)
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
                            dummy, F_gz_u, F_GZ_ip1,Minx,Maxx,Miny,Maxy)

         call inp_dst_levels (dstlev, Ver_vgdobj, ip1_target,&
                              F_ssu0, F_ssu0Ls)

         call vertint2 ( F_u,dstlev,G_nk, ur,srclev,nka,&
                         l_minx,l_maxx,l_miny,l_maxy, 1,l_ni,1,l_nj,&
                         levtype=Inp_levtype_S )

         call inp_src_levels (srclev, nka, ip1_list, Inp_vgd_src, F_ssvr,&
                            dummy, F_gz_v, F_GZ_ip1,Minx,Maxx,Miny,Maxy)

         call inp_dst_levels (dstlev, Ver_vgdobj, ip1_target,&
                              F_ssv0, F_ssv0Ls)

         call vertint2 ( F_v,dstlev,G_nk, vr,srclev,nka,&
                         l_minx,l_maxx,l_miny,l_maxy, 1,l_ni,1,l_nj,&
                         levtype=Inp_levtype_S )

      else

         call inp_src_levels (srclev, nka, ip1_list, Inp_vgd_src, F_ssqr,&
                            dummy, F_gz_q, F_GZ_ip1,Minx,Maxx,Miny,Maxy)

         call inp_dst_levels (dstlev, Ver_vgdobj, ip1_target,&
                              F_ssq0, F_ssq0Ls)

         call vertint2 ( F_u,dstlev,G_nk, ur,srclev,nka,&
                         l_minx,l_maxx,l_miny,l_maxy, 1,l_ni,1,l_nj,&
                         levtype=Inp_levtype_S )
         call vertint2 ( F_v,dstlev,G_nk, vr,srclev,nka,&
                         l_minx,l_maxx,l_miny,l_maxy, 1,l_ni,1,l_nj,&
                         levtype=Inp_levtype_S)

      end if

      deallocate (ip1_list,ip1_target,ur,vr,srclev,dstlev)
      nullify    (ip1_list,ip1_target,ur,vr,srclev,dstlev)
!
!---------------------------------------------------------------------
!
      return
      end subroutine inp_hwnd

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
