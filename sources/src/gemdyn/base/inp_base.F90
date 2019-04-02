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
   use inp_mod
   use tdpack
   implicit none
#include <arch_specific.hf>
#include <rmnlib_basics.hf>

   public

contains

!**s/r inp_get - Read variable F_var_S valid at F_datev, perform hor.
!                interpolation to F_hgrid_S hor. grid and perform
!                vertical interpolation to F_vgrid_S vertical grid

      integer function inp_get ( F_var_S, F_hgrid_S, F_ver_ip1    ,&
                                 F_vgd_src, F_vgd_dst             ,&
                                 F_sfc_src, F_sfc_dst, F_sfcLS_dst,&
                                 F_dest , Minx,Maxx,Miny,Maxy     ,&
                                 F_nk   , F_inttype_S )

      character(len=*)          , intent(in) :: F_var_S,F_hgrid_S
      character(len=*), optional, intent(in) :: F_inttype_S
      integer                   , intent(in) :: Minx,Maxx,Miny,Maxy, F_nk
      integer, dimension(:)  , pointer, intent(in) :: F_ver_ip1
      type(vgrid_descriptor) , intent(in) :: F_vgd_src, F_vgd_dst
      real, dimension (:,:), pointer, intent(in) :: &
                                      F_sfc_src,F_sfc_dst,F_sfcLS_dst
      real, dimension(Minx:Maxx,Miny:Maxy,F_nk), intent(out):: F_dest

!     local variables
      character(len=12) :: inttype
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

      inp_get= inp_read ( F_var_S, F_hgrid_S, wrkr, ip1_list, nka )

      if (inp_get < 0) then
         if (associated(ip1_list)) deallocate (ip1_list)
         if (associated(wrkr    )) deallocate (wrkr    )
         return
      end if

      allocate ( srclev(l_minx:l_maxx,l_miny:l_maxy,nka) ,&
                 dstlev(l_minx:l_maxx,l_miny:l_maxy,G_nk) )

      call inp_3dpres ( F_vgd_src,  ip1_list, F_sfc_src, dummy      ,&
                        srclev, 1,  nka, F_inlog_S='in_log')
      call inp_3dpres ( F_vgd_dst, F_ver_ip1, F_sfc_dst, F_sfcLS_dst,&
                        dstlev, 1, G_nk, F_inlog_S='in_log')

      inttype= 'cubic'
      if (present(F_inttype_S)) inttype= F_inttype_S
      call vertint2 ( F_dest,dstlev,G_nk, wrkr,srclev,nka       ,&
                      l_minx,l_maxx,l_miny,l_maxy, 1,l_ni,1,l_nj,&
                      varname=F_var_S, inttype= inttype)

      deallocate (ip1_list,wrkr,srclev,dstlev)
!
!---------------------------------------------------------------------
!
      return
      end function inp_get

!**s/r inp_read - Parallel read of variable F_var_S and horizontal
!                 interpolation to F_hgrid_S Arakawa grid

      integer function inp_read ( F_var_S, F_hgrid_S, F_dest, &
                                  F_ip1, F_nka, F_hint_S )

      use glb_pil
      use hgc
      implicit none

      character(len=*)                     ,intent(in)  :: F_var_S,F_hgrid_S
      character(len=*),            optional,intent(in)  :: F_hint_S
      integer                           ,intent(out) :: F_nka
      integer, dimension(:    ), pointer,intent(out) :: F_ip1
      real   , dimension(:,:,:), pointer,intent(inout) :: F_dest


!     local variables
      integer, external :: RPN_COMM_shuf_ezdist, &
                           samegrid_gid, samegrid_rot, inp_is_real_wind
      character(len=1) typ,grd
      character(len=4) nomvar,var
      character(len=12) lab,interp_S
      logical, dimension (:), allocatable :: zlist_o
      integer, parameter :: nlis = 1024
      integer i,err, nz, n1,n2,n3, nrec, liste(nlis),lislon,cnt
      ! Remove the following line by 2021
      integer ni_dest,nj_dest,ut1_is_urt1
      integer subid,nicore,njcore,datev
      integer mpx,local_nk,irest,kstart, src_gid, vcode, ip1
      integer dte, det, ipas, p1, p2, p3, g1, g2, g3, g4, bit, &
              dty, swa, lng, dlf, ubc, ex1, ex2, ex3
      integer, dimension(:  ), allocatable :: zlist
      real   , dimension(:,:), allocatable :: wk1,wk2,wk3
      real   , dimension(:  ), pointer     :: posx,posy
      real*8 add, mult
      common /bcast_i / lislon,nz
!
!---------------------------------------------------------------------
!
      inp_read= -1
      F_nka= -1 ; local_nk= 0
      add= 0.d0 ; mult= 1.d0
      ! Remove the following line by 2021
      ut1_is_urt1 = -1
      if (associated(F_ip1 )) deallocate (F_ip1 )
      if (associated(F_dest)) deallocate (F_dest)
      nullify (F_ip1,F_dest)

      nomvar = F_var_S ; ip1= -1
      select case (F_var_S)
         case ('OROGRAPHY')
            if (Inp_kind == 2  ) nomvar= '@NUL'
            if (Inp_kind == 1  ) nomvar= 'GZ'
            if (Inp_kind == 5  ) nomvar= 'GZ'
            if (Inp_kind == 105) nomvar= 'FIS0'
            if ( nomvar == 'GZ' ) ip1= 93423264
            if (( nomvar == 'GZ' ) .and. (Inp_kind == 1  )) ip1= 12000
            if ( nomvar == 'GZ' ) mult= 10.d0 * grav_8
         case ('SFCPRES')
            if (Inp_kind == 2  ) nomvar= '@NUL'
            if (Inp_kind == 1  ) nomvar= 'P0'
            if (Inp_kind == 5  ) nomvar= 'P0'
            if (Inp_kind == 105) nomvar= 'ST1'
            if ( nomvar == 'P0' ) mult= 100.d0
         case ('TEMPERATURE')
            if (Inp_kind == 2  ) nomvar= 'TT'
            if (Inp_kind == 1  ) nomvar= 'TT'
            if (Inp_kind == 5  ) nomvar= 'TT'
            if (Inp_kind == 105) nomvar= 'TT1'
            if ( nomvar == 'TT' ) add= tcdk_8
         case ('GEOPOTENTIAL')
            if (Inp_kind == 2  ) nomvar= 'GZ'
            if (Inp_kind == 1  ) nomvar= '@NUL'
            if (Inp_kind == 5  ) nomvar= 'GZ'
            if (Inp_kind == 105) nomvar= '@NUL'
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

         allocate (F_ip1(lislon))
         if (lislon > 1) then
            call sort_ip1 (liste,F_ip1,lislon)
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
         allocate (wk3(ni_dest*nj_dest,nz+1))
         allocate (wk2(G_ni*G_nj,nz+1))
         allocate (wk1(n1*n2,max(local_nk,1)))

         cnt= 0
         do i= kstart, kstart+local_nk-1
            cnt= cnt+1
            err= fstluk (wk1(1,cnt), liste(i), n1,n2,n3)
            ! Remove the following line by 2021
            ut1_is_urt1 = inp_is_real_wind (wk1(1,cnt),n1*n2,nomvar)
         end do

         if (local_nk > 0) then

            if (F_hgrid_S == 'Q') then
               posx => geomh_lonQ
               posy => geomh_latQ
            end if
            if (F_hgrid_S == 'U') then
               posx => geomh_lonF
               posy => geomh_latQ
            end if
            if (F_hgrid_S == 'V') then
               posx => geomh_lonQ
               posy => geomh_latF
            end if
            if (F_hgrid_S == 'F') then
               posx => geomh_lonF
               posy => geomh_latF
            end if

            dstf_gid = ezgdef_fmem (ni_dest, nj_dest, 'Z', 'E', &
                    Hgc_ig1ro, Hgc_ig2ro, Hgc_ig3ro, Hgc_ig4ro, &
                    posx, posy)

            interp_S= 'CUBIC'
            if (present(F_hint_S)) interp_S= F_hint_S

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
            write(6,1001) 'Interpolating: ',trim(nomvar),', nka= ',&
               lislon,', valid: ',Inp_datev,' on ',F_hgrid_S, ' grid'
         end if

         do i=1,local_nk
            err = ezsint(wk3(1,i), wk1(1,i))
            call reshape_wk ( wk3(1,i),wk2(1,i),ni_dest,nj_dest,&
                              G_ni,G_nj,G_halox,G_haloy )
         end do
         if (err == 2) then
            write(6,1002) 'EXTRApolating: ',trim(nomvar),', nka= ',&
               lislon,', valid: ',Inp_datev,' on ',F_hgrid_S, ' grid'
         end if
         err = ezsetopt ( 'USE_1SUBGRID', 'NO' )
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

         inp_read= 0
         if (F_nka >= 1) then
            if (Inp_iome < 0) allocate ( F_ip1(F_nka) )
            call rpn_comm_bcast ( F_ip1, F_nka, "MPI_INTEGER", &
                                  Inp_iobcast, "grid", err )
         end if

         allocate (zlist(nz)) ; zlist= -1
         do i=1, local_nk
            zlist(i)= i + kstart - 1
         end do

         allocate ( F_dest(l_minx:l_maxx,l_miny:l_maxy,lislon), &
                    zlist_o(lislon) )

         zlist_o= .FALSE.

         err = RPN_COMM_shuf_ezdist ( Inp_comm_setno, Inp_comm_id, &
                           wk2, nz, F_dest, lislon, zlist, zlist_o )

         deallocate (wk2,zlist,zlist_o)

         F_dest(1:l_ni,1:l_nj,:) = F_dest(1:l_ni,1:l_nj,:) * mult + add
         if (nomvar == 'ST1') &
         F_dest(1:l_ni,1:l_nj,:)= Inp_pref_a_8 * &
                                  exp(F_dest(1:l_ni,1:l_nj,:))

      else

         inp_read= -1
         if (Inp_iome >= 0) write(6,'(7a)') ' FIELD: ',trim(F_var_S),&
                     ':',trim(nomvar),' valid: ',Inp_datev, 'NOT FOUND'

      end if

 1001 format (3a,i3,5a)
 1002 format (3a,i3,5a)
!
!---------------------------------------------------------------------
!
      return
      end function inp_read

!**s/r inp_read_uv - Read horizontal winds UU,VV F valid at F_datev
!                    and perform vectorial horizontal interpolation
!                    on proper Arakawa grid u and v respectively

      subroutine inp_read_uv ( F_u, F_v, F_target_S, F_ip1, F_nka )

      use hgc
      implicit none

      character(len=*)                     , intent(in)  :: F_target_S
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
      integer liste_u(nlis),liste_v(nlis),lislon
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

         allocate (F_ip1(lislon))
         if (lislon > 1) then
            call sort_ip1 (liste_u,F_ip1,lislon)
            call sort_ip1 (liste_v,F_ip1,lislon)
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

               write(6,1001) 'Interpolating: UU, nka= ',&
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
                 write(6,1002) 'EXTRApolating: UU, nka= ',&
                             lislon,', valid: ',Inp_datev,' on U grid'
               end if

               write(6,1001) 'Interpolating: VV, nka= ',&
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
                 write(6,1002) 'EXTRApolating: VV, nka= ',&
                             lislon,', valid: ',Inp_datev,' on V grid'
               end if
               deallocate (uv)

            else

               posxu => geomh_lonQ
               posyu => geomh_latQ

               write(6,1001) 'Interpolating: UV, nka= ',&
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
                 write(6,1002) 'EXTRApolating: UV, nka= ',&
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

         if (Inp_iome >= 0) write(6,'(3a)') &
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
      call handle_error_l(istat==VGD_OK,'inp_3dpres','Error computing pressure')
!
!---------------------------------------------------------------------
!
      return
      end subroutine inp_3dpres

!**s/r inp_hwnd - Read horizontal winds UU,VV F valid at F_datev
!                 and perform vectorial horizontal interpolation
!                 and  vertical interpolation to momentum levels

      subroutine inp_hwnd ( F_u,F_v, F_vgd_src,F_vgd_dst, F_stag_L    ,&
                             F_ssqr,F_ssur,F_ssvr, F_ssq0,F_ssu0,F_ssv0,&
                             F_ssq0LS,F_ssu0LS,F_ssv0LS                ,&
                             Minx,Maxx,Miny,Maxy, F_nk )
      use ver
      implicit none

      logical                , intent(in)  :: F_stag_L
      integer                , intent(in)  :: Minx,Maxx,Miny,Maxy, F_nk
      type(vgrid_descriptor) , intent(in)  :: F_vgd_src, F_vgd_dst
      real, dimension (:,:), pointer, intent(in) :: &
                             F_ssqr,F_ssur,F_ssvr, F_ssq0,F_ssu0,F_ssv0,&
                             F_ssq0LS,F_ssu0LS,F_ssv0LS
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

         call inp_3dpres ( F_vgd_src,  ip1_list, F_ssur, dummy    ,&
                           srclev, 1,  nka, F_inlog_S='in_log')
         call inp_3dpres ( F_vgd_dst, ip1_target, F_ssu0, F_ssu0Ls,&
                           dstlev, 1, G_nk, F_inlog_S='in_log')

         call vertint2 ( F_u,dstlev,G_nk, ur,srclev,nka,&
                         l_minx,l_maxx,l_miny,l_maxy, 1,l_ni,1,l_nj )

         call inp_3dpres ( F_vgd_src,  ip1_list, F_ssvr, dummy    ,&
                           srclev, 1,  nka, F_inlog_S='in_log')
         call inp_3dpres ( F_vgd_dst, ip1_target, F_ssv0, F_ssv0Ls,&
                           dstlev, 1, G_nk, F_inlog_S='in_log')

         call vertint2 ( F_v,dstlev,G_nk, vr,srclev,nka,&
                         l_minx,l_maxx,l_miny,l_maxy, 1,l_ni,1,l_nj )

      else

         call inp_3dpres ( F_vgd_src,  ip1_list, F_ssqr, dummy    ,&
                           srclev, 1,  nka, F_inlog_S='in_log')
         call inp_3dpres ( F_vgd_dst, ip1_target, F_ssq0, F_ssq0Ls,&
                           dstlev, 1, G_nk, F_inlog_S='in_log')

         call vertint2 ( F_u,dstlev,G_nk, ur,srclev,nka,&
                         l_minx,l_maxx,l_miny,l_maxy, 1,l_ni,1,l_nj )
         call vertint2 ( F_v,dstlev,G_nk, vr,srclev,nka,&
                         l_minx,l_maxx,l_miny,l_maxy, 1,l_ni,1,l_nj )

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
