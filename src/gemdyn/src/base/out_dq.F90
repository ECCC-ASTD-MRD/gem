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

!**s/r out_dq - Compute and output divergence and vorticity

      subroutine out_dq (levset,set)
      use vertical_interpolation, only: vertint2
      use gmm_vt1
      use gmm_pw
      use gem_options
      use out_options
      use glb_ld
      use out3
      use levels
      use svro_mod
      use out_meta
      use outd
      use ver
      use rmn_gmm
      use outgrid
      use set_level_mod, only: set_level_usr_val,set_level_ERROR
      implicit none

      integer levset,set

      logical write_diag_lev
      integer i,istat,kind,nko,pndd,pnqq,pnqr,pnxx,gridset
      integer, dimension(:), pointer :: indo
      real, dimension(:    ), pointer :: rf
      real, dimension(:,:,:), pointer :: cible, usr_src
      real, dimension(:,:,:), allocatable:: uu_pres,vv_pres, div,vor,qr
!
!----------------------------------------------------------------------
!
      nullify(indo,rf,cible,usr_src)
      pndd=0 ; pnqq=0 ; pnqr=0 ; write_diag_lev = .false.

      do i=1,Outd_var_max(set)
        if (Outd_var_S(i,set) == 'DD') pndd=i
        if (Outd_var_S(i,set) == 'QQ') pnqq=i
        if (Outd_var_S(i,set) == 'QR') pnqr=i
      enddo

      if (pndd+pnqq+pnqr == 0) return

      if (Level_typ_S(levset) == 'M') then

         Out_stag_S= 'MM '//OutGrid_hgrid_usr(Outd_grid(set))%usr_grid_index_S
         kind= Level_kind_ip1
         allocate (indo( min(Level_max(levset),Level_momentum) ))
         call out_slev ( Level(1,levset), Level_max(levset), &
                          Level_momentum , indo,nko,write_diag_lev )

         call rpn_comm_xch_halo (ut1,l_minx,l_maxx,l_miny,l_maxy,&
            l_niu,l_nj,G_nk,G_halox,G_haloy,G_periodx,G_periody,l_ni,0)
         call rpn_comm_xch_halo (vt1,l_minx,l_maxx,l_miny,l_maxy,&
            l_ni,l_njv,G_nk,G_halox,G_haloy,G_periodx,G_periody,l_ni,0)

         if (pndd > 0) then
            allocate ( div(l_minx:l_maxx,l_miny:l_maxy,G_nk) )
            call cal_div ( div, ut1, vt1 , Outd_filtpass(pndd,set),&
                           Outd_filtcoef(pndd,set)                ,&
                           l_minx,l_maxx,l_miny,l_maxy, G_nk )
            if ( .not. OUTs_server_L) then
            gridset = Outd_grid(set)
            call out_href ( 'Mass_point', &
                    OutGrid_x0 (gridset), OutGrid_x1 (gridset), 1, &
                    OutGrid_y0 (gridset), OutGrid_y1 (gridset), 1 )
            endif
            call out_fstecr( div, l_minx,l_maxx,l_miny,l_maxy        ,&
                              Ver_hyb%m,'DD  ',Outd_convmult(pndd,set),&
                              Outd_convadd(pndd,set),kind,-1          ,&
                              G_nk, indo,nko, Outd_nbit(pndd,set),.false.)
            deallocate (div)
         endif

         if ((pnqq > 0).or.(pnqr > 0)) then

            Out_stag_S= 'FM '//OutGrid_hgrid_usr(Outd_grid(set))%usr_grid_index_S
            allocate ( vor(l_minx:l_maxx,l_miny:l_maxy,G_nk),&
                        qr(l_minx:l_maxx,l_miny:l_maxy,G_nk) )
            if(pnqq > 0)then
               pnxx=pnqq
            else
               pnxx=pnqr
            endif
            call cal_vor ( qr, vor, ut1, vt1 , Outd_filtpass(pnxx,set),&
                           Outd_filtcoef(pnxx,set),(pnqq > 0)        ,&
                           l_minx,l_maxx,l_miny,l_maxy, G_nk )
            if ( .not. OUTs_server_L) then
            gridset = Outd_grid(set)
            call out_href ( 'F_point', &
                    OutGrid_x0 (gridset), OutGrid_x1 (gridset), 1, &
                    OutGrid_y0 (gridset), OutGrid_y1 (gridset), 1 )
            endif
            if (pnqq > 0) &
            call out_fstecr( vor, l_minx,l_maxx,l_miny,l_maxy        ,&
                              Ver_hyb%m,'QQ  ',Outd_convmult(pnqq,set),&
                              Outd_convadd(pnqq,set),kind,-1          ,&
                              G_nk, indo,nko, Outd_nbit(pnqq,set),.false.)
            if (pnqr > 0) &
            call out_fstecr( qr, l_minx,l_maxx,l_miny,l_maxy         ,&
                              Ver_hyb%m,'QR  ',Outd_convmult(pnqr,set),&
                              Outd_convadd(pnqr,set),kind,-1          ,&
                              G_nk, indo,nko, Outd_nbit(pnqr,set),.false.)
            deallocate (vor, qr)

          endif

          deallocate(indo)
          
      else ! Output on pressure, heights AGL or user levels

         ! Note: The pointers cible, usr_src, indo, and rf
         ! will point to memory allocated within set_level_usr_val.
         ! This "internal" memory allocation will persist 
         ! throughout the model integration and will be expanded if necessary. 
         ! Therefore, cible, usr_src, indo, and rf must not
         ! be deallocated, but they can be nullified if needed.

         Out_stag_S='M???'
         if( set_level_usr_val(cible, indo, rf, kind, nko, usr_src, &
              levset,Level,Level_max, 'MOMENTUM',&
              l_minx,l_maxx,l_miny,l_maxy, G_nk, &
              Out_stag_S, Level_typ_S(levset), Outd_grid(set)) &
              == set_level_ERROR )then
            print*,'TODO in out_qd handle error gracefully 1'
            return
         end if

         allocate (uu_pres(l_minx:l_maxx,l_miny:l_maxy,nko),&
                  vv_pres(l_minx:l_maxx,l_miny:l_maxy,nko) )

         call vertint2 ( uu_pres,cible,nko, ut1,usr_src,G_nk      ,&
                         l_minx,l_maxx,l_miny,l_maxy, 1,l_niu,1,l_nj,&
                         inttype=Out3_vinterp_type_S )
         call vertint2 ( vv_pres,cible,nko, vt1,usr_src,G_nk      ,&
                         l_minx,l_maxx,l_miny,l_maxy, 1,l_ni,1,l_njv,&
                         inttype=Out3_vinterp_type_S )

         call rpn_comm_xch_halo (uu_pres,l_minx,l_maxx,l_miny,l_maxy,&
            l_niu,l_nj,nko,G_halox,G_haloy,G_periodx,G_periody,l_ni,0)
         call rpn_comm_xch_halo (vv_pres,l_minx,l_maxx,l_miny,l_maxy,&
            l_ni,l_njv,nko,G_halox,G_haloy,G_periodx,G_periody,l_ni,0)

         if (pndd > 0) then
            allocate ( div(l_minx:l_maxx,l_miny:l_maxy,nko) )
            call cal_div ( div, uu_pres, vv_pres  ,&
                           Outd_filtpass(pndd,set),&
                           Outd_filtcoef(pndd,set),&
                           l_minx,l_maxx,l_miny,l_maxy, nko )
            if ( .not. OUTs_server_L) then
            gridset = Outd_grid(set)
            call out_href ( 'Mass_point', &
                    OutGrid_x0 (gridset), OutGrid_x1 (gridset), 1, &
                    OutGrid_y0 (gridset), OutGrid_y1 (gridset), 1 )
            endif
            call out_fstecr( div, l_minx,l_maxx,l_miny,l_maxy, &
                              rf,'DD  ',Outd_convmult(pndd,set),&
                              Outd_convadd(pndd,set), kind,-1  ,&
                              nko, indo, nko, Outd_nbit(pndd,set),.false.)
            deallocate (div)
         endif

         if ((pnqq > 0).or.(pnqr > 0)) then

            allocate ( vor(l_minx:l_maxx,l_miny:l_maxy,nko),&
                        qr(l_minx:l_maxx,l_miny:l_maxy,nko) )
            if(pnqq > 0)then
               pnxx=pnqq
            else
               pnxx=pnqr
            endif
            call cal_vor ( qr, vor, uu_pres, vv_pres          ,&
                           Outd_filtpass(pnxx,set)            ,&
                           Outd_filtcoef(pnxx,set),(pnqq > 0),&
                           l_minx,l_maxx,l_miny,l_maxy, nko )
            Out_stag_S(1:1)= 'F'
            if ( .not. OUTs_server_L) then
            gridset = Outd_grid(set)
            call out_href ( 'F_point', &
                    OutGrid_x0 (gridset), OutGrid_x1 (gridset), 1, &
                    OutGrid_y0 (gridset), OutGrid_y1 (gridset), 1 )
            endif
            if (pnqq > 0) &
            call out_fstecr( vor, l_minx,l_maxx,l_miny,l_maxy, &
                              rf,'QQ  ',Outd_convmult(pnqq,set),&
                              Outd_convadd(pnqq,set), kind,-1  ,&
                              nko, indo, nko, Outd_nbit(pnqq,set),.false.)
            if (pnqr > 0) &
            call out_fstecr ( qr, l_minx,l_maxx,l_miny,l_maxy, &
                              rf,'QR  ',Outd_convmult(pnqr,set),&
                              Outd_convadd(pnqr,set), kind,-1  ,&
                              nko, indo, nko, Outd_nbit(pnqr,set),.false.)
            deallocate (vor, qr)
          endif
          
          deallocate (uu_pres,vv_pres)

      endif
!
!----------------------------------------------------------------------
!
      return
      end
