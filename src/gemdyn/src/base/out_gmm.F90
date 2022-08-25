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

!**s/r out_gmm - output GMM fields

      subroutine out_gmm (levset, set)
      use dyn_fisl_options
      use glb_ld
      use out3
      use levels
      use outd
      use ver
      use gmm_table
      use outgrid
      implicit none

      integer, intent(in) :: levset,set

      type(gmm_metadata) :: tmp_meta
      character(len=4) nomvar
      logical write_diag_lev
      integer nko,i,ii,gridset,istat,NK,WNK,ind0(1),indo(G_nk+1)
      real, pointer, dimension(:,:,:) :: tr3
      real, pointer, dimension(:,:  ) :: tr2
      real, pointer, dimension(:    ) :: level_type
      real hyb0(1)
!
!     ---------------------------------------------------------------
!
      if ( Level_typ_S(levset) == 'P') return

      hyb0(1)=0.0
      ind0(1)=1

!     Setup the indexing for output
      call out_slev ( Level(1,levset), Level_max(levset),G_nk,indo,nko,write_diag_lev)
      indo(nko+1)= G_nk+1

      do ii= 1, Outd_var_max(set)
      do  i= 1, gmm_cnt

         if ( trim(Outd_varnm_S(ii,set)) == trim(GMM_tbl%vname(i)) ) then
            gridset = Outd_grid(set)
            level_type => Ver_hyb%t
            if (GMM_tbl%cn(i)(1:1) == 'M') level_type => Ver_hyb%m

            select case (GMM_tbl%ara(i))
            case('UU')
               call out_href ( 'U_point', &
                    OutGrid_x0 (gridset), OutGrid_x1 (gridset), 1, &
                    OutGrid_y0 (gridset), OutGrid_y1 (gridset), 1 )
            case('VV')
               call out_href ( 'V_point', &
                    OutGrid_x0 (gridset), OutGrid_x1 (gridset), 1, &
                    OutGrid_y0 (gridset), OutGrid_y1 (gridset), 1 )
            case('FF')
               call out_href ( 'F_point', &
                    OutGrid_x0 (gridset), OutGrid_x1 (gridset), 1, &
                    OutGrid_y0 (gridset), OutGrid_y1 (gridset), 1 )
            case default
               call out_href ( 'Mass_point', &
                    OutGrid_x0 (gridset), OutGrid_x1 (gridset), 1, &
                    OutGrid_y0 (gridset), OutGrid_y1 (gridset), 1 )
            end select

            nullify(tr2,tr3)
            istat = gmm_getmeta(GMM_tbl%vname(i),tmp_meta)
            nomvar= GMM_tbl%fst(i)
            if (tmp_meta%l(3)%high <= 1) then
               istat = gmm_get(GMM_tbl%vname(i),tr2,tmp_meta)
               call out_fstecr (tr2, tmp_meta%l(1)%low,tmp_meta%l(1)%high,&
                                     tmp_meta%l(2)%low,tmp_meta%l(2)%high,&
                                     hyb0,nomvar,Outd_convmult(ii,set)   ,&
                                     Outd_convadd(ii,set),Level_kind_ip1 ,&
                                     -1,1,ind0,1, Outd_nbit(ii,set),.false. )
            else
               istat = gmm_get(GMM_tbl%vname(i),tr3,tmp_meta)
               NK= tmp_meta%l(3)%high ; WNK= nko
               if ((tmp_meta%l(3)%high>G_nk).and.write_diag_lev) WNK= nko+1
               call out_fstecr (tr3, tmp_meta%l(1)%low,tmp_meta%l(1)%high,&
                                     tmp_meta%l(2)%low,tmp_meta%l(2)%high,&
                                  level_type,nomvar,Outd_convmult(ii,set),&
                                      Outd_convadd(ii,set),Level_kind_ip1,&
                               -1,NK,indo,WNK, Outd_nbit(ii,set),.false. )
            end if

            cycle
         end if

      end do
      end do
!
!     ---------------------------------------------------------------
!
      return
      end subroutine out_gmm
