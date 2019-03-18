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

      subroutine ilmc_lam_loop_core_n ()

      do shift=1,SIDE_max_0

!$omp parallel                                             &
!$omp private(k,ii,o_shoot,u_shoot,ratio_max,ratio_min,    &
!$omp         m_sweep_max,m_sweep_min,sweep_max,sweep_min, &
!$omp         m_use_max,m_use_min,i_rd,j_rd,k_rd,i,j,sweep)&
!$omp shared (sweep_rd,shift,il_c,ir_c,jl_c,jr_c,new,      &
!$omp         air_mass_m)
!$omp do

      do k=F_k0+shift-1,F_nk,SIDE_max_0

         do j=jl_c,jr_c
         do i=il_c,ir_c

            o_shoot = max(0.,   new(i,j,k) - F_max(i,j,k))
            u_shoot = max(0., F_min(i,j,k) -   new(i,j,k))

            if (o_shoot==0..and.u_shoot==0.) cycle

            do sweep = 1,F_ILMC_sweep_max

               if (o_shoot==0..and.u_shoot==0.) cycle

               m_sweep_max = 0.
               m_sweep_min = 0.

               do ii = 1,sweep_rd(sweep,i,j,k)%cell

                  i_rd = sweep_rd(sweep,i,j,k)%i_rd(ii)
                  j_rd = sweep_rd(sweep,i,j,k)%j_rd(ii)
                  k_rd = sweep_rd(sweep,i,j,k)%k_rd(ii)

                  if (     ( (i_rd>=il_c.and.i_rd<=ir_c) .and. (j_rd>=jl_c.and.j_rd<=jr_c) ) ) cycle

                  sweep_max(ii) = max(0., F_max(i_rd,j_rd,k_rd) -   new(i_rd,j_rd,k_rd))
                  m_sweep_max   = sweep_max(ii) * air_mass_m(i_rd,j_rd,k_rd) + m_sweep_max

                  sweep_min(ii) = max(0.,   new(i_rd,j_rd,k_rd) - F_min(i_rd,j_rd,k_rd))
                  m_sweep_min   = sweep_min(ii) * air_mass_m(i_rd,j_rd,k_rd) + m_sweep_min

               end do

               m_use_max = min(m_sweep_max,o_shoot*air_mass_m(i,j,k))
               m_use_min = min(m_sweep_min,u_shoot*air_mass_m(i,j,k))

               !Positive violation: Add Mass to surrounding cells up to Local Maximum
               !---------------------------------------------------------------------
               if (m_use_max>0.) then

                   reset(k,3) = reset(k,3) + 1

                   new(i,j,k) = max(new(i,j,k) - real(m_use_max/air_mass_m(i,j,k)),F_min(i,j,k))

                   ratio_max = m_use_max/m_sweep_max

                   do ii = 1,sweep_rd(sweep,i,j,k)%cell

                      i_rd = sweep_rd(sweep,i,j,k)%i_rd(ii)
                      j_rd = sweep_rd(sweep,i,j,k)%j_rd(ii)
                      k_rd = sweep_rd(sweep,i,j,k)%k_rd(ii)

                      if (     ( (i_rd>=il_c.and.i_rd<=ir_c) .and. (j_rd>=jl_c.and.j_rd<=jr_c) ) ) cycle

                      if (sweep_max(ii)>0.) new(i_rd,j_rd,k_rd) = min(new(i_rd,j_rd,k_rd) + real(sweep_max(ii)*ratio_max),F_max(i_rd,j_rd,k_rd))

                   end do

               !Negative violation: Subtract Mass from surrounding cells up to Local Minimum
               !----------------------------------------------------------------------------
               elseif (m_use_min>0.) then

                   reset(k,3) = reset(k,3) + 1

                   new(i,j,k) = min(new(i,j,k) + real(m_use_min/air_mass_m(i,j,k)),F_max(i,j,k))

                   ratio_min = m_use_min/m_sweep_min

                   do ii = 1,sweep_rd(sweep,i,j,k)%cell

                      i_rd = sweep_rd(sweep,i,j,k)%i_rd(ii)
                      j_rd = sweep_rd(sweep,i,j,k)%j_rd(ii)
                      k_rd = sweep_rd(sweep,i,j,k)%k_rd(ii)

                      if (     ( (i_rd>=il_c.and.i_rd<=ir_c) .and. (j_rd>=jl_c.and.j_rd<=jr_c) ) ) cycle

                      if (sweep_min(ii)>0.) new(i_rd,j_rd,k_rd) = max(new(i_rd,j_rd,k_rd) - real(sweep_min(ii)*ratio_min),F_min(i_rd,j_rd,k_rd))

                   end do

               endif

               o_shoot = max(0.,   new(i,j,k) - F_max(i,j,k))
               u_shoot = max(0., F_min(i,j,k) -   new(i,j,k))

               if (o_shoot==0..and.u_shoot==0.) cycle

            end do ! Do sweep

            if (o_shoot==0..and.u_shoot==0.) cycle

         end do
         end do

      end do
!$omp end do
!$omp end parallel

      end do ! do shift

      end subroutine ilmc_lam_loop_core_n
