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

      subroutine ilmc_lam_loop ()

      do k=Adz_k0t,l_nk

         do j=jl,jr
         do i=il,ir

            o_shoot = max(0., F_ilmc(i,j,k) - Adz_post(F_ns)%max(i,j,k))
            u_shoot = max(0., Adz_post(F_ns)%min(i,j,k) - F_ilmc(i,j,k))

            if (o_shoot==0..and.u_shoot==0.) cycle

            !Positive violation
            !------------------
            if (o_shoot>0.) then

               do sweep = 1,Adz_ILMC_sweep_max

                  m_sweep_max = 0.

                  do ii = 1,sweep_rd(sweep,i,j,k)%cell

                     i_rd = sweep_rd(sweep,i,j,k)%i_rd(ii)
                     j_rd = sweep_rd(sweep,i,j,k)%j_rd(ii)
                     k_rd = sweep_rd(sweep,i,j,k)%k_rd(ii)

                     sweep_max(ii) = max(0., Adz_post(F_ns)%max(i_rd,j_rd,k_rd) -  F_ilmc(i_rd,j_rd,k_rd))
                     m_sweep_max   = sweep_max(ii) * air_mass_m(i_rd,j_rd,k_rd) + m_sweep_max

                  end do

                  m_use_max = min(m_sweep_max,o_shoot*air_mass_m(i,j,k))

                  !Add Mass to surrounding cells up to Local Maximum
                  !-------------------------------------------------
                  if (m_use_max>0.) then

                     reset(k,5) = reset(k,5) + 1

                     F_ilmc(i,j,k) = max(F_ilmc(i,j,k) - real(m_use_max/air_mass_m(i,j,k)),Adz_post(F_ns)%min(i,j,k))

                     o_shoot = max(0., F_ilmc(i,j,k) - Adz_post(F_ns)%max(i,j,k))

                     ratio_max = m_use_max/m_sweep_max

                     do ii = 1,sweep_rd(sweep,i,j,k)%cell

                        i_rd = sweep_rd(sweep,i,j,k)%i_rd(ii)
                        j_rd = sweep_rd(sweep,i,j,k)%j_rd(ii)
                        k_rd = sweep_rd(sweep,i,j,k)%k_rd(ii)

                        if (sweep_max(ii)>0.) then

                           F_ilmc(i_rd,j_rd,k_rd) = min(F_ilmc(i_rd,j_rd,k_rd) + &
                                                    real(sweep_max(ii)*ratio_max),Adz_post(F_ns)%max(i_rd,j_rd,k_rd))

                        end if

                     end do

                  end if

               end do !DO sweep u_shoot

            !Negative violation
            !------------------
            else

               do sweep = 1,Adz_ILMC_sweep_max

                  m_sweep_min = 0.

                  do ii = 1,sweep_rd(sweep,i,j,k)%cell

                     i_rd = sweep_rd(sweep,i,j,k)%i_rd(ii)
                     j_rd = sweep_rd(sweep,i,j,k)%j_rd(ii)
                     k_rd = sweep_rd(sweep,i,j,k)%k_rd(ii)

                     sweep_min(ii) = max(0., F_ilmc(i_rd,j_rd,k_rd) - Adz_post(F_ns)%min(i_rd,j_rd,k_rd))
                     m_sweep_min   = sweep_min(ii) * air_mass_m(i_rd,j_rd,k_rd) + m_sweep_min

                  end do

                  m_use_min = min(m_sweep_min,u_shoot*air_mass_m(i,j,k))

                  !Subtract Mass from surrounding cells up to Local Minimum
                  !--------------------------------------------------------
                  if (m_use_min>0.) then

                     reset(k,5) = reset(k,5) + 1

                     F_ilmc(i,j,k) = min(F_ilmc(i,j,k) + real(m_use_min/air_mass_m(i,j,k)),Adz_post(F_ns)%max(i,j,k))

                     u_shoot = max(0., Adz_post(F_ns)%min(i,j,k) - F_ilmc(i,j,k))

                     ratio_min = m_use_min/m_sweep_min

                     do ii = 1,sweep_rd(sweep,i,j,k)%cell

                        i_rd = sweep_rd(sweep,i,j,k)%i_rd(ii)
                        j_rd = sweep_rd(sweep,i,j,k)%j_rd(ii)
                        k_rd = sweep_rd(sweep,i,j,k)%k_rd(ii)

                        if (sweep_min(ii)>0.) then

                           F_ilmc(i_rd,j_rd,k_rd) = max(F_ilmc(i_rd,j_rd,k_rd) - &
                                                    real(sweep_min(ii)*ratio_min),Adz_post(F_ns)%min(i_rd,j_rd,k_rd))

                        end if

                     end do

                  end if

               end do !DO sweep u_shoot

            end if

         end do
         end do

      end do

      end subroutine ilmc_lam_loop
