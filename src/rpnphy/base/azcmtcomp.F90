!-------------------------------------- LICENCE BEGIN -------------------------
!Environment Canada - Atmospheric Science and Technology License/Disclaimer,
!                     version 3; Last Modified: May 7, 2008.
!This is free but copyrighted software; you can use/redistribute/modify it under the terms
!of the Environment Canada - Atmospheric Science and Technology License/Disclaimer
!version 3 or (at your option) any later version that should be found at:
!http://collaboration.cmc.ec.gc.ca/science/rpn.comm/license.html
!
!This software is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
!without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!See the above mentioned License/Disclaimer for more details.
!You should have received a copy of the License/Disclaimer along with this software;
!if not, you can write to: EC-RPN COMM Group, 2121 TransCanada, suite 500, Dorval (Quebec),
!CANADA, H9P 1J3; or send e-mail to service.rpn@ec.gc.ca
!-------------------------------------- LICENCE END ---------------------------

!/@*
subroutine azcmtcomp(nstep, dt2, mf, u0, uu, ud, w, adiag, aodiag, &
     bup, bdown, tum_u, tum_d, tum_w)
   !@Object 
   !@Arguments
   !#TODO: specify the intent(in), intent(out) or intent(inout)
   integer, intent(in) :: nstep, mf
   real, intent(in)  :: dt2
   real, intent(in)  :: u0(mf), uu(mf), ud(mf), w(mf), adiag(mf), aodiag(mf), bup(mf), bdown(mf)
   real, intent(out) :: tum_u(mf), tum_d(mf), tum_w(mf)

   !@Author P.Vaillancourt, 2021
   !*@/

   integer :: nk, n, kpm
   real :: bu(mf), bd(mf), aterm(mf), um(mf)
   real :: tum(mf), tum_n(mf), tum_un(mf), tum_dn(mf), tum_wn(mf)
   real :: tterm(mf), uterm(mf), dterm(mf), wterm(mf)

   !# Same thing (i.e. tendency partition) but in explicit loop form
   !#
   !# Recall that:
   !#   nstep      = # of substeps                 (e.g. 6)
   !#   dt2        = sub-timestep                  (e.g. 600s)
   !#   f          = number of active levels       (e.g. 45)
   !#   u0(nk)     = initial profile of U          (e.g. u00)
   !#   uu(nk)     = updraft profile               (e.g. uu)
   !#   ud(nk)     = downdraft profile             (e.g. ud)
   !#   w(nk)      = omega                         (e.g. omga)
   !#   Adiag(nk)  = diagonal part of matrix A     (e.g. ct1)
   !#   Aodiag(nk) = off-diagonal part of matrix A (e.g. ct2)
   !#   Bup(nk)    = updraft term                  (e.g. ct3)
   !#   Bdown(nk)  = downdraft term                (e.g. ct4)
   !#   Bu(nk)     = Bup(nk)  *Adiag(nk)*uu(nk)
   !#   Bd(nk)     = Bdown(nk)*Adiag(nk)*ud(nk)

   do nk = 1,mf
      Bu(nk) = Bup(nk)*Adiag(nk)*uu(nk)
      Bd(nk) = Bdown(nk)*Adiag(nk)*ud(nk)
   enddo

   DO_NSTEP: do n=1,nstep
      IF_N1: if (n == 1) then !# if calculating things for 1st substep

         do nk = 1,mf

            !# Use tridiagonal matrix A to compute transport term = A * u0
            kpm = max(1, min(mf, nk + int(sign(1., w(nk)))))  !# neighbor index
            Aterm(nk) = u0(nk)*Adiag(nk) + u0(kpm)*Aodiag(nk) !# transport term

            um(nk)    = Aterm(nk) + Bu(nk) + Bd(nk)      !# U+ after 1 substep

            tum(nk)   = (um(nk) - u0(nk))/dt2            !# 1-substep total tendency
            tum_u(nk) = Bup(nk)  *(uu(nk) - um(nk))/dt2  !# 1-substep updraft tendency
            tum_d(nk) = Bdown(nk)*(ud(nk) - um(nk))/dt2  !# 1-substep downdraft tendency
            tum_w(nk) = tum(nk) - tum_u(nk) - tum_d(nk)  !# 1-substep transport tendency

            if (nstep > 1) then
               !# if more than 1 substep, initialize some work arrays 
               !# using tendencies from 1st substep, 
               !# to be used in 2nd substep calculations
               tum_n(nk)  = tum(nk)   !# work array (1st contribution) for total tendency
               tum_un(nk) = tum_u(nk) !# work array (1st contribution) for updraft tendency  
               tum_dn(nk) = tum_d(nk) !# work array (1st contribution) for downdraft tendency
               tum_wn(nk) = tum_w(nk) !# work array (1st contribution) for transport tendency
            endif

         enddo

      else  !# IF_N1 - if calculating things for follwing substeps (when n > 1)

         !# Use tridiagonal matrix A to propagate tendency terms:
         !# n-th tendency term = A * ( (n-1)-th tendency term )

         do nk = 1,mf
            kpm = max(1, min(mf, nk + int(sign(1., w(nk)))))    !# neighbor index
            tterm(nk) = tum_n(nk) *Adiag(nk) + tum_n(kpm) *Aodiag(nk) !# total term
            uterm(nk) = tum_un(nk)*Adiag(nk) + tum_un(kpm)*Aodiag(nk) !# updraft term
            dterm(nk) = tum_dn(nk)*Adiag(nk) + tum_dn(kpm)*Aodiag(nk) !# downdraft term
            wterm(nk) = tum_wn(nk)*Adiag(nk) + tum_wn(kpm)*Aodiag(nk) !# transport term
         enddo

         !# NOTE: the loop above and the loop below should NOT be merged !

         do nk = 1,mf  
            tum_n(nk)  = tterm(nk)              !# n-th contribution and 
            tum(nk)    = tum(nk)   + tum_n(nk)  !# accumulator for total tendency

            tum_un(nk) = uterm(nk)              !# n-th contribution and 
            tum_u(nk)  = tum_u(nk) + tum_un(nk) !# accumulator for updraft tendency

            tum_dn(nk) = dterm(nk)              !# n-th contribution and 
            tum_d(nk)  = tum_d(nk) + tum_dn(nk) !# accumulator for downdraft tendency

            tum_wn(nk) = wterm(nk)              !# n-th contribution and 
            tum_w(nk)  = tum_w(nk) + tum_wn(nk) !# accumulator for transport tendency
         enddo

      endif IF_N1
   enddo DO_NSTEP

   !# finally, divide accumulators by nstep to get tendencies over tauc

   do nk = 1,mf
      tum(nk)   = tum(nk)  /nstep  !# final total tendency
      tum_u(nk) = tum_u(nk)/nstep  !# final updraft tendency 
      tum_d(nk) = tum_d(nk)/nstep  !# final downdraft tendency
      tum_w(nk) = tum_w(nk)/nstep  !# final transport tendency
   enddo

   return
end subroutine azcmtcomp

