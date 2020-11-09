!-------------------------------------- LICENCE BEGIN ------------------------
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
!-------------------------------------- LICENCE END --------------------------

!/@*
function phydebu2(p_ni, p_nj, p_nk, F_path_S) result(F_istat)
   use iso_c_binding
   use rpn_comm_itf_mod
   use phy_status, only: PHY_ERROR, PHY_OK, phy_error_L
   use phy_options
   use phybus, only: entbus, perbus, dynbus, volbus
   use phybusalloc, only: phybusalloc1
   use module_mp_p3, only: p3_init
   use ghg_mod, only: ghg_init
   use mixing_length, only: ML_CLOSURES
   use ens_perturb, only: ens_spp_map, ENS_OK
   implicit none
!!!#include <arch_specific.hf>
   !@Object Init physics at the beginning of each execution of the model
   !@Arguments
   character(len=*), intent(in) :: F_path_S !# data/tables dir
   integer, intent(in) :: p_ni, p_nj, p_nk  !# horiz and vert dimensions
   !@return
   integer :: F_istat
   !@Author B. Bilodeau (Spring 1994)
   !@Revisions
   ! 001      M. Gagnon   (Jul 95) - Added validation code for radniv
   ! 002      M. Desgagne (Nov 95) - Unified physics interface
   ! 003      L. Lefaivre (Nov 95/Feb 96) - Initialize ETRMIN and Z0MIN
   !                                 with values passed from dynamicsa
   ! 004      B. Dugas (Sep 96) - Coherence check between CLIMAT, RADFIX
   ! 005      G. Pellerin (Nov 95) - Added switches for deep convection
   !                                 KUOSTD,KUOSYM,KUOSUN with CONSUN
   ! 006      G. Pellerin (Nov 96) - Insert common tables for RAS option
   ! 007      B. Bilodeau (Apr 97) - Insert comdeck for CLASS
   ! 008      M. Desgagne (Spring 97) - Microphysics
   ! 009      B. Bilodeau (Jan 98) - Connect FCP and KFC with CONSUN
   ! 010      Y. Delage (Feb 98) - Addition of HMIN in "surfcon.cdk"
   ! 011      B. Bilodeau (Jun 98) - RADFILES and FOMICHEV
   ! 012      M. Desgagne (Oct 98) - call back to rdradf_d (from dynamics)
   ! 013      B. Bilodeau (Dec 98) - New "entry" bus
   ! 014      M. Desgagne (Dec 98) - Correct bug in calculation of moyhr
   ! 015      B. Bilodeau (Oct 99) - CW_RAD
   ! 016      B. Bilodeau (Oct 2000) - Move consistency tests at the end
   !                         of the subroutine to correct FOMIC-REDUC bug
   ! 017      B. Bilodeau (Nov 2000) - Replace call to radini, turbini,
   !                                   gwdini and convini by call to phy_ini.
   !                                   Eliminate call to ptcalc.
   ! 018      S. Belair and B. Bilodeau (May 2001)
   !                                 - New density for fresh snow.
   ! 019      B. Bilodeau (Mar 2001) - OPTIX
   ! 020      B. Dugas (Jan 2002) - FOMIC and REDUC are now compatible
   ! 021      B. Bilodeau (Mar 2002) - Correct bug in calculation of nspliti
   !                                   and add dzsedi.cdk
   ! 022      A-M. Leduc (Jan 2003)  - SHLCVT becomes SHLCVT(1) or SHLCVT(2)
   ! 023      B. Dugas (Feb 2003)    - share small_sedimentation_dt and
   !                                   cldopt_mode comdecks with SAVE_OPTIONS
   ! 024      B. Bilodeau (Feb 2003) - AS2, BETA2 and KKL2 parameters
   !                                   Remove ALAT and BLAT
   ! 025      B. Dugas (Mar 2003)    - Add STRATOS parametre
   ! 026      A. Plante (June 2003)  - Add VARMTN (mountains.cdk)
   ! 027      B. Dugas (July 2003)   - Add CRITLAC parametre
   ! 028      A. Plante (sep 2003)   - Add key pcptype rule
   ! 029      Y. Delage (Apr 2004)   - Reactivate land surface module CLASS
   !                                 - Default values of parameters in common SURFCON
   !                                    now defined in surfcon_ini.cdk
   ! 030      B. Bilodeau (Jul 2004) - Add Z0TLAT
   ! 031      L. Spacek (Aug 2004)   - cloud clean-up
   !                                   elimination of ISTCOND=2,6,7,8 ICONVEC=4
   ! 032      B. Bilodeau (Oct 2004) - Add protective code for dzsedi
   !
   ! 033      S. Valcke (Apr 2005)   - COUPLING and IMPFLX incompatible
   ! 034      B. Dugas (Aug 2005)    - Initialize commons in Block DATA PHYDEBU4_DATA
   ! 035      D. Talbot (may 2006)   - Add option cccmarad
   ! 036      J. Cole  (May 2006)    - Implement the ISCCP cloud simulator
   ! 037      B. Dugas (Dec 2006)    - Remove all reference to DEBUT
   ! 038      J. Milbrandt (Dec 2006) - Added options for 5 versions of Milbrandt-Yau scheme
   ! 039      M. Desgagne (July 2006) - Revised interface. Change name to phy_debu.
   ! 039      B. Bilodeau (Feb 2007) - Cleanup and creation of check_options
   ! 040      M. Desgagne (Mar 2008) - optional ozone file
   ! 041      A-M. Leduc  (Feb 2009) - add TRIGLAT
   ! 042      A. Plante (Sept. 2011) - allocate standard Pressure profil std_p_prof for non_oro gwd.
   ! 043      J. Milbrandt (Mar 2015) - add call to P3_INIT, for initialization for P3 microphysics
   !@Notes
   !          phy_debu does the following :
   !          1) it initializes a few constants necessary
   !             for the execution of the physics package.
   !          2) it reads the radiation files if necessary.
   !          3) it constructs the 3 main buses dictionaries.
   !*@/
#include <msg.h>
#include <rmnlib_basics.hf>

   include "clefcon.cdk"
   include "machcon.cdk"

   logical, save :: okinit = .false.

   character(len=1024) :: fichier, path
   integer :: i, nv, myproc, ier
   !---------------------------------------------------------------------
   F_istat = PHY_ERROR

   if (jdateo == 0 .or. delt == 0.) then
      call msg(MSG_ERROR,'(phydebu) VARIABLES: jdateo,delt NOT INITIALIZED')
      return
   endif

   ! INITIALISATION DE VARIABLES POUR CLEF
   ! - - - - - - - - - - - - - - - - - - -
   ETRMIN = ETRMIN2
   EXPLIM = 75.
   TANLIM = exp(12. * ALOG(2.))

   ! CONSTANTES NUMERIQUES DANS LA FERMETURE DU MODELE CLEF
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - -
   !     REF : THERRY ET LACARRERE
   !           ANDRE ET AL.
   !           BOUGEAULT
   !           MAILHOT ET BENOIT , JAS 1982
   !           WYNGAARD ET AL.
   CLEFC1 = 3.75/1.75
   CLEFC4 = 4.5
   CLEFC6 = 4.85
   CLEFC7 = 1.0-0.125*CLEFC6
   CLEFC8 = 6.5
   CLEFCB = 0.4
   CLEFAE = 3.0*CLEFC4/CLEFC8

   ! RADIATION
   ! - - - - -
   ! reduction des niveaux pour les calculs radiatifs
   ! pour option "NEWRAD"
   NV = 0
   do i=1,LEVMAX
      if (radnivl(i) == 0) exit
      if (radnivl(i) < 0 .or. radnivl(i) > p_nk) then
         call msg(MSG_ERROR,'(phydebu) radnivl - WRONG LEVELS IN THE LIST')
         return
      endif
      NV = NV+1
   enddo

   if (NV > 0 .and. NV < p_nk) then
      do i=NV,1,-1
         radnivl(i+1) = radnivl(i)
      enddo
      radnivl(1) = NV
   else if (NV == p_nk) then
      radnivl(1) = 0
   else if (NV > p_nk) then
      call msg(MSG_ERROR,'(phydebu) radnivl - TOO MANY LEVELS IN THE LIST')
      return
   endif

   if (RADNIVL(1) == 0 .or. RADNIVL(1) == P_NK) then
      REDUC = .false.
   else
      REDUC = .true.
      if( radnivl(2) /= 1 ) then
         call msg(MSG_ERROR,'(phydebu) radnivl - THE LIST MUST BEGIN WITH LEVEL 1')
         return
      endif
      do i=2,radnivl(1)
         if( radnivl(i) >= radnivl(i+1) ) then
            call msg(MSG_ERROR,'(phydebu) radnivl - THE LIST MUST BE SORTED IN ASCENDING ORDER')
            return
         endif
      enddo
   endif
 
   !# lecture des tableaux de radiation
   IF_RADIA: if (any(radia(1:8) == (/'NEWRAD  ','CCCMARAD'/))) then
      if (.not.okinit) then

         call rpn_comm_rank(RPN_COMM_GRID, myproc, ier)
         fichier = trim(F_path_S)//'ozone_clim.fst'
         call litozon(fichier, myproc)

         fichier = trim(F_path_S) //'/rad_table.fst'
         call litblrad(fichier, myproc)
         if (phy_error_L) return

         if (simisccp) then
            phy_error_L = .true.
            !# compute the table needed to generate variability
            !  for the  stochastic cloud generator
!!$            call TABULATE_XCW()

            !# read in data blocks for ISCCP simulator code
!!$            call READ_ISCCPDATA()
            if (phy_error_L) return
         endif

         !# read GHG concentration factor file
         path = trim(F_path_S)//'/CLIMATO' !#ghg-table-1950-2015_v1'
         ier = ghg_init(path, jdateo, myproc)
         if (.not. RMN_IS_OK(ier)) then
            call msg(MSG_ERROR,'(phydebu) Problem in ghg_init')
            return
         endif

         okinit = .true.
      endif
   endif IF_RADIA

   ! MICROPHYSICS (preliminary calculations for P3 scheme)
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - -
   if (stcond == 'MP_P3') then
      path = trim(F_path_S)//'MODEL_INPUT/'
      call p3_init(path, p3_ncat, stat=ier)
      if (ier < 0) then
         call msg_toall(MSG_ERROR,'(phydebu) Problem in p3_init')
         return
      endif
   endif

   ! Set mixing length index mapping for SPP
   if (ens_spp_map('longmel', ML_CLOSURES(:)%name, ML_CLOSURES(:)%key) /= ENS_OK) then
      call msg_toall(MSG_ERROR,'(phydebu) Problem mapping mixing length names')
      return
   endif
   
   ! CONSTRUCTION OF THE 4 MAIN BUSES DICTIONARIES:
   ! BUSENT, BUSDYN, BUSPER and BUSVOL
   ! - - - - - - - - - - - - - - - - -
   call phybusinit(p_ni, p_nk)
   if (phy_error_L) then
      call msg(MSG_ERROR,'(phydebu) Problem in phybusinit')
      return
   endif

   nullify(entbus, perbus, dynbus, volbus)
   ier = phybusalloc1(p_nj, entbus, perbus, dynbus, volbus)
   if (phy_error_L .or. .not.RMN_IS_OK(ier)) return

   ! moyhr (acchr) est la periode de moyennage (accumulation) des diagnostics.
   ! conversion : nombre d'heures --> nombre de pas de temps.

   moyhr = nint(moyhr * 3600./delt)
   acchr = nint(acchr * 3600./delt)

   F_istat = PHY_OK
   !----------------------------------------------------------------------
   return
end function phydebu2
