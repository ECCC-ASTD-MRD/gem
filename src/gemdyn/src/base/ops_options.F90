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
module Ops_options
   implicit none
   public
   save

   !# Main selector for Operational configurations combo defaults
   character(len=32) :: Ops_configuration_S = ''
   namelist /Ops_cfgs/ Ops_configuration_S

contains

!**s/r Ops_nml - Read namelist Ops_cfgs

      integer function Ops_nml (F_unf)
      use lun
      use clib_itf_mod
      implicit none

      integer, intent(in) :: F_unf

      character(len=64) :: nml_S
      logical nml_must
      integer err
!
!-------------------------------------------------------------------
!
! boiler plate - start
      if ( F_unf < 0 ) then
         Ops_nml= 0
         if ( Lun_out >= 0) write (Lun_out,nml=Ops_cfgs)
         return
      end if

      Ops_nml= -1 ; nml_must= .false. ; nml_S= 'Ops_cfgs'

      rewind(F_unf)
      read (F_unf, nml=Ops_cfgs, end= 1001, err=1003)
      Ops_nml= 0 ; goto 1000
 1001 if (Lun_out >= 0) write (Lun_out, 6005) trim(nml_S)
      if (.not.nml_must) then
         Ops_nml= 1
         if (Lun_out >= 0) write (Lun_out, 6002) trim(nml_S)
      end if
      goto 1000
 1003 if (Lun_out >= 0) write (Lun_out, 6007) trim(nml_S)

 1000 if (Ops_nml < 0 ) return
      if ((Lun_out>=0).and.(Ops_nml==0)) write (Lun_out, 6004) trim(nml_S)
      Ops_nml= 1

 6002 format (' Skipping reading of namelist ',A)
 6004 format (' Reading of namelist ',A,' is successful')
 6005 format (' Namelist ',A,' NOT AVAILABLE')
 6007 format (/,' NAMELIST ',A,' IS INVALID'/)
! boiler plate - end

      err = clib_toupper ( Ops_configuration_S )
!
!-------------------------------------------------------------------
!
      return
      end function Ops_nml
      
      subroutine Ops_defaults()
      use dynkernel_options
      use dyn_fisl_options
      use step_options
      use HORgrid_options
      use VERgrid_options
      use hvdif_options
      use init_options
      use inp_options
      use gem_options
      use lam_options
      use out_options
      use adz_options
      use sol_options
      use inp_mod
      use out_mod
      implicit none
!
!-------------------------------------------------------------------
!
      if ( Ops_configuration_S(1:15) == 'HRDPS2P5_202207' ) then

         Grd_typ_S  = 'LU'    
         Grd_ni     = 2536    ;     Grd_nj    = 1286
         Grd_dx     = 0.0225  ;    Grd_dy    = 0.0225
         Grd_lonr   = 184.961246948689;Grd_latr  = 2.19875000044703
         Grd_xlon1  = -100.00 ;     Grd_xlat1 = 53.00
         Grd_xlon2  =  -85.00 ;     Grd_xlat2 = 50.00
         Grd_maxcfl = 4       

         Step_dt        = 60.
         Fcst_gstat_S     = '1h'
         Fcst_nesdt_S     = '1h'
         Fcst_end_S     = '48'   

         dynamics_Kernel_S = 'DYNAMICS_FISL_H'
         dynamics_hydro_l  = .false.
         
         Lam_blendoro_L = .false.  
         Lam_gbpil_T=3
         Lam_blend_T=3
         
         hyb_H(1:62)= (/&
              43756.4023,  39318.1914,  35090.8477,  31060.8438,  27136.6133,& 
              24161.5703,  21995.9473,  20447.1094,  19351.3457,  18580.7734,& 
              18042.8926,  17607.1973,  17175.1387,  16743.1016,  16310.9150,& 
              15878.8105,  15446.1660,  15013.3779,  14580.7666,  14147.7617,& 
              13714.9189,  13281.5586,  12848.4268,  12414.7881,  11981.2051,& 
              11547.2939,  11114.3330,  10681.3379,  10246.2051,   9809.4922,& 
              9372.0869,   8935.4023,   8500.8115,   8069.0381,   7641.2344,& 
              7218.6914,   6802.2773,   6393.0088,   5991.8540,   5599.2866,& 
              5216.2485,   4842.9409,   4480.1753,   4128.3555,   3787.5525,& 
              3458.4089,   3140.7646,   2835.0325,   2541.2358,   2259.0730,& 
              1988.9929,   1730.5981,   1485.0991,   1254.8451,   1041.4413,& 
              845.6649,    667.7510,    507.6251,    364.7609,    238.2119,& 
              126.7142,     42.3482/)
         hyb_rcoef(1:2) = (/1., 16./)

         Sol_type_S      = 'ITERATIVE_3D'
         Sol_precond3D_S = 'RAS'
         Cstv_tstr_8     = 240.0      
         Schm_itcn       = 2          
         Schm_itnlh      = 2          
         Schm_wload_L    = .true.     
         Schm_hzdadw_L   = .true.     
         Schm_advec      = 2
         Sol_one_transpose_L = .false.

         Hzd_pwr         = 2          
         Hzd_lnr         = 0.4        
         Hzd_pwr_theta   = 6          
         Hzd_lnr_theta   = 0.02       
         
         Vspng_coeftop   = 67000.     
         Vspng_nk        = 0          

         Out3_nbitg      = 12        
         Out3_cliph_L     = .true.   
         Out3_linbot     =  3        
         
         Out3_etik_s     = 'CTRL'
         Out3_close_interval_S='1h'  
         Out3_postproc_fact=6

         Out3_lieb_levels(1:50)= (/&
              5000., 4900., 4800., 4700., 4600., 4500., 4400., 4300., 4200., 4100.,&
              4000., 3900., 3800., 3700., 3600., 3500., 3400., 3300., 3200., 3100.,&
              3000., 2900., 2800., 2700., 2600., 2500., 2400., 2300., 2200., 2100.,&
              2000., 1900., 1800., 1700., 1600., 1500., 1400., 1300., 1200., 1100.,&
              1000.,  900.,  800.,  700.,  600.,  500.,  400.,  300.,  200.,  100./)

         Out3_npes    = 80           
         Out3_npex    = -1           
         Out3_npey    = -1
         Out3_npes    = 63            
         INP_NPES=85

      endif
      
      if ( Ops_configuration_S(1:14) == 'GDPS15K_202207' ) then

         Grd_typ_S       = 'GY'       
         Grd_nj          = 683     
         Grd_xlat1       = 57.   ; Grd_xlon1       = 250.
         Grd_xlat2       = 56.   ; Grd_xlon2       = 291.
         Grd_overlap     = 2.0   ; Grd_maxcfl      = 10

         dynamics_Kernel_S = 'DYNAMICS_FISL_H'
         dynamics_hydro_l  = .true.

         Fcst_start_S = "0H"    ; Fcst_end_S = "144H"
         Fcst_gstat_S = "24h"
         Step_dt      = 450.

         adz_BC_min_max_L   = .false.
         adz_ILMC_min_max_L = .false.
         adz_slt_winds      = .true.

         Tr3d_list_S(1:2) =  (/'HU,wload=0,mono=2,mass=1,hzd=0,min=0.',&
                               'QC,wload=1,mono=2,mass=1,hzd=0,min=0.'/)

         Sol_type_S      = 'ITERATIVE_3D'
         Sol_precond3D_S = 'RAS'
         cstv_tstr_8     = 280.0
         sol_yyg_eps     = 1.e-04
         Schm_nblendyy   = 1
         Schm_itcn       = 2          
         Schm_itnlh      = 2
         Schm_hzdadw_L   = .true.
         Schm_wload_L    = .true.
         Schm_psadj      = 1

         Hzd_pwr         = 4          ; Hzd_lnr       = 0.04
         Hzd_pwr_theta   = 6          ; Hzd_lnr_theta = 0.01
         Vspng_coeftop   =  440000.   ; Vspng_nk        = 6
         Eq_sponge(1:8)  = (/50.0, 46.6, 38.3, 28.2, 18.4, 10.4, 4.6, 1.1/)
         P_lmvd_weigh_high_lat = 0.0  
         hyb_H(1:84)= (/&
 0.64765078E+05,  0.61484387E+05,  0.58097922E+05,  0.54647410E+05,  0.51260059E+05,  0.48051715E+05,  0.45091891E+05,&
 0.42408621E+05,  0.39991402E+05,  0.37818480E+05,  0.35865402E+05,  0.34107391E+05,  0.32520846E+05,  0.31080762E+05,&
 0.29757918E+05,  0.28537400E+05,  0.27406707E+05,  0.26354709E+05,  0.25371539E+05,  0.24448520E+05,  0.23581697E+05,&
 0.22768879E+05,  0.22006967E+05,  0.21292291E+05,  0.20620951E+05,  0.19988980E+05,  0.19391668E+05,  0.18824627E+05,&
 0.18284535E+05,  0.17768248E+05,  0.17272793E+05,  0.16795383E+05,  0.16333382E+05,  0.15884299E+05,  0.15445753E+05,&
 0.15015443E+05,  0.14591126E+05,  0.14171107E+05,  0.13754707E+05,  0.13341396E+05,  0.12930710E+05,  0.12522283E+05,&
 0.12115817E+05,  0.11711102E+05,  0.11307979E+05,  0.10906264E+05,  0.10502634E+05,  0.10095734E+05,  0.96856943E+04,&
 0.92727051E+04,  0.88570342E+04,  0.84390186E+04,  0.80190830E+04,  0.75977124E+04,  0.71754697E+04,  0.67529868E+04,&
 0.63318667E+04,  0.59153848E+04,  0.55064302E+04,  0.51074028E+04,  0.47201973E+04,  0.43462192E+04,  0.39864304E+04,&
 0.36414043E+04,  0.33113760E+04,  0.29963164E+04,  0.26959851E+04,  0.24104956E+04,  0.21411509E+04,  0.18889203E+04,&
 0.16543748E+04,  0.14377098E+04,  0.12387908E+04,  0.10572026E+04,  0.89229810E+03,  0.74325543E+03,  0.60785242E+03,&
 0.48407559E+03,  0.37209671E+03,  0.27309692E+03,  0.18836304E+03,  0.11854238E+03,  0.63249203E+02,  0.21093637E+02/)
         hyb_rcoef(1:2) = (/3., 15./)
         
         Init_balgm_L = .true.  ; Init_dftr_L     = .false.
         Init_dfwin_L = .true.  ; Init_dflength_S = '6h'
         Init_dfpl_S  = '6h'    ; iau_stats_l     = .true.
         
         Inp_npes = 81

         Out3_etik_s     = 'CTRL'
         Out3_lieb_levels(1:50)= (/&
               5000., 4900., 4800., 4700., 4600., 4500., 4400., 4300., 4200., 4100.,&
               4000., 3900., 3800., 3700., 3600., 3500., 3400., 3300., 3200., 3100.,&
               3000., 2900., 2800., 2700., 2600., 2500., 2400., 2300., 2200., 2100.,&
               2000., 1900., 1800., 1700., 1600., 1500., 1400., 1300., 1200., 1100.,&
               1000.,  900.,  800.,  700.,  600.,  500.,  400.,  300.,  200.,  100./)
         Out3_nbitg      = 12    ; Out3_cliph_L     = .true.
         Out3_postproc_fact= 48
         Out3_linbot       =  3
         Out3_npes         = 40

      endif
!
!-------------------------------------------------------------------
!
      return
      end subroutine Ops_defaults
      

end module Ops_options
