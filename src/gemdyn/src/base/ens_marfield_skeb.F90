!---------------------------------- LICENCE BEGIN -------------------------------
! GEM - Library of kernel routines for the GEM numerical atmospheric model
! Copyright (C) 1990-2010 - Division de Recherche en Prevision Numerique
!                          Environnement Canada
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

!**   s/r ens_marfield_skeb if ( Step_kount == 1 ) then- define a markov chain field for SKEB
!
!
      subroutine ens_marfield_skeb(fgem)
!
      use cstv
      use ens_gmm_dim
      use ens_gmm_var
      use ens_options
      use ens_param
      use glb_ld
      use gmm_itf_mod
      use HORgrid_options
      use gem_options
      use init_options
      use lun
      use path
      use ptopo
      use out_mod
      use tdpack
      use step_options
      use, intrinsic :: iso_fortran_env
      implicit none
#include <arch_specific.hf>

      real, dimension(l_ni,l_nj), intent(out) :: fgem
!
!author: Rabah Aider R.P.N-A
!
!arguments     none
!
#include <rmnlib_basics.hf>

       real,    external :: gasdev
       real(kind=REAL64)  :: pl
!
! nlat, nlon                 dimension of the Gaussian grid
! idum                       Semence du générateur de nombres aléatoires
!
      logical :: write_markov_l
      integer :: nlat, nlon, lmin,lmax
      integer :: sig, l ,m, n, i, j, dim,  indx, ier, gmmstat, gdyy
      real    :: fstd, tau , sumsp , fact, fact2, offi,offj
      real    :: xfi(l_ni),yfi(l_nj)
      real(kind=REAL64)  :: rad2deg_8, pri_8
      logical, save :: init_done=.false.
!
! paidum   pointer vers l'etat du generateur sauvegarde idum
      integer, pointer :: paiv,paiy,paiset,pagset,paidum
!
! dt   Pas de temps du modèle (secondes)
! tau  Temps de décorrélation du champ aléatoire f(i,j) (secondes)
! eps  EXP(-dt/tau/2.146)
      real(kind=REAL64)    :: dt, eps ,fmax, fmin
      real(kind=REAL64),    dimension(:), allocatable :: pspectrum , fact1, fact1n
      real(kind=REAL64),    dimension(:), allocatable  :: wrk1
      real(kind=REAL64),    dimension(:,:,:),allocatable :: cc
      real  ,    dimension(:,:),allocatable :: f
      integer :: unf0, unf1, err, errop, itstep_s, iperiod_iau
      character(len=1024) fn0, fn1
!
!---------------------------------------------------------------------
!
      dt=real(Cstv_dt_8)
      rad2deg_8=180.0d0/pi_8
      itstep_s=step_dt*step_kount
      iperiod_iau = iau_period
      write_markov_l=(itstep_s==iperiod_iau)

! Look for the spectral coefficients, legendre polymome & random numbers
!
      gmmstat = gmm_get(gmmk_ar_s,ar_s,meta2d_ar_s)
      if (GMM_IS_ERROR(gmmstat))write(*,6000)'ar_s'

      gmmstat = gmm_get(gmmk_ai_s,ai_s,meta2d_ai_s)
      if (GMM_IS_ERROR(gmmstat))write(*,6000)'ai_s'

      gmmstat = gmm_get(gmmk_br_s,br_s,meta2d_br_s)
      if (GMM_IS_ERROR(gmmstat))write(*,6000)'br_s'
      gmmstat = gmm_get(gmmk_bi_s,bi_s,meta2d_bi_s)
      if (GMM_IS_ERROR(gmmstat))write(*,6000)'bi_s'
      gmmstat = gmm_get(gmmk_dumdum_s,dumdum,meta2d_dum)
      if (GMM_IS_ERROR(gmmstat))write(*,6000)'dumdum'

      gmmstat = gmm_get(gmmk_pls_s,pls,meta3d_pls)
      if (GMM_IS_ERROR(gmmstat))write(*,6000)'pls'

!
      if (.not.init_done) then
        if (Lun_out > 0) then
        write( Lun_out,1000)
      end if
        init_done=.true.
      end if

      lmin = Ens_skeb_trnl
      lmax = Ens_skeb_trnh
      nlat = Ens_skeb_nlat
      nlon = Ens_skeb_nlon
      fstd = Ens_skeb_std
      fmin = Ens_skeb_min
      fmax = Ens_skeb_max
      tau  = Ens_skeb_tau/2.146
      eps  = exp(-dt/tau)

      if ( Step_kount == 1 ) then

! Compute Associated Legendre polynomials to use in each timestep
         pls=0.D0
         do l=lmin,lmax
            fact=DSQRT((2.D0*DBLE(l)+1.D0)/(4.D0*pi_8))
            do m=0,l
               sig=(-1.D0)**(l+m)

               do j=1,nlat/2
                  call pleg (l, m, j, nlat, pl)
                  pls(j,lmax-l+1,m+1)=pl*fact
                  pls(nlat-j+1,lmax-l+1,m+1)=pl*fact*sig
               end do
            end do
         end do
         !Initialise spectral coeffs and stochastic params
         if(Ens_recycle_mc) then
            !Read saved stochastic numbers and spectral coeffs ar,br.ai,bi
            if (ptopo_myproc==0 .and. ptopo_couleur==0) then
                  unf0=1
                  fn0= trim(Path_input_S)//'/MODEL_INPUT'//'/MRKV_SKEB.bin'
                  open ( unf0,file=trim(fn0),status='OLD', &
                             form='unformatted',iostat=errop )
                  if (errop == 0) then
                     write(output_unit,2000) 'READING', trim(fn0)
                     do i=1,36
                        read (unf0)dumdum(i,1)
                     end do
                     do l=lmin,lmax
                        do m=1,l+1
                           read(unf0)ar_s(lmax-l+1,m)
                           read(unf0)ai_s(lmax-l+1,m)
                           read(unf0)br_s(lmax-l+1,m)
                           read(unf0)bi_s(lmax-l+1,m)
                        end do
                     end do
                     close(unf0)
                  else
                     write (output_unit, 3000) trim(fn0)
                     call gem_error ( err,'read_markov_skeb', 'problem reading file' )
                  endif
            endif
            dim=ens_skeb_l*ens_skeb_m
            call RPN_COMM_bcast (dumdum,36,"MPI_INTEGER",0,"MULTIGRID", err)
            call RPN_COMM_bcast (ar_s,dim,"MPI_REAL",0, "MULTIGRID", err)
            call RPN_COMM_bcast (ai_s,dim,"MPI_REAL",0, "MULTIGRID", err)
            call RPN_COMM_bcast (br_s,dim,"MPI_REAL",0, "MULTIGRID", err)
            call RPN_COMM_bcast (bi_s,dim,"MPI_REAL",0, "MULTIGRID", err)
         else
            fstd=Ens_skeb_std
            tau = Ens_skeb_tau/2.146
            eps=exp(-dt/tau)
! Bruit blanc en nombre d'ondes
            allocate ( pspectrum(lmin:lmax) , fact1(lmin:lmax) )
            do l=lmin,lmax
               pspectrum(l)=1.D0
            end do
!Normalisation du spectre pour que la variance du champ aléatoire soit std**2
            sumsp=0.D0
            do l=lmin,lmax
               sumsp=sumsp+pspectrum(l)
            end do
               pspectrum=pspectrum/sumsp

            do l=lmin,lmax
               fact1(l)=fstd*SQRT(4.*pi/real((2*l+1))*pspectrum(l))
            end do
! Random fonction generator
            dumdum(:,1)=0
            paiv  => dumdum(1,1)
            paiy  => dumdum(33,1)
            paiset=> dumdum(34,1)
            pagset=> dumdum(35,1)
            paidum=> dumdum(36,1)
            paidum=-Ens_mc_seed
! Valeurs initiales des coefficients spectraux
            ar_s(:,:)=0.d0
            br_s(:,:)=0.d0
            ai_s(:,:)=0.d0
            bi_s(:,:)=0.d0

            do l=lmin,lmax
               br_s(lmax-l+1,1)=fact1(l)*gasdev(paiv,paiy,paiset,pagset,paidum)
               ar_s(lmax-l+1,1)=br_s(lmax-l+1,1)
               do m=2,l+1
                  br_s(lmax-l+1,m)=fact1(l)*gasdev(paiv,paiy,paiset,pagset,paidum)/sqrt(2.)
                  ar_s(lmax-l+1,m)=br_s(lmax-l+1,m)
                  bi_s(lmax-l+1,m)=fact1(l)*gasdev(paiv,paiy,paiset,pagset,paidum)/sqrt(2.)
                  ai_s(lmax-l+1,m)=bi_s(lmax-l+1,m)
               end do
            end do

            deallocate (pspectrum, fact1)
         end if
      end if

!  Begin Markov chains

! Spectrum choice
      allocate ( pspectrum(lmin:lmax) )
      allocate ( fact1(lmin:lmax),fact1n(lmin:lmax) )

      do l=lmin,lmax
         pspectrum(l)=1.D0
      end do
      sumsp=0.D0
      do l=lmin,lmax
         sumsp=sumsp+pspectrum(l)
      end do
      pspectrum=pspectrum/sumsp
      fact2 =(1.-eps*eps)/SQRT(1.+eps*eps)

! Random generator function
      paiv  => dumdum(1,1)
      paiy  => dumdum(33,1)
      pagset=> dumdum(34,1)
      paiset=> dumdum(35,1)
      paidum=> dumdum(36,1)


! Save random numbers and coefficient ar,ai,br,bi
      if (write_markov_l) then
         if (ptopo_couleur == 0  .and. ptopo_myproc == 0) then
            fn1=trim(Out_dirname_S)//'/'// 'MRKV_SKEB.bin'
            unf1=1
            open ( unf1,file=trim(fn1),status='NEW', &
                 form='unformatted',iostat=errop )
            if ( errop == 0 ) then
               write(output_unit,2000) 'WRITING', trim(fn1)
               do i=1,36
                  write (unf1)dumdum(i,1)
               end do
               do l=lmin,lmax
                  do m=1,l+1
                     write(unf1)ar_s(lmax-l+1,m)
                     write(unf1)ai_s(lmax-l+1,m)
                     write(unf1)br_s(lmax-l+1,m)
                     write(unf1)bi_s(lmax-l+1,m)
                  end do
               end do
               close(unf1)
            else
               write(output_unit,4000) 'WRITING', trim(fn1)
            endif
         endif
      endif

      do l=lmin,lmax
         fact1n(l)=fstd*SQRT(4.*pi_8/real((2*l+1)) &
                   *pspectrum(l))*SQRT((1.-eps*eps))

         br_s(lmax-l+1,1) = eps*br_s(lmax-l+1,1)  &
                          + gasdev(paiv,paiy,paiset,pagset,paidum)*fact1n(l)
         ar_s(lmax-l+1,1) = eps*ar_s(lmax-l+1,1)  + br_s(lmax-l+1,1)*fact2
         do m=2,l+1
            br_s(lmax-l+1,m) = eps*br_s(lmax-l+1,m) &
                             + gasdev(paiv,paiy,paiset,pagset,paidum)*fact1n(l)/SQRT(2.)
            ar_s(lmax-l+1,m) = eps*ar_s(lmax-l+1,m)+br_s(lmax-l+1,m)*fact2
            bi_s(lmax-l+1,m) = eps*bi_s(lmax-l+1,m) &
                             + gasdev(paiv,paiy,paiset,pagset,paidum)*fact1n(l)/SQRT(2.)
            ai_s(lmax-l+1,m) = eps*ai_s(lmax-l+1,m)+bi_s(lmax-l+1,m)*fact2
         end do
      end do
      deallocate (pspectrum, fact1, fact1n)

      allocate(cc(2 , nlat, lmax+1))
      allocate(wrk1( nlat * (nlon+2)))
      allocate(f(nlon, nlat) )

      cc(1:2,1:nlat,1:lmax+1)=0.D0

      do m=1,lmax+1
         do j=1,nlat
            cc(1,j,m)=0.d0
            cc(2,j,m)=0.d0
            cc(1,j,m)=cc(1,j,m) + Dot_product(pls(j,1:lmax-lmin+1,m),ar_s(1:lmax-lmin+1,m))
            cc(2,j,m)=cc(2,j,m) + Dot_product(pls(j,1:lmax-lmin+1,m),ai_s(1:lmax-lmin+1,m))
         end do
      end do

! Fourier Transform (inverse)

      wrk1(:)=0.0
      n=-1
      do i=1,nlat
         do j=1,lmax+1
            n = n + 2
            wrk1(n)   = cc(1,i,j)
            wrk1(n+1) = cc(2,i,j)
         end do
         n=n+nlon-2*lmax
      end do

      call itf_fft_set(nlon,'PERIODIC',pri_8)
      call itf_fft_drv(wrk1,1,nlon+2,nlat,1)
      n=0
      do j=1,nlat
         do i=1,nlon+2
            n = n + 1
            if (i <= nlon) then
               f(i,j) = wrk1(n)
            end if
         end do
      end do

      deallocate(cc,wrk1)

!*    Interpolation to the processors grids

      offi = Ptopo_gindx(1,Ptopo_myproc+1)-1
      offj = Ptopo_gindx(3,Ptopo_myproc+1)-1

      do i=1,l_ni
         indx = offi + i
         xfi(i) = G_xg_8(indx)*rad2deg_8
      end do
      do i=1,l_nj
         indx = offj + i
         yfi(i) = G_yg_8(indx)*rad2deg_8
      end do

      gdyy = ezqkdef(nlon,nlat,'A', 0,0,0,0,0)

         ier = ezdefset(Grd_local_gid, gdyy)
         ier = ezsetopt('INTERP_DEGREE', 'LINEAR')
         ier = ezsint(fgem,f)

      if(Ens_stat)then
         call glbstat (fgem,'MCSK','',&
           1,l_ni,1,l_nj,1,1,1,G_ni,1,G_nj,1,1)
      end if

      deallocate(f)



 1000 format( &
           /,'INITIALIZE SCHEMES CONTROL PARAMETERS (S/R ENS_MARFIELD_SKEB)', &
           /,'======================================================')

 2000 format (/' MARKOV: ',a,' FILE ',a)

 3000 format(' S/R  ENS_MARFIELD_SKEB : problem in opening MRKV_PARAM_SKEB file)', &
              /,'======================================================')

 4000 format(' S/R  ENS_MARFIELD_SKEB : problem with registering MRKV_PARAM_SKEB file)', &
              /,'======================================================')

 6000 format('ens_marfield_skeb at gmm_get(',A,')')

      return

contains

 subroutine pleg(l, m, jlat, nlat, pls )
      use, intrinsic :: iso_fortran_env
      implicit none

      integer l,m ,i,j ,jlat ,nlat
      real(kind=REAL64)   pls
      real(kind=REAL64)  factor , x  ,lat, theta
      real(kind=REAL64) , dimension(0:l+1) :: pl
      real(kind=REAL64), parameter :: ZERO=0.0D0  , ONE_8=1.0d0 , TWO_8=2.0d0

!-------------------------------------------------------------------------
!
      if ( m < 0 .OR. m > l ) then
         print*, ' error :  m must non-negative and m <=l '
         stop
      end if

      lat=(-90.D0+90.D0/DBLE(nlat)+DBLE(jlat-1)*180.D0/DBLE(nlat))*pi_8/180.D0
      theta=pi_8/2.D0-lat
      x=DCOS(theta)

      pl=ZERO
      if ( m <= l ) then
         pl(m) = ONE_8
         factor = ONE_8
         do i = 1, m
            pl(m) = -pl(m)*factor*sqrt(1.d0 - x**2)/ &
                   dsqrt(dble((l+i)*(l-m+i)))
            factor = factor + 2.d0
         end do
         pls=pl(m)
      end if

      if ( m + 1 <= l ) then
         pls = x * dble ( 2 * m + 1 ) * pl(m)
         pl(m+1)=pls
      end if

      do j = m + 2, l
          pl(j) = ( x * dble (2*j-1) * pl(j-1) &
                     + dble (-j-m+1) * pl(j-2) ) &
                     / dble (j-m)
      end do

      pls=pl(l)

  end subroutine pleg

end subroutine ens_marfield_skeb

