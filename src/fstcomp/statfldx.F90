!> Compute the average, variance, minimum and maxium of the provided
!> field and prints the result.
!> @param nomvar Used to provide extra info on the printout
!> @param typvar Unused!
!> @param ip1 IP1 of the field f
!> @param ip2 IP2 of the field f.  Only used as extra info on the printout.
!> @param ip3 IP3 of the field f.  Only used as extra info on the printout.
!> @param date Unused!
!> @param etiket Only used as extra info on the printout.
!> @param f Field for which we wish to compute statistics
!> @param ni Number of points of the field on the X axis
!> @param nj Number of points of the field on the Y axis
!> @param nk Number of vertical levels of the field
      subroutine statfldx(nomvar, typvar, ip1, ip2, ip3, date, etiket, &
         f, ni, nj, nk)

      ! Part of Fortran 2003 Standard
      use IEEE_ARITHMETIC
      implicit none

      character*4, intent(in) :: nomvar
      character*2, intent(in) :: typvar
      ! ip1 muist have an inout intent since it's used to call convip_plus
      ! which expects it to be inout
      integer, intent(inout) :: ip1
      integer, intent(in) ::  ip2, ip3, date
      character*12, intent(in) :: etiket
      integer, intent(in) :: ni, nj, nk
      real, intent(in) :: f(ni, nj, nk)

      ! Local variables
      integer i, j, k
      real*8 sum, moy, var, rmin, rmax
      integer imin, jmin, kmin, imax, jmax, kmax, kind, dat2, dat3
      character*15 Level
      real rlevel
! --------------------------------------------------------------------

      ! Calcul de la moyenne
      sum = 0.0
      do k = 1,nk
         do j = 1, nj
            do i = 1, ni
               sum = sum + f(i, j, k)
            end do
         end do
      end do
      moy = sum / float(ni * nj * nk)

      ! Calcul de la variance
      sum = 0.0
      do k = 1, nk
         do j = 1, nj
            do i = 1, ni
               sum = sum + ((f(i,j,k) - moy)*(f(i,j,k) - moy))
            end do
         end do
      end do
      var = sqrt (sum / float(ni * nj * nk))

      ! Trouver le minimum et le maximum.
      imin = 1
      jmin = 1
      kmin = 1
      imax = 1
      jmax = 1
      kmax = 1
      rmax = f(1, 1, 1)
      rmin = f(1, 1, 1)

      do k = 1, nk
         do j = 1, nj
            do i = 1, ni
               if (f(i,j,k) .gt. rmax) then
                  rmax  = f(i,j,k)
                  imax = i
                  jmax = j
                  kmax = k
               endif
               if (f(i,j,k) .lt. rmin) then
                  rmin  = f(i,j,k)
                  imin = i
                  jmin = j
                  kmin = k
               endif
            end do
         end do
      end do

      CALL convip_plus(ip1, rlevel, kind, -1, level, .true.)

      write(6,10) nomvar, etiket, level, ip2, ip3, &
           moy, var, imin, jmin+(kmin-1)*nj, rmin, &
           imax, jmax+(kmax-1)*nj, rmax

 10   format ('  <', a4, '>', 1x, a12, a15, 1x, i8, 1x, i8, 1x, &
          ' Mean:', e15.8, ' Stdev:', e15.8, &
          '  Min:[(', i4, ',', i4, '):', &
          e11.4, ']', ' Max:[(', i4, ',', i4, '):', &
          e11.4, ']')

      ! On essaie de detecter la presence de Nan
      if (.not. IEEE_IS_NORMAL(moy)) then
         print *, '**** NaN detected'
      end if

      do k = 1, nk
         do j = 1, nj
            do i = 1, ni
               if (.not. IEEE_IS_NORMAL(f(i,j,k))) then
                  write (6,20) i,j,k
               end if
            end do
         end do
      end do

 20 format(' ','**** NaN at grid point(', i4.4,',',i4.4,',',i3.3,')')

      return
      end
