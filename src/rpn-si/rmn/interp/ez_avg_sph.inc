   subroutine ez_avg_sph(zout, xx, yy, lats_dst, ni_dst, nj_dst, zin, ni_src, nj_src,extension)
   implicit none

   integer udst, extension
   integer ni_src, nj_src, ni_dst, nj_dst
   real, dimension(ni_src, nj_src) :: zin
   real, dimension(ni_dst, nj_dst) :: zout,xx,yy
   real, dimension(nj_dst) :: lats_dst
   integer avg

   real, parameter :: degre_a_radian = 0.017453295199
   real, dimension(:), allocatable :: x, y, y_low, y_high, amplif_lats
   real, dimension(:,:), allocatable :: ztmp, x_low, x_high

   integer i, j, nix, njx, nkx, imid, jmid
   integer lcl_avg

   real dxcore, dycore, dxcoarse, dycoarse
   real total_area, xtile_min, ytile_min, xtile_max, ytile_max, xfrac, yfrac, area
   integer ez_cherche,ii,jj
   integer istart, iend, jstart, jend, compression_code, usr_datyp
   integer low_bound_x, high_bound_x, low_bound_y, high_bound_y

   if (extension == 0) then
      low_bound_x = 1
      high_bound_x = ni_src
      allocate(ztmp(low_bound_x:high_bound_x, nj_src))
      do j=1,nj_src
         do i=1,ni_src
            ztmp(i,j) = zin(i,j)
         enddo
      enddo
   else if (extension == 1) then
      low_bound_x =  -ni_src
      high_bound_x = 2 * ni_src - 2
      allocate(ztmp(low_bound_x:high_bound_x, nj_src))
      do j=1,nj_src
         do i=1,ni_src-1
            ztmp(i,j) = zin(i,j)
            ztmp(1-ni_src+i,j) = zin(i,j)
            ztmp(ni_src-1+i,j) = zin(i,j)
         enddo
      enddo
   else if (extension == 2) then
      low_bound_x = 1 - ni_src
      high_bound_x = 2 * ni_src
      allocate(ztmp(low_bound_x:high_bound_x, nj_src))
      do j=1,nj_src
         do i=1,ni_src
            ztmp(i,j) = zin(i,j)
            ztmp(-ni_src+i,j) = zin(i,j)
            ztmp(ni_src+i,j) = zin(i,j)
         enddo
      enddo
   else
      print *, 'Invalid extension'
   endif


   allocate(x_low(ni_dst,nj_dst))
   allocate(y_low(nj_dst))
   allocate(x_high(ni_dst,nj_dst))
   allocate(y_high(nj_dst))
   allocate(x(ni_dst),y(nj_dst),amplif_lats(nj_dst))

   do j=1,nj_dst
      amplif_lats(j) = 1.0/cos(lats_dst(j)*degre_a_radian)
   enddo

   x = xx(:,1)
   y = yy(1,:)
   if (x(1) > real(ni_src-1)) x(1) = 1.0

   do j=1,nj_dst
      x_low(1,j) = x(1)-0.5*(x(2)-x(1))*amplif_lats(j)
      do i=2,ni_dst
         x_low(i,j) = x(i)-0.5*(x(i)-x(i-1))*amplif_lats(j)
      enddo
   enddo

   y_low(1) = max(1.0,y(1)-0.5*(y(2)-y(1)))
   do j=2,nj_dst
      y_low(j) =y(j)- 0.5*(y(j)-y(j-1))
   enddo

   do j=1,nj_dst
      x_high(ni_dst,j) = x(ni_dst)+0.5*(x(ni_dst)-x(ni_dst-1))*amplif_lats(j)
      do i=1,ni_dst-1
         x_high(i,j) = x(i)+0.5 * (x(i+1)-x(i))*amplif_lats(j)
      enddo
   enddo

   y_high(nj_dst) = min(1.0*nj_src,y(nj_dst)+0.5*(y(nj_dst)-y(nj_dst-1)))
   do j=1,nj_dst-1
      y_high(j) = y(j)+0.5 * (y(j+1)-y(j))
   enddo

   do j=2,nj_dst-1
      jstart = int(y_low(j))
      jend   = nint(y_high(j))
      if ((0.5+real(jstart)) < y_low(j)) jstart = jstart + 1
!      if (real(jend) < y_high(j)) jend = jend + 1
      do i=1,ni_dst
         zout(i,j) = 0.0
         total_area = 0.0
         istart = int(x_low(i,j))
         iend   = nint(x_high(i,j))
         if ((0.5+real(istart)) < x_low(i,j)) istart = istart + 1
!         if (real(iend) < x_high(j)) iend = iend + 1
         do jj=jstart, jend
            ytile_min = real(jj)-0.5
            ytile_max = real(jj)+0.5
            yfrac = 1.0
            if (ytile_min < y_low(j)) then
               yfrac = ytile_max - y_low(j)
            endif
            if (ytile_max > y_high(j)) then
               yfrac = y_high(j) - ytile_min
            endif
            do ii=istart, iend
               xtile_min = real(ii)-0.5
               xtile_max = real(ii)+0.5
               xfrac = 1.0
               if (xtile_min < x_low(i,j)) then
                  xfrac = xtile_max - x_low(i,j)
               endif
               if (xtile_max > x_high(i,j)) then
                  xfrac = x_high(i,j) - xtile_min
               endif
               area =  xfrac*yfrac
               total_area = total_area + area
               zout(i,j) = zout(i,j) + ztmp(ii,jj) * area
            enddo
         enddo
         if (total_area /= 0.0) zout(i,j) = zout(i,j)/total_area
      enddo
   enddo

! Moyenne 1e rangee

  j = 1
  jstart = int(y_low(1))
  jend = nint(y_high(1))
  if (real(jstart) > y_low(1)) jstart = jstart - 1
  do i=1,ni_dst
      zout(i,j) = 0.0
      total_area = 0.0
      istart = int(x_low(i,j))
      iend   = nint(x_high(i,j))
      if ((0.5+real(istart)) < x_low(i,j)) istart = istart + 1
      if (real(iend) < x_high(i,j)) iend = iend + 1
      do jj=jstart, jend

         if (jj == 1) then
           ytile_min = 1.0
         else
           ytile_min = real(jj)-0.5
         endif

         ytile_max = real(jj)+0.5
         yfrac = ytile_max - ytile_min
         if (ytile_min < y_low(j)) then
            yfrac = ytile_max - y_low(j)
         endif
         if (ytile_max > y_high(j)) then
            yfrac = y_high(j) - ytile_min
         endif

         do ii=istart, iend
            xtile_min = real(ii)-0.5
            xtile_max = real(ii)+0.5
            xfrac = 1.0
            if (xtile_min < x_low(i,j)) then
               xfrac = xtile_max - x_low(i,j)
            endif
            if (xtile_max > x_high(i,j)) then
               xfrac = x_high(i,j) - xtile_min
            endif
            area =  xfrac*yfrac
            total_area = total_area + area
            zout(i,j) = zout(i,j) + ztmp(ii,jj) * area
         enddo
      enddo
         if (total_area /= 0.0) zout(i,j) = zout(i,j)/total_area
   enddo


! Moyenne rangee du haut

   j = nj_dst
   jstart = int(y(nj_dst))
   jend = nint(y_high(nj_dst))
   do i=1,ni_dst
      zout(i,j) = 0.0
      total_area = 0.0
      istart = int(x_low(i,j))
      if ((0.5+real(istart)) < x_low(i,j)) istart = istart + 1
      iend   = nint(x_high(i,j))
      do jj=jstart, jend
         if (jj == nj_src) then
            ytile_max = real(nj_src)
         else
            ytile_max = real(jj)+0.5
         endif
         ytile_min = real(jj)-0.5
         yfrac = ytile_max - ytile_min
         if (ytile_min < y_low(j)) then
            yfrac = ytile_max - y_low(j)
         endif
         if (ytile_max > y_high(j)) then
            yfrac = y_high(j) - ytile_min
         endif

         do ii=istart, iend
            xtile_min = real(ii)-0.5
            xtile_max = real(ii)+0.5
            xfrac = 1.0
            if (xtile_min < x_low(i,j)) then
               xfrac = xtile_max - x_low(i,j)
            endif
            if (xtile_max > x_high(i,j)) then
               xfrac = x_high(i,j) - xtile_min
            endif
            area =  xfrac*yfrac
            total_area = total_area + area
            zout(i,j) = zout(i,j) + ztmp(ii,jj) * area
         enddo
      enddo
         if (total_area /= 0.0) zout(i,j) = zout(i,j)/total_area
   enddo

   deallocate(x_low, y_low, x_high, y_high, x, y, amplif_lats, ztmp)

   return
  end