!    program testlat
!    real lons(2,2), lats(2,2)
!    real uv, lat, lon
!
!    integer pt_in_quad, res
!
!    lons(1,1) = 285.0
!    lons(2,1) = 288.0
!    lons(2,2) = 287.0
!    lons(1,2) = 286.0
!
!    lats(1,1) = 45.0
!    lats(2,1) = 43.0
!    lats(2,2) = 48.0
!    lats(1,2) = 47.0
!
!    lat = 45.5
!    lon = 285.5
!
!    res = pt_in_quad(lon,lat,lons(1,1),lats(1,1),lons(2,1), lats(2,1), &
!                      lons(1,2),lats(1,2),lons(2,2), lats(2,2))
!
!
!    call ez_uvfllc2d(u,v,lat, lon, lats, lons)
!    print *, u,v,lat, lon,lats, lons
!    print *
!    print *,'res = ', res
!
!    lat = lats(1,1)
!    lon = lons(1,1)
!    call ez_uvfllc2d(u,v,lat, lon, lats, lons)
!    print *, u,v,lat, lon,lats, lons
!    res = pt_in_quad(lon,lat,lons(1,1),lats(1,1),lons(2,1), lats(2,1), &
!                      lons(1,2),lats(1,2),lons(2,2), lats(2,2))
!    print *,'res = ', res
!    print *
!
!    lat = lats(2,1)
!    lon = lons(2,1)
!    call ez_uvfllc2d(u,v,lat, lon, lats, lons)
!    print *, u,v,lat, lon,lats, lons
!    res = pt_in_quad(lon,lat,lons(1,1),lats(1,1),lons(2,1), lats(2,1), &
!                      lons(1,2),lats(1,2),lons(2,2), lats(2,2))
!    print *,'res = ', res
!    print *
!
!    lat = lats(1,2)
!    lon = lons(1,2)
!    call ez_uvfllc2d(u,v,lat, lon, lats, lons)
!    print *, u,v,lat, lon,lats, lons
!    res = pt_in_quad(lon,lat,lons(1,1),lats(1,1),lons(2,1), lats(2,1), &
!                      lons(1,2),lats(1,2),lons(2,2), lats(2,2))
!    print *,'res = ', res
!    print *
!
!    lat = lats(2,2)
!    lon = lons(2,2)
!    call ez_uvfllc2d(u,v,lat, lon, lats, lons)
!    print *, u,v,lat, lon,lats, lons
!    res = pt_in_quad(lon,lat,lons(1,1),lats(1,1),lons(2,1), lats(2,1), &
!                      lons(1,2),lats(1,2),lons(2,2), lats(2,2))
!    print *,'res = ', res
!    print *
!
!    stop
!    end program testlat

   subroutine ez_uvfllc2d(u, v, x, y, x0, y0, x1, y1, x2, y2, x3, y3)
   implicit none
   real :: u, v, lat, lon, x, y

   real a,b,c,d,e, f, g,h, i
   real dx1, dx2, dy1, dy2, som_x ,som_y
   real x0, y0, x1, y1, x2, y2, x3, y3, q

   real, dimension(3, 3) ::  invmat

   som_x = x0 - x1 + x2 - x3
   som_y = y0 - y1 + y2 - y3

   dx1 = x1 - x2
   dx2 = x3 - x2
   dy1 = y1 - y2
   dy2 = y3 - y2

   g = (dy2*som_x-som_y*dx2)/(dx1*dy2-dx2*dy1)
   h = (dx1*som_y-som_x*dy1)/(dx1*dy2-dx2*dy1)
   i = 1.0

   a = x1 - x0 + g * x1
   b = x3 - x0 + h * x3
   c = x0
   d = y1 - y0 + g * y1
   e = y3 - y0 + h * y3
   f = y0

   invmat(1,3) = e * i - f * h
   invmat(2,3) = c * h - b * i
   invmat(3,3) = b * f - c * e
   invmat(1,2) = f * g - d * i
   invmat(2,2) = a * i - c * g
   invmat(3,2) = c * d - a * f
   invmat(1,1) = d * h - e * g
   invmat(2,1) = b * g - a * h
   invmat(3,1) = a * e - b * d

   q = invmat(1,1) * x + invmat(2,1) * y + invmat(3,1)
   if (q /= 0.0) then
      u = (invmat(1,3) * x + invmat(2,3) * y + invmat(3,3)) / q
      v = (invmat(1,2) * x + invmat(2,2) * y + invmat(3,2)) / q
      if (abs(u) < 0.01) u = abs(u)
      if (abs(v) < 0.01) v = abs(u)
   else
      u = -1.0
      v = -1.0
   endif

   end subroutine ez_uvfllc2d