   subroutine inside_or_outside(masque,x, y, out_lat,out_lon, gdin_lat, gdin_lon, ni_src, nj_src, wgts, idxs, num_wgts)
   implicit none

   integer :: masque, ni_src, nj_src, num_wgts
   real out_lat, out_lon, gdin_lat(ni_src, nj_src), gdin_lon(ni_src, nj_src)
   real :: wgts(num_wgts), x, y
   integer :: idxs(num_wgts,2), i,j, ix
   logical pt_in_quad
   external pt_in_quad

   ix = minloc(wgts, 1)
   i = min(max(2,idxs(ix,1)),ni_src-1)
   j = min(max(2,idxs(ix,2)),nj_src-1)

   !  Cas 1 - Coin superieur droit

   if (pt_in_quad(out_lon, out_lat, gdin_lon(i-1,j-1), gdin_lat(i-1,j-1), &
         gdin_lon(i,j-1), gdin_lat(i,j-1), gdin_lon(i,j), gdin_lat(i,j), &
         gdin_lon(i-1,j), gdin_lat(i-1,j))) then
         masque = 1
         call ez_uvfllc2d(x, y,out_lon, out_lat, &
            gdin_lon(i-1,j-1), gdin_lat(i-1,j-1), &
            gdin_lon(i,j-1), gdin_lat(i,j-1), &
            gdin_lon(i,j), gdin_lat(i,j), &
            gdin_lon(i-1,j), gdin_lat(i-1,j))
         x = x + (i-1)
         y = y + (j-1)

   !  Cas 2 - Coin superieur gauche

   else if (pt_in_quad(out_lon, out_lat, gdin_lon(i,j-1), gdin_lat(i,j-1), &
         gdin_lon(i+1,j-1), gdin_lat(i+1,j-1), gdin_lon(i+1,j), gdin_lat(i+1,j), &
         gdin_lon(i,j), gdin_lat(i,j))) then
         masque = 1
         call ez_uvfllc2d(x, y,out_lon, out_lat, &
            gdin_lon(i,j-1), gdin_lat(i,j-1), &
            gdin_lon(i+1,j-1), gdin_lat(i+1,j-1), &
            gdin_lon(i+1,j), gdin_lat(i+1,j), &
            gdin_lon(i,j), gdin_lat(i,j))
         x = x + (i)
         y = y + (j-1)

!  Cas 3 - Coin inferieur droit

   else if (pt_in_quad(out_lon, out_lat, gdin_lon(i-1,j), gdin_lat(i-1,j), &
         gdin_lon(i,j), gdin_lat(i,j), gdin_lon(i,j+1), gdin_lat(i,j+1), &
         gdin_lon(i-1,j+1), gdin_lat(i-1,j+1))) then
         masque = 1
         call ez_uvfllc2d(x, y,out_lon, out_lat, &
            gdin_lon(i-1,j), gdin_lat(i-1,j), &
            gdin_lon(i,j), gdin_lat(i,j), &
            gdin_lon(i,j+1), gdin_lat(i,j+1), &
            gdin_lon(i-1,j+1), gdin_lat(i-1,j+1))
         x = x + (i-1)
         y = y + (j)

!  Cas 4 - Coin inferieur gauche

   else if (pt_in_quad(out_lon, out_lat, gdin_lon(i,j), gdin_lat(i,j), &
         gdin_lon(i+1,j), gdin_lat(i+1,j), gdin_lon(i+1,j+1), gdin_lat(i+1,j+1), &
         gdin_lon(i,j+1), gdin_lat(i,j+1))) then
         masque = 1
         call ez_uvfllc2d(x, y,out_lon, out_lat, &
            gdin_lon(i,j), gdin_lat(i,j), &
            gdin_lon(i+1,j), gdin_lat(i+1,j), &
            gdin_lon(i+1,j+1), gdin_lat(i+1,j+1), &
            gdin_lon(i,j+1), gdin_lat(i,j+1))
         x = x + (i)
         y = y + (j)
   else
      masque = 0
      x = -1.0
      y = -1.0
   endif

   return

   end subroutine inside_or_outside
