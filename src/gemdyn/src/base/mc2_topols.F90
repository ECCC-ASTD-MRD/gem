      subroutine mc2_topols (F_dst,F_src,minx,maxx,miny,maxy,F_np)
      use glb_ld

      implicit none

      integer, intent(IN) :: minx,maxx,miny,maxy, F_np
      real, dimension(minx:maxx,miny:maxy), intent(IN ) :: F_src
      real, dimension(minx:maxx,miny:maxy), intent(OUT) :: F_dst

      integer n,nflt
      real, dimension(minx:maxx,miny:maxy) :: w1,w2
!
!-----------------------------------------------------------------------
!
      w1= F_src ; nflt= F_np
      
      do n=1,nflt/2
         call gem_xch_halo (w1, l_minx, l_maxx, l_miny, l_maxy, 1)
         call dc_topop (w2,w1, minx,maxx,miny,maxy,l_ni,l_nj)
         call gem_xch_halo (w2, l_minx, l_maxx, l_miny, l_maxy, 1)
         call dc_topop (w1,w2, minx,maxx,miny,maxy,l_ni,l_nj)
      end do
      
      call gem_xch_halo (w1, l_minx, l_maxx, l_miny, l_maxy, 1)
      F_dst= w1
!
!-----------------------------------------------------------------------
!
      return
      end
      
      subroutine dc_topop (d,s,lminx,lmaxx,lminy,lmaxy,ldni,ldnj)
      use glb_ld
      implicit none
      integer, intent(IN) :: lminx,lmaxx,lminy,lmaxy,ldni,ldnj
      real, intent(IN ) :: s(lminx:lmaxx,lminy:lmaxy)
      real, intent(OUT) :: d(lminx:lmaxx,lminy:lmaxy)

      integer i,j,hx,hy
!
!-----------------------------------------------------------------------
!
      hx= 1-lminx ; hy= hx

      do j=2-hy,ldnj+hy-1
      do i=2-hx,ldni+hx-1
         d(i,j) = s(i,j)/2.+(s(i-1,j)+s(i+1,j)+s(i,j-1)+s(i,j+1))/8.
      end do
      end do
      
      if (l_south.and.l_west)  &  
      d(1-hx,1-hy)   = s(1-hx,1-hy)/2.  + s(2-hx,1-hy)/4.  + s(1-hx,2-hy)/4.
      
      if (l_north.and.l_west) &  
      d(1-hx,ldnj+hy)  = s(1-hx,ldnj+hy)/2.  + s(1-hx,ldnj+hy-1)/4. + s(2-hx,ldnj+hy)/4.
      
      if (l_south.and.l_east)  & 
      d(ldni+hx,1-hy)  = s(ldni+hx,1-hy)/2.  + s(ldni+hx-1,1-hy)/4.  + s(ldni+hx,2-hy)/4.
      
      if (l_north.and.l_east)  & 
      d(ldni+hx,ldnj+hy) = s(ldni+hx,ldnj+hy)/2. + s(ldni+hx-1,ldnj+hy)/4. + s(ldni+hx,ldnj+hy-1)/4.
      
      if (l_west) then
         do j=2-hy,ldnj+hy-1
            d(1-hx,j) = s(1-hx,j)/2.  + s(2-hx,j)/4. + (s(1-hx,j-1)  + s(1-hx,j+1))/8.
         end do
      endif
      
      if (l_east) then
         do j=2-hy,ldnj+hy-1
            d(ldni+hx,j)= s(ldni+hx,j)/2.+s(ldni+hx-1,j)/4.+(s(ldni+hx,j-1)+s(ldni+hx,j+1))/8.
         end do 
      endif
      
      if (l_south) then
         do i=2-hx,ldni+hx-1
            d(i,1-hy)  = s(i,1-hy)/2.  + s(i,2-hy)/4. + (s(i-1,1-hy)  + s(i+1,1-hy))/8.
         end do
      endif
      
      if (l_north) then
         do i=2-hx,ldni+hx-1
            d(i,ldnj+hy) = s(i,ldnj+hy)/2.+s(i,ldnj+hy-1)/4.+(s(i-1,ldnj+hy)+s(i+1,ldnj+hy))/8.
          end do
      endif
!
!-----------------------------------------------------------------------
!
      return
      end
