      subroutine RPN_COMM_xch_halo_8(g,minx,maxx,miny,maxy,ni,nj,nk,halox,haloy,periodx,periody,gni,npol_row)
      implicit none
!
!     exchange a halo with N/S/E/W neighbours for 64 bit items
!     use RPN_COMM_xch_halo (32 bit items) with fudged dimensions along x
!
      integer minx,maxx,miny,maxy,ni,nj,nk,halox,haloy
      integer gni,npol_row
      logical periodx,periody
!
      integer g(2*minx-1:2*maxx,miny:maxy,nk)
!     integer*8  g(minx:maxx,miny:maxy,nk)
!
      if(npol_row /=0) then
        print *,"ERROR: (RPN_COMM_xch_halo_8)",   &
             "npol_row must be zero when calling RPN_COMM_xch_halo_8"
        return  ! must figure out something sensible to signal errors
      endif
!
      call RPN_COMM_xch_halo(g,2*minx-1,2*maxx,miny,maxy,  &
                   2*ni,nj,nk,2*halox,haloy,periodx,periody,  &
                   2*gni,0)
      return
      end subroutine RPN_COMM_xch_halo_8