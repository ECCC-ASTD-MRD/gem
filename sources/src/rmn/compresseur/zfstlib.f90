! 
!=============================================================================
  subroutine ibicubic_int4(izo,ni,nj,step,ajus_x,ajus_y)
  implicit none
  integer npts,ni,nj,nic,njc,step,ajus_x,ajus_y
  integer izo(ni,nj)
  real*8 y1,y2,y3,y4
  real*8 one, two,three,six, third, sixth,half,my_half
  real*8 z,z11,z12,z13,z14,z21,z22,z23,z24,z31,z32,z33,z34,z41,z42,z43,z44,fac1,fac2,unsurfac2
  integer m,n,i,j,ii,jj,ic,jc,iref,jref,nimax,njmax,nilim,njlim
  parameter (one = 1.0D0)
  parameter (two = 2.0D0)
  parameter (three = 3.0D0)
  parameter (six = 6.0D0)
  parameter (sixth = one/six)
  parameter (third = one/three)
  parameter (half = 0.5D0)
  parameter (my_half = 0.5001D0)
  parameter (fac1 = 108.0D0) 
  parameter (fac2 = 1944.0D0) 
  parameter (unsurfac2 = 1.0D0/fac2) 
  real*8 icubic, icubic1, icubic2, cubic, dx,dy,z1,z2,z3,z4
  integer my_nint
!  cubic(z1,z2,z3,z4,dx)=((((z4-z1)*sixth + 0.5*(z2-z3))*dx  + 0.5*(z1+z3)-z2)*dx  + z3-sixth*z4-0.5*z2-third*z1)*dx+z2
  icubic(z1,z2,z3,z4,dx)=z2+(dx*(6*(dx*(2*(dx*((z4-z1)+3*(z2-z3)))+18*((z1+z3)-2*z2)))+fac1*(6*z3-z4-3*z2-2*z1)))*unsurfac2
!  icubic1(z1,z2,z3,z4)=z2+(1*(6*(1*(2*(1*((z4-z1)+3*(z2-z3)))+18*((z1+z3)-2*z2)))+fac1*(6*z3-z4-3*z2-2*z1)))/fac2
!  icubic2(z1,z2,z3,z4)=z2+(2*(6*(2*(2*(2*((z4-z1)+3*(z2-z3)))+18*((z1+z3)-2*z2)))+fac1*(6*z3-z4-3*z2-2*z1)))/fac2
  my_nint(z)=(z+sign(my_half,z))
!  nimax = min(ni-ajus_x,ni-2)
  if (ajus_x == 0) then
    nimax = ni - 3
    nilim = nimax - 3
  else if (ajus_x == 1) then
    nimax = ni - 4
    nilim = nimax
  else
    nimax = ni - 5
    nilim = nimax
  endif
  if (ajus_y == 0) then
    njmax = nj - 3
    njlim = njmax - 3
  else if (ajus_y == 1) then
    njmax = nj - 4
    njlim = njmax
  else
    njmax = nj - 5
    njlim = njmax
  endif
  do j=1,nj-ajus_y,step
     do i=1,nimax,step
      iref = min(nilim,max(4,i))
      z12 = dble(izo(iref-step,j ))
      z22 = dble(izo(iref,j ))
      z32 = dble(izo(iref+step,j  ))
      z42 = dble(izo(min(ni,iref+2*step),j  ))
!     m=0, n=0
!     m=1, n=0
      dx = dble(i+1-iref)
!      izo(i+1,j)=my_nint(cubic(z12,z22,z32,z42,dx)) 
      izo(i+1,j)=my_nint(icubic(z12,z22,z32,z42,dx))
!     m=2, n=0
      dx = dble(i+2-iref)
      izo(i+2,j)=my_nint(icubic(z12,z22,z32,z42,dx))
     enddo
  enddo
  if (ajus_x == 2) then
    do j=1,nj-ajus_y,step
      izo(ni-1,j)=my_nint(half*(dble(izo(ni,j))+dble(izo(ni-2,j))))
    enddo
  endif
  do j=0,ajus_y
     do i=1,nimax,step
      iref = min(nilim,max(4,i))
      z12 = dble(izo(iref-step,nj-j ))
      z22 = dble(izo(iref,nj-j ))
      z32 = dble(izo(iref+step,nj-j  ))
      z42 = dble(izo(min(ni,iref+2*step),nj-j))
!     m=0, n=0
!     m=1, n=0
      dx = dble(i+1-iref)
      izo(i+1,nj-j)=my_nint(icubic(z12,z22,z32,z42,dx))
!     m=2, n=0
      dx = dble(i+2-iref)
      izo(i+2,nj-j)=my_nint(icubic(z12,z22,z32,z42,dx))
     enddo
    if (ajus_x == 2) then
       izo(ni-1,nj-j)=my_nint(half*(dble(izo(ni,nj-j))+dble(izo(ni-2,nj-j))))
    endif
  enddo
  do j=1,njmax,step
     jref = min(njlim,max(4,j))
     do i=1,ni
      z21 = dble(izo(i,jref-step))
      z22 = dble(izo(i,jref  ))
      z23 = dble(izo(i,jref+step))
      z24 = dble(izo(i,min(nj,jref+2*step)))
!     m=0, n=1
      dy = dble(j+1-jref)
      izo(i,j+1)=my_nint(icubic(z21,z22,z23,z24,dy))
!     m=0, n=2
      dy = (dble(j+2-jref))
      izo(i,j+2)=my_nint(icubic(z21,z22,z23,z24,dy))
     enddo
  enddo
  if (ajus_y == 2) then
    do i=1,ni
      izo(i,nj-1)=my_nint(half*(dble(izo(i,nj))+dble(izo(i,nj-2))))
    enddo
  endif
  return
  end
!=============================================================================
  subroutine ibicubic_int3(izo,ni,nj,izc,nic,njc,step)
  implicit none
  integer npts,ni,nj,nic,njc,step
  integer izo(ni,nj)
  integer izc(nic,njc)
  real zc(nic,njc)
  real*8 y1,y2,y3,y4
  real*8 one, three,six, third, sixth
  real*8 z11,z12,z13,z14,z21,z22,z23,z24,z31,z32,z33,z34,z41,z42,z43,z44
  integer m,n,i,j,ii,jj,ic,jc
  parameter (one = 1.0D0)
  parameter (three = 3.0D0)
  parameter (six = 6.0D0)
  parameter (sixth = one/six)
  parameter (third = one/three)
  real*8 cubic, dx,dy,z1,z2,z3,z4
  cubic(z1,z2,z3,z4,dx)=((((z4-z1)*sixth + 0.5*(z2-z3))*dx  + 0.5*(z1+z3)-z2)*dx  + z3-sixth*z4-0.5*z2-third*z1)*dx+z2
  do j=1,njc
     do i=1,nic
      zc(i,j) = real(izc(i,j))
     enddo
  enddo
  do jc=1,njc-2
     do ic=1,nic-2
      i = min(nic-2,max(2,ic))
      j = min(njc-2,max(2,jc))
      z11 = dble(zc(i-1,j-1))
      z12 = dble(zc(i-1,j  ))
      z13 = dble(zc(i-1,j+1))
      z14 = dble(zc(i-1,j+2))
      z21 = dble(zc(i,j-1))
      z22 = dble(zc(i,j  ))
      z23 = dble(zc(i,j+1))
      z24 = dble(zc(i,j+2))
      z31 = dble(zc(i+1,j-1))
      z32 = dble(zc(i+1,j  ))
      z33 = dble(zc(i+1,j+1))
      z34 = dble(zc(i+1,j+2))
      z41 = dble(zc(i+2,j-1))
      z42 = dble(zc(i+2,j  ))
      z43 = dble(zc(i+2,j+1))
      z44 = dble(zc(i+2,j+2))
!     m=0, n=0
!       dy = 0.0
!       dx = 0.0
!       y1=cubic(z11,z21,z31,z41,dx)
!       y2=cubic(z12,z22,z32,z42,dx)
!       y3=cubic(z13,z23,z33,z43,dx)
!       y4=cubic(z14,z24,z34,z44,dx)
      izo(step*(ic-1)+1,step*(jc-1)+1)=nint(z22)
!     m=1, n=0
      dy = 0
      dx = third
      y2=cubic(z12,z22,z32,z42,dx)
      izo(step*(ic-1)+2,step*(jc-1)+1)=nint(y2)
!     m=2, n=0
      dx = 2.0 * third
      y2=cubic(z12,z22,z32,z42,dx)
      izo(step*(ic-1)+3,step*(jc-1)+1)=nint(y2)
!     m=0, n=1
      dy = third
      dx = 0.0
      izo(step*(ic-1)+1,step*(jc-1)+2)=nint(cubic(z21,z22,z23,z24,dy))
!     m=1, n=1
      dx = third
      y1=cubic(z11,z21,z31,z41,dx)
      y2=cubic(z12,z22,z32,z42,dx)
      y3=cubic(z13,z23,z33,z43,dx)
      y4=cubic(z14,z24,z34,z44,dx)
      izo(step*(ic-1)+2,step*(jc-1)+2)=nint(cubic(y1,y2,y3,y4,dy))
!     m=2, n=1
      dx = 2.0 * third
      y1=cubic(z11,z21,z31,z41,dx)
      y2=cubic(z12,z22,z32,z42,dx)
      y3=cubic(z13,z23,z33,z43,dx)
      y4=cubic(z14,z24,z34,z44,dx)
      izo(step*(ic-1)+3,step*(jc-1)+2)=nint(cubic(y1,y2,y3,y4,dy))
!     m=0, n=2
      dy = 2.0*third
      dx = 0.0
      izo(step*(ic-1)+1,step*(jc-1)+3)=nint(cubic(z21,z22,z23,z24,dy))
!     m=1, n=2
      dx = third
      y1=cubic(z11,z21,z31,z41,dx)
      y2=cubic(z12,z22,z32,z42,dx)
      y3=cubic(z13,z23,z33,z43,dx)
      y4=cubic(z14,z24,z34,z44,dx)
      izo(step*(ic-1)+2,step*(jc-1)+3)=nint(cubic(y1,y2,y3,y4,dy))
!     m=2, n=2
      dx = 2.0 * third
      y1=cubic(z11,z21,z31,z41,dx)
      y2=cubic(z12,z22,z32,z42,dx)
      y3=cubic(z13,z23,z33,z43,dx)
      y4=cubic(z14,z24,z34,z44,dx)
      izo(step*(ic-1)+3,step*(jc-1)+3)=nint(cubic(y1,y2,y3,y4,dy))
     enddo
  enddo
  return
  end
!=============================================================================
subroutine fill_coarse_grid(zc, nicoarse, njcoarse, z, ni, nj, istep)
  implicit none
  integer nicoarse, njcoarse, ni, nj, istep
  integer zc(nicoarse, njcoarse), z(ni, nj)
  integer i, j
  if (ni > 1 .and. nj > 1) then
     do j=1,njcoarse-1
        do i=1,nicoarse-1
           zc(i,j) = z(istep*(i-1)+1,istep*(j-1)+1)
        enddo
     enddo
     do j=1,njcoarse-1
        zc(nicoarse,j) = z(ni,istep*(j-1)+1)
     enddo
     do i=1,nicoarse-1
        zc(i,njcoarse) = z(istep*(i-1)+1,nj)
     enddo
     zc(nicoarse, njcoarse) = z(ni,nj)
     return
  endif
  if (nj == 1) then
     do i=1,nicoarse-1
        zc(i,1) = z(istep*(i-1)+1,1)
     enddo
     zc(nicoarse,1) = z(ni,1)
  endif
  if (ni == 1) then
     do j=1,njcoarse-1
        zc(1,j) = z(1,istep*(j-1)+1)
     enddo
     zc(1,njcoarse) = z(1,nj)
  endif
end subroutine fill_coarse_grid
subroutine fill_coarse_nodes(z, ni, nj, zc, nicoarse, njcoarse, istep)
  implicit none
  integer nicoarse, njcoarse, ni, nj, istep
  integer zc(nicoarse, njcoarse), z(ni, nj)
  integer i, j
  if (ni > 1 .and. nj > 1) then
     do j=1,njcoarse-1
        do i=1,nicoarse-1
           z(istep*(i-1)+1,istep*(j-1)+1) = zc(i,j)
        enddo
     enddo
     do j=1,njcoarse-1
        z(ni,istep*(j-1)+1) = zc(nicoarse,j)
     enddo
     do i=1,nicoarse-1
        z(istep*(i-1)+1,nj) = zc(i,njcoarse)
     enddo
     z(ni, nj) = zc(nicoarse,njcoarse)
     return
  endif
end subroutine fill_coarse_nodes
subroutine fill_last_colrows(px, py, z, ni, nj, nicoarse, njcoarse, istep)
  implicit none
  integer ni,nj,nicoarse,njcoarse, istep, z(ni,nj)
  real px(ni,nj), py(ni,nj) , rstep, dx, dy
  integer i,j, istart, jstart,nintervalles
  rstep = 1.0 * istep
  istart = (nicoarse-1)*istep
  jstart = (njcoarse-1)*istep
!  print *, istart, jstart
  if (jstart .ne. nj) then
     nintervalles = nj - jstart
     dy = 1.0 / nintervalles
     do j=jstart, nj
        do i=1,ni
           px(i,j) = 1.0+1.0*(i-1)/rstep
           py(i,j) = 1.0*(njcoarse-1) + (j-jstart)*dy
!           print *, i, j, px(i,j), py(i,j)
       enddo
     enddo
  endif
!  print *, 'rangee du haut completee...'
  if (istart .ne. ni) then
     nintervalles = ni - istart
     dx = 1.0 / nintervalles
     do i=istart, ni
        do j=1,nj
           px(i,j) = 1.0*(nicoarse-1) + (i-istart)*dx
           py(i,j) = 1.0+1.0*(j-1)/rstep
!           print *, i, j, px(i,j), py(i,j)
        enddo
     enddo
  endif
!  print *, 'colonne de droite completee...'
  do j=jstart,nj
     do i=istart,ni
        px(i,j) = 1.0*(nicoarse-1) + (i-istart)*dx
        py(i,j) = 1.0*(njcoarse-1) + (j-jstart)*dy
!           print *, i, j, px(i,j), py(i,j)
     enddo
  enddo
!  print *, 'coin complete...'
  return	
end subroutine fill_last_colrows
