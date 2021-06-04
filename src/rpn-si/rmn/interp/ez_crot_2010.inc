  subroutine ez_crot_2010( r, ri, Grd_xlon1, Grd_xlat1, Grd_xlon2, Grd_xlat2 )
  implicit none

  integer i,j
  real r(3,3), ri(3,3)
  real Grd_xlat1, Grd_xlon1,Grd_xlon2, Grd_xlat2
     
  real*8 Grd_rot_8(3,3),xyz1_8(3),xyz2_8(3),xyz3_8(3),b_8
  real   xyz1(3), xyz2(3)

!   Grd_xlat1 = 45.
!   Grd_xlon1 = 0.
!   Grd_xlat2 = 0.5
!   Grd_xlon2 = 90.5

!  Grd_xlon1 = Grd_xlon1 + 180.0
!  Grd_xlon2 = Grd_xlon2 + 180.0

  call ez_crot(r, ri, Grd_xlon1, Grd_xlat1, Grd_xlon2, Grd_xlat2)

  print *, 'r et ri original'
  print *, '-- r  ----------'
  print *, r
  print *, '-- ri ----------'
  print *, ri
  print *, '################'

  call ez_lac_8 ( xyz1_8, Grd_xlon1, Grd_xlat1, 1 )
  call ez_lac_8 ( xyz2_8, Grd_xlon2, Grd_xlat2, 1 )           

  xyz3_8=xyz2_8

!    Compute r_3= r_1Xr_2 normalized (New North Pole)

  xyz3_8(1)=xyz1_8(2)*xyz2_8(3)-xyz1_8(3)*xyz2_8(2)
  xyz3_8(2)=xyz1_8(3)*xyz2_8(1)-xyz1_8(1)*xyz2_8(3)
  xyz3_8(3)=xyz1_8(1)*xyz2_8(2)-xyz1_8(2)*xyz2_8(1)

! Normalize

  b_8=sqrt(xyz3_8(1)**2 +xyz3_8(2)**2+xyz3_8(3)**2)
  xyz3_8(1)=xyz3_8(1)/b_8
  xyz3_8(2)=xyz3_8(2)/b_8
  xyz3_8(3)=xyz3_8(3)/b_8

!    Compute r_2= r_3Xr_1

  xyz2_8(1)=xyz3_8(2)*xyz1_8(3)-xyz3_8(3)*xyz1_8(2)
  xyz2_8(2)=xyz3_8(3)*xyz1_8(1)-xyz3_8(1)*xyz1_8(3)
  xyz2_8(3)=xyz3_8(1)*xyz1_8(2)-xyz3_8(2)*xyz1_8(1)

!  Now r_1, r_2 and r_3 are the rotated axes

  Grd_rot_8(1,1)=xyz1_8(1)
  Grd_rot_8(1,2)=xyz1_8(2)
  Grd_rot_8(1,3)=xyz1_8(3)

  Grd_rot_8(2,1)=xyz2_8(1)
  Grd_rot_8(2,2)=xyz2_8(2)
  Grd_rot_8(2,3)=xyz2_8(3)

  Grd_rot_8(3,1)=xyz3_8(1)
  Grd_rot_8(3,2)=xyz3_8(2)
  Grd_rot_8(3,3)=xyz3_8(3)

  r = real(Grd_rot_8)

  print *, Grd_xlon1, Grd_xlat1, Grd_xlon2, Grd_xlat2
  print *, Grd_rot_8
  print *, '--------------------'
  print *, r

  do j=1,3
    do i=1,3
      ri(i,j) = r(j,i)
    enddo
  enddo

  print *, '--------------------'
  print *, ri

  return
  end subroutine ez_crot_2010
