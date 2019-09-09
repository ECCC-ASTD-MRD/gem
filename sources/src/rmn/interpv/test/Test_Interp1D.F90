!!s/r Test_Interp1D - test routines in the 1D interpolation package
program Test_Interp1D
!
!AUTHOR
!     J.W. Blezius JULY 2002
!
!REVISION
! v1_0    Blezius J.W.          - new code
! v1_1    Blezius J.W. SEP 2002 - make all input dim's the same; ditto for output
! v1_2    Blezius J.W. OCT 2002 - add test for extrapolation in the case where
!                                 the vertical levels are in DEcreasing order
! v1_3    Blezius J.W. OCT 2003 - add the extended interface
!
!OBJECT
!        To verify the 1D interpolation package.
!
!ARGUMENTS
  implicit none
!
!NOTES
!        Although, for those functions that do not require the extended
!        interface, this test calls only the non-extended interface, the extended
!        interface has also been implicitly tested by virtue of the fact that the
!        non-extended shell calls the extended interface itself.  Nonetheless,
!        just to ensure that the extended interface is visible, one such function
!        is called.
!
!!

! number of points, horizontally
  integer :: hNumPts=2
# define hDim    6

! number of points, vertically, in target
# define nPtsTgt 5

! number of points, vertically, in source
# define nPtsSrc 4


! variables used for testing Extrap1D_Surface_X
  real ft(hDim),ilmo(hDim),hBound(hDim),z0(hDim),f(hDim)
  real hmin,factn
  parameter ( hmin = 40., factn =1.2)


  external FluxStub, slfun_uvStub
  integer :: numExtArraysIn
  integer :: numExtArraysOut
  real, dimension(hDim, 7) :: ExtArraysIn
  real, dimension(hDim, nPtsTgt) :: ExtArraysOut

  logical :: pass = .true., &
             unitPass
  integer :: testCount = 0
  integer :: i, v

  ! definition of the horizontal grid
  real, parameter, dimension(hDim) :: pSurf=(/100000., 100000., 100000., 100000., 100000., 100000./)

  ! definition of the initial grid:  the actual vertical levels
  real, dimension(nPtsSrc) :: levelsSrc=(/0.64, 1.25, 2.44, 2.97/)
  real, dimension(hDim, nPtsSrc) :: levelsIn ! filled-in cube

  ! definition of the target grid:  the actual vertical levels
  real, dimension(nPtsTgt) :: levelsTgt=(/1.13, 2.62, 2.79, 0.5, 3.1/)
  real, dimension(hDim, nPtsTgt) :: levelsOut ! filled-in cube

                                        ! indices to the levelsSrc which are just
                                        ! below the levelsTgt
  integer, dimension(hDim, nPtsTgt) :: interpIndex


  ! input state and derivative
  real, dimension(hDim, nPtsSrc) :: stateIn, &
                                       derivIn, &
                                       stateAscending, stateDescending, &
                                       derivAscending, derivDescending

  ! output state and derivative
  real, dimension(hDim, nPtsTgt) :: stateOut, &
                                       derivOut

  ! extra components for the wind
  real, dimension(hDim) :: v_windIn, latitude, dudz, dvdz, tha
  real, dimension(hDim, nPtsTgt) :: v_windOut
  real :: thmax

  real :: lapseRateDown = 0.4, &
          lapseRateUp = -0.5


! descending
  stateDescending(1,4) = 0.5971954   ! sin 0.64
  stateDescending(1,3) = 0.9489846   ! sin 1.25
  stateDescending(1,2) = 0.6454      ! sin 2.44
  stateDescending(1,1) = 0.1708      ! sin 2.97

  derivDescending(1,4) = 0.8020958   ! cos 0.64
  derivDescending(1,3) = 0.3153224   ! cos 1.25
  derivDescending(1,2) =-0.7638      ! cos 2.44
  derivDescending(1,1) =-0.9853      ! cos 2.97


  stateDescending(2,4) = 0.744544    ! tan 0.64
  stateDescending(2,3) = 3.009570    ! tan 1.25
  stateDescending(2,2) =-0.8450      ! tan 2.44
  stateDescending(2,1) =-0.1733      ! tan 2.97

  derivDescending(2,4) = 1.554346    ! 1/cos2 0.64
  derivDescending(2,3) = 10.057510   ! 1/cos2 1.25
  derivDescending(2,2) = 1.7140      ! 1/cos2 2.44
  derivDescending(2,1) = 1.03003     ! 1/cos2 2.97


! ascending
  do v=1,4
    do i=1,2
      stateAscending(i,v) = stateDescending(i,nPtsSrc+1-v)
      derivAscending(i,v) = derivDescending(i,nPtsSrc+1-v)
    end do
  end do

  
  ! Create cubes of vertical levels
  do v=1,nPtsSrc
    do i=1,hNumPts
      levelsIn(i,v) = levelsSrc(nPtsSrc+1-v) ! DESCENDING ORDER
    end do
  end do

  do v=1,nPtsTgt
    do i=1,hNumPts
      levelsOut(i,v) = levelsTgt(v)
    end do
  end do

  stateIn = stateDescending
  derivIn = derivDescending
!!$

  !
  ! SOURCE LEVELS IN DESCENDING ORDER
  !
  call Interp1D_FindPos(hNumPts, nPtsSrc, nPtsTgt, &
                        hDim, hDim, &

                        levelsIn, interpIndex, levelsOut &
                       )



  call Interp1D_CubicWithDerivs (hNumPts, nPtsSrc, nPtsTgt, &
                                 hDim, hDim, &

                                 levelsIn, stateIn, derivIn, &

                                 interpIndex, levelsOut, stateOut, derivOut, &

                                 .false., .false., &
                                 lapseRateDown, lapseRateUp &
                                )

  call TestReportClear(levelsSrc, stateAscending, derivAscending, levelsTgt, &
                      stateOut, derivOut, hnumpts, &
                      'Cubic with Derivatives - DEscending source levels', pass,&
                      testCount)








  call Extrap1D_LapseRate (hNumPts, nPtsSrc, nPtsTgt, &
                                  hDim, hDim, &

                                  levelsIn, stateIn, derivIn, &

                                  interpIndex, levelsOut, stateOut, derivOut, &

                                  .true., .true., &
                                  lapseRateDown, lapseRateUp &
                                 )

  ! Test this case individually
  unitPass = .true.
  if(     abs(stateOut(1,4) -  0.5411954)  > 1e-7 &
     .or. abs(stateOut(1,5) -  0.1058001)  > 1e-7 &
     .or. abs(stateOut(2,4) -  0.6885440)  > 1e-7 &
     .or. abs(stateOut(2,5) -(-0.2382999)) > 1e-7 &
    ) then
    unitPass = .false.
    pass = .false.
  endif

  call ReportAndClear(unitPass, stateOut, derivOut, hnumpts, &
                      'extrapolation by lapse rate - DEscending source levels', &
                      testCount)







  !
  ! SOURCE LEVELS IN ASCENDING ORDER
  !

  ! Change cube of source vertical levels to ASCENDING order
  do v=1,nPtsSrc
    do i=1,hNumPts
      levelsIn(i,v) = levelsSrc(v)
    end do
  end do

  stateIn = stateAscending
  derivIn = derivAscending

  call Interp1D_FindPos(hNumPts, nPtsSrc, nPtsTgt, &
                        hDim, hDim, &

                        levelsIn, interpIndex, levelsOut &
                       )



  call Interp1D_CubicWithDerivs (hNumPts, nPtsSrc, nPtsTgt, &
                                 hDim, hDim, &

                                 levelsIn, stateIn, derivIn, &

                                 interpIndex, levelsOut, stateOut, derivOut, &

                                 .false., .false., &
                                 lapseRateDown, lapseRateUp &
                                )

  call TestReportClear(levelsSrc, stateIn, derivIn, levelsTgt, stateOut, &
                       derivOut, hnumpts, &
                       'Cubic with Derivatives - AScending source levels', pass,&
                       testCount)







  call Interp1D_NearestNeighbour (hNumPts, nPtsSrc, nPtsTgt, &
                                  hDim, hDim, &

                                  levelsIn, stateIn, derivIn, &

                                  interpIndex, levelsOut, stateOut, derivOut, &

                                  .false., .false., &
                                  lapseRateDown, lapseRateUp &
                                 )

  ! Test this case individually
  unitPass = .true.
  do i=1,hNumPts
    if(     abs(stateOut(i,1) - stateIn(i,2)) > 1e-15 &
       .or. abs(stateOut(i,2) - stateIn(i,3)) > 1e-15 &
       .or. abs(stateOut(i,3) - stateIn(i,4)) > 1e-15 &
       .or. abs(stateOut(i,4) - stateIn(i,1)) > 1e-15 &
       .or. abs(stateOut(i,5) - stateIn(i,4)) > 1e-15 &
      ) then
      unitPass = .false.
      pass = .false.
    endif
  end do

  call ReportAndClear(unitPass, stateOut, derivOut, hnumpts,'nearest neighbour',&
                      testCount)








  call Interp1D_Linear (hNumPts, nPtsSrc, nPtsTgt, &
                                  hDim, hDim, &

                                  levelsIn, stateIn, derivIn, &

                                  interpIndex, levelsOut, stateOut, derivOut, &

                                  .false., .false., &
                                  lapseRateDown, lapseRateUp &
                                 )

  call TestReportClear(levelsSrc, stateIn, derivIn, levelsTgt, stateOut, &
                       derivOut, hnumpts, &
                       'linear', pass, testCount)








  call Interp1D_CubicLagrange (hNumPts, nPtsSrc, nPtsTgt, &
                                  hDim, hDim, &

                                  levelsIn, stateIn, derivIn, &

                                  interpIndex, levelsOut, stateOut, derivOut, &

                                  .false., .false., &
                                  lapseRateDown, lapseRateUp &
                                 )

  call TestReportClear(levelsSrc, stateIn, derivIn, levelsTgt, stateOut, &
                      derivOut, hnumpts, &
                      'cubic Lagrange', pass, testCount)








  call Extrap1D_LapseRate (hNumPts, nPtsSrc, nPtsTgt, &
                                  hDim, hDim, &

                                  levelsIn, stateIn, derivIn, &

                                  interpIndex, levelsOut, stateOut, derivOut, &

                                  .true., .true., &
                                  lapseRateDown, lapseRateUp &
                                 )

  ! Test this case individually
  unitPass = .true.
  if(     abs(stateOut(1,4) -  0.5411954)  > 1e-7 &
     .or. abs(stateOut(1,5) -  0.1058001)  > 1e-7 &
     .or. abs(stateOut(2,4) -  0.6885440)  > 1e-7 &
     .or. abs(stateOut(2,5) -(-0.2382999)) > 1e-7 &
    ) then
    unitPass = .false.
    pass = .false.
  endif

  call ReportAndClear(unitPass, stateOut, derivOut, hnumpts, &
                      'extrapolation by lapse rate', testCount)








  ! Test one example of the extended interface
! external FluxStub   ! at the top of this file -- any routine would do

  numExtArraysIn = 0
  numExtArraysOut = 0
  call Interp1D_CubicLagrange_X (hNumPts, nPtsSrc, nPtsTgt, &
                                  hDim, hDim, &

                                  levelsIn, stateIn, derivIn, &

                                  interpIndex, levelsOut, stateOut, derivOut, &

                                  .false., .false., &
                                  lapseRateDown, lapseRateUp, &

                                  FluxStub, &
                                  numExtArraysIn, numExtArraysOut, &
                                  ExtArraysIn, ExtArraysOut &
                                 )

  call TestReportClear(levelsSrc, stateIn, derivIn, levelsTgt, stateOut, &
                      derivOut, hnumpts, &
                      'cubic Lagrange_X', pass, testCount)







! TO TEST Extrap1D_Abort, UNCOMMENT THESE LINES:
! TO TEST Extrap1D_Abort, UNCOMMENT THESE LINES:
!!$  write(6,*)" "
!!$  write(6,*)"Testing Extrap1D_Abort.  With these data it should abort ..."
!!$  call Extrap1D_Abort (hNumPts, nPtsSrc, nPtsTgt, &
!!$                                  hDim, hDim, &
!!$
!!$                                  levelsIn, stateIn, derivIn, &
!!$
!!$                                  interpIndex, levelsOut, stateOut, derivOut, &
!!$
!!$                                  .true., .true., &
!!$                                  lapseRateDown, lapseRateUp &
!!$                                 )
!!$  write(6,*)"... Oops!:   Extrap1D_LapseRate did not abort"
!!$  write(6,*)" "
!!$  call exit(13)








  write(6,*)" "
  write(6,*)"Testing Extrap1D_Abort. With these data it should NOT abort ..."
  do i=1,hNumPts
    levelsOut(i,4) = 0.65
    levelsOut(i,5) = 2.9
  end do
  call Extrap1D_Abort (hNumPts, nPtsSrc, nPtsTgt, &
                                  hDim, hDim, &

                                  levelsIn, stateIn, derivIn, &

                                  interpIndex, levelsOut, stateOut, derivOut, &

                                  .true., .true., &
                                  lapseRateDown, lapseRateUp &
                                 )
  write(6,*)"... it didn't"









  ! Simulate temperature or humidity values at lowest source level above the
  ! surface (stateIn(i,2)) and at the surface (stateIn(i,1)), roughness length
  ! (z0), inverse of Monin-Obukhov length (ilmo) and height of the boundary layer
  ! (hBound).  Set the height (in m) of the lowest source level above the surface
  ! (levelsIn(i,2)) and of the surface itself (levelsIn(i,1)).

  hNumPts=6

  do i=1,hNumPts
     stateIn(i,1) = 0.
     stateIn(i,2) = 100.
     levelsIn(i,1) = 0.
     levelsIn(i,2) = 100.
     hBound(i)=300.
     z0(i)=1.0
     ilmo(i)=-0.01+i*0.01
  enddo

  ! Set the destination levels (levelsOut(i,v)).
  do v=1,nPtsTgt
    do i=1,hNumPts
      levelsOut(i,v)=5*(v-1.)
    enddo
  enddo

  ! Perform a quality control on hBound
  do i=1,hNumPts
     hBound(i)=max(hmin,hBound(i),(stateIn(i,2)+2*z0(i))*factn)
  enddo

  ! Calculate the vertical gradient of temperature or humidity at the lowest
  ! level above the surface along with the flux (ft).
  call FluxStub(f,stateIn(1,2),z0,ilmo,hBound,hNumPts) 
  do i=1,hNumPts
    ft(i)=(stateIn(i,2)-stateIn(i,1))/f(i)
  enddo

  ! Agglomerate the input / output arrays into the single input / output
  ! extension arrays
  numExtArraysIn = 4
  numExtArraysOut = 0
  ExtArraysIn(:,1) = z0
  ExtArraysIn(:,2) = ilmo
  ExtArraysIn(:,3) = hBound
  ExtArraysIn(:,4) = ft

  derivOut=0.0


  call Interp1D_FindPos(hNumPts, 2, nPtsTgt, &
                        hDim, hDim, &

                        levelsIn, interpIndex, levelsOut &
                       )

  call Extrap1D_Surface_X (hNumPts, 2, nPtsTgt, &
                                  hDim, hDim, &

                                  levelsIn, stateIn, derivIn, &

                                  interpIndex, &
                                  levelsOut, stateOut, derivOut, &

                                  .true., .true., &
                                  lapseRateDown, lapseRateUp, &

                                  FluxStub, &
                                  numExtArraysIn, numExtArraysOut, &
                                  ExtArraysIn, ExtArraysOut &
                                 )

  ! Test this case individually
  unitPass = .true.
  if(     abs(stateOut(1,4) - 60.0761986)  > 1e-7 &
     .or. abs(stateOut(1,5) - 65.9684296)  > 1e-7 &
     .or. abs(stateOut(2,2) - 23.9199448)  > 2e-6 &
     .or. abs(stateOut(2,3) - 33.9498825)  > 4e-6 &
    ) then
    unitPass = .false.
    pass = .false.
  endif

  call ReportAndClear(unitPass, stateOut, derivOut, hnumpts, &
                      'extrapolation near surface', testCount)









  ! Simulate wind components at lowest source level above the surface
  ! (stateIn(i,2)) and at the surface (stateIn(i,1)), latitude (in radians),
  ! roughness length (z0), inverse of Monin-Obukhov length (ilmo) and height of
  ! the boundary layer (hBound).  Set the height (in m) of the lowest source
  ! level above the surface (levelsIn(i,2)) and of the surface itself
  ! (levelsIn(i,1)).
  do i=1,hNumPts
     stateIn(i,2) = 10.
     v_windIn(i)  = 0.
     latitude(i)  = 1.0
     levelsIn(i,1) = 0.
     levelsIn(i,2) = 80.
     z0(i)=1.0
     ilmo(i)=-0.05+i*0.05
     hBound(i)=100./max(ilmo(i),1.e-9)
  enddo

  ! Set the destination levels (levelsOut(i,v)).
  do v=1,nPtsTgt
    do i=1,hNumPts
      levelsOut(i,v)=12.*v
    enddo
  enddo

  ! Calculate the vertical gradient of the wind components at the lowest
  ! level above the surface along with the stress (ft) and wind angle.
  call sltop_uvStub(ft,dudz,dvdz,tha,thmax,stateIn(:,2),v_windIn, &
                    levelsIn(:,2),z0,ilmo,hBound,latitude,hNumPts)

  ! Clear stateIn, just to be sure that the routines do NOT use it
  stateIn = 0.

  ! Agglomerate the input / output arrays into the single input / output
  ! extension arrays
  numExtArraysIn = 6
  numExtArraysOut = 2*nPtsTgt
  ExtArraysIn(:,1) = z0
  ExtArraysIn(:,2) = ilmo
  ExtArraysIn(:,3) = hBound
  ExtArraysIn(:,4) = ft
  ExtArraysIn(:,5) = tha
  ExtArraysIn(:,6) = latitude


  call Interp1D_FindPos(hNumPts, 2, nPtsTgt, &
                        hDim, hDim, &

                        levelsIn, interpIndex, levelsOut &
                       )

  call Extrap1D_SurfaceWind_X (hNumPts, 2, nPtsTgt, &
                                  hDim, hDim, &

                                  levelsIn, stateIn, derivIn, &

                                  interpIndex, &
                                  levelsOut, stateOut, derivOut, &

                                  .true., .true., &
                                  thmax, lapseRateUp, &

                                  slfun_uvStub, &
                                  numExtArraysIn, numExtArraysOut, &
                                  ExtArraysIn, ExtArraysOut &
                                 )

  ! Test this case individually
  unitPass = .true.
!
! N.B.:  Some of these values change by 2e-6, depending on the optimization
!
! ExtArraysOut contains the y-components of the wind
  if(     abs(stateOut    (1,4) - 8.8562183)  > 2e-6 &
     .or. abs(ExtArraysOut(1,4) - 0.0000000)  > 2e-6 &
     .or. abs(stateOut    (1,5) - 9.3546963)  > 2e-6 &
     .or. abs(ExtArraysOut(1,5) - 0.0000000)  > 2e-6 &
     .or. abs(stateOut    (2,2) - 5.4435577)  > 5e-6 &
     .or. abs(ExtArraysOut(2,2) - 0.1090328)  > 2e-6 &
     .or. abs(stateOut    (2,3) - 6.7105956)  > 5e-6 &
     .or. abs(ExtArraysOut(2,3) - 0.1056033)  > 2e-6 &
    ) then
    unitPass = .false.
    pass = .false.
  endif

  call ReportAndClear(unitPass, stateOut, ExtArraysOut, hnumpts, &
                      'wind extrapolation near surface', testCount)







  ! Report the overall test status
  if(pass) then
    write(6,*)' '
    write(6,*)'*  T E S T   P A S S E D  *'
  else
    write(6,*)'* * * * * * * * * * * * * *'
    write(6,*)'*                         *'
    write(6,*)'*  T E S T   F A I L E D  *'
    write(6,*)'*                         *'
    write(6,*)'* * * * * * * * * * * * * *'
  end if
  write(6,*)'   number of tests completed:  ', testCount

end program Test_Interp1D








!!s/r TestReportClear - test and report the results, and clear them for the next
!                       trial
subroutine TestReportClear(levelIn, stateIn, derivIn, level, state, deriv, &
                           hnumpts, title, passOverall, testCount)
!
!AUTHOR
!     J.W. Blezius JULY 2002
!
!REVISION
! v1_0    Blezius J.W.          - new code
!
!OBJECT
!        To avoid repeating code.
!
!ARGUMENTS
  implicit none
  real, dimension(nPtsSrc), intent(in) :: levelIn
  real, dimension(hDim, nPtsSrc), intent(in) :: stateIn, derivIn
  real, dimension(nPtsTgt), intent(in) :: level
  real, dimension(hDim, nPtsTgt), intent(inout) :: state, deriv
  integer, intent(in) :: hnumpts
  character(*), intent(in) :: title
  logical, intent(inout) :: passOverall
  integer, intent(inout) :: testCount
!
!NOTES
!  This routine assumes that the levelIn, stateIn, and derivIn arrays are in
!  ascending order
!
!!

  integer :: i, v
  logical :: pass

  pass = .true.

  !
  ! Check clamped results (non-derivative and derivative)
  !
  do v=1, nPtsTgt
    if( level(v) < levelIn(1) ) then
      ! the state should be clamped to the minimum input state
      do i=1,hNumPts
        if( abs(state(i,v) - stateIn(i,1)) > 1.e-15 ) pass = .false.
        if(       abs(deriv(i,v) - derivIn(i,1)) > 1.e-15 &
            .and. abs(deriv(i,v)) > 1.e-15 &
          ) pass = .false.
      end do

    else if( level(v) > levelIn(nPtsSrc) ) then
      ! the state should be clamped to the minimum input state
      do i=1,hNumPts
        if( abs(state(i,v) - stateIn(i,nPtsSrc)) > 1.e-15 ) pass = .false.
        if( deriv(i,v) /= derivIn(i,nPtsSrc) .and. abs(deriv(i,v)) > 1.e-15 ) &
                                                                   pass = .false.
      end do

  !
  ! Check all other non-derivative results
  !
    else ! verify for correct interpolation
      ! i=1 is sin
      if( abs( state(1,v) - sin(level(v)) ) > 0.06) pass = .false.

      ! i=2 is tan      d/dx tan = cos**(-2)
      if(abs( state(2,v) - tan(level(v)) ) > 0.5/cos(level(v))**2) pass = .false.
    end if
  end do

  !
  ! Check all other derivative results
  !
  if( abs(deriv(1,1)) > 1e-15 ) then
    ! the derivatives have been estimated:  verify them
    do v=1, nPtsTgt
      ! i=1 derivative is cos
      if( abs( deriv(1,v) - cos(level(v)) ) > 0.1) pass = .false.

      ! i=2 derivative is cos**(-2):  its derivative varies as tan
      if( abs( deriv(2,v) - 1./cos(level(v))**2 ) > 0.7*abs(tan(level(v)))) &
                                                                   pass = .false.
    end do
  end if !derivatives estimated

  passOverall = passOverall .and. pass

  call ReportAndClear(pass, state, deriv, hnumpts, title, testCount)
  
end subroutine TestReportClear







!!s/r Report - report the results, and clear them for the next trial
subroutine ReportAndClear(pass, state, deriv, hnumpts, title, testCount)
!
!AUTHOR
!     J.W. Blezius JULY 2002
!
!REVISION
! v1_0    Blezius J.W.          - new code
!
!OBJECT
!        To avoid repeating code.
!
!ARGUMENTS
  implicit none
  logical, intent(in) :: pass
  real, dimension(hDim, nPtsTgt), intent(inout) :: state, deriv
  integer, intent(in) :: hnumpts
  character(*), intent(in) :: title
  integer, intent(inout) :: testCount
!
!!
  integer :: i, v

  character(len=19), parameter :: spass  = '  PASS             ', &
                                  sfail  = '  **** F A I L ****'
  character(len=19) :: string

  testCount = testCount + 1
  if(pass) then
    string = spass
  else
    string = sfail
  endif


  write(6,*)
  write(6,*) title, string
  write(6,*)'stateOut='
  do i=1,hNumPts
    write(6,'((10f12.7))')(state(i,v), v=1,nPtsTgt)
  end do
  write(6,*)'derivOut='
  do i=1,hNumPts
    write(6,'((10f12.7))')(deriv(i,v), v=1,nPtsTgt)
  end do

  ! Clear the destination tables
  do v=1, nPtsTgt
    do i=1, hNumPts
      state(i,v)=0.
      deriv(i,v)=0.
    end do
  end do

end subroutine ReportAndClear







!!s/r FluxStub - Testing stub for the flux routine
subroutine FluxStub(f,zz,z0,ilmo,h,n)
!
!AUTHOR
!     J.W. Blezius OCT 2003
!
!REVISION
! v1_3    Blezius J.W.          - new code
!
!OBJECT
!        To fill the role, in part, of the physics library during testing.
!
!ARGUMENTS
  implicit none
      REAL F(N),ZZ(N),Z0(N),ILMO(N),H(N)
      INTEGER N
!
!NOTES
!        This routine was supplied by Yves Delage, who called it fh.
!
!!
      INTEGER J
      REAL BETA,CI,AS,ASX,FACTN,HH

      PARAMETER ( BETA = 1.0 , &
                  CI   = 40. , &
                  AS   = 12. , &
                  ASX  = 5.0 , &
                  FACTN= 1.2  )

      real LZZ0,Y,Y0,RAC3
      real a,b,c,d,psi,z
      real unsl,hi

!************************************************************************
!**  fonctions de couche de surface pour le cas stable                 **
!************************************************************************

      d  (unsl) = 4*AS*BETA*unsl
      c  (hi)   = d(unsl)*hi - hi**2
      b  (hi)   = d(unsl) - 2*hi
      a  (z,hi) = sqrt(1 + b(hi)*z - c(hi)*z**2)
      psi(z,hi) = 0.5 * (a(z,hi)-z*hi-log(1+b(hi)*z*0.5+a(z,hi))- &
                  b(hi)/(2*sqrt(c(hi)))*asin((b(hi)-2*c(hi)*z)/d(unsl)))

!   Limites de validite: unsl >= 0 (cas stable ou neutre)
!                        c > 0 (hi < d)
!                        z*hi < 1
!   Ces 2 dernieres conditions imposees a l'aide du facteur 'factn'
!
!   Reference :  Y. Delage, BLM, 82 (p23-48) (Eq.33-37)
!***********************************************************************
      RAC3=SQRT(3.)

      DO J=1,N
      LZZ0=LOG(ZZ(J)/Z0(J)+1)
      IF(ILMO(J).LE.0.) THEN
!---------------------------------------------------------------------
!                      UNSTABLE CASE
           Y=(1-BETA*CI*(ZZ(J)+Z0(J))*ILMO(J))**(1./3)
           Y0=(1-BETA*CI*Z0(J)*ILMO(J))**(1./3)
           F(J)=BETA*(LZZ0+1.5*ALOG((Y0**2+Y0+1)/(Y**2+Y+1))+RAC3* &
              ATAN(RAC3*2*(Y-Y0)/((2*Y0+1)*(2*Y+1)+3)))
      ELSE
!---------------------------------------------------------------------
!                        STABLE CASE
           unsl=ilmo(j)
        hi=1/MAX(H(J),factn/d(ILMO(J)))
           !write(6,*)'j=',j,'lzz0=',lzz0
           !write(6,*)'     min(psi...)=',min(psi(ZZ(J)+Z0(J),hi)-psi(Z0(J),hi),&
           !                ASX*BETA*ILMO(J)*ZZ(J))
           f(j)=BETA*(LZZ0+min(psi(ZZ(J)+Z0(J),hi)-psi(Z0(J),hi), &
                           ASX*BETA*ILMO(J)*ZZ(J)))
      ENDIF
!---------------------------------------------------------------------
      END DO

      return
      end subroutine FluxStub








!!s/r sltop_uvStub - Testing stub for the sltop_uv routine
subroutine sltop_uvStub(nss,dudz,dvdz,angtop,angmax, &
                        utop,vtop,ztop,z0,ilmo,h,lat,n)
!
!AUTHOR
!     J.W. Blezius OCT 2003
!
!REVISION
! v1_3    Blezius J.W.          - new code
!
!OBJECT
!        To fill the role, in part, of the physics library during testing.
!
!ARGUMENTS
  implicit none
      REAL NSS(N),DUDZ(N),DVDZ(N),ANGTOP(N), &
           UTOP(N),VTOP(N),ZTOP(N),Z0(N),ILMO(N),H(N),LAT(N)
      INTEGER N
!
!          - Output -
! NSS     normalised surface stress
! DUDZ    slope of the U component of wind profile at top of SL
! DVDZ    slope of the V component of wind profile at top of SL
! ANGTOP  wind direction at top of SL
! ANGMAX  maximum wind direction change between surface and H
!
!          - Input -
! UTOP    U component of wind at the top of surface layer
! VTOP    V component of wind at the top of surface layer
! ZTOP    height of the top of the surface layer
! Z0      roughness length for wind
! ILMO    inverse of MONIN-OBUKHOV lenth
! H       height of boundary layer (for stable case only)
! LAT     latitude in radians
! N       number of horizontal points to process
!
!
!NOTES
!        This routine was supplied by Yves Delage, who called it sltop_uv.
!
!!
      external slfun_uvStub
      INTEGER J
      
      REAL      FM(N)
      REAL AS,ASX,CI,BS,BETA,FACTN,HMIN,ANGMAX,SPEED

! Initilisation of constants in the common SURFCON

      AS    = 12.
      ASX   = 13.
      CI    = 40.
      BS    = 1.0
      BETA  = 1.0
      FACTN = 1.2
      HMIN  = 40.
      ANGMAX= 0.85

      call slfun_uvStub(fm,ztop,z0,ilmo,h,n)
      DO J=1,N
        speed=sqrt(utop(j)**2+vtop(j)**2)
        angtop(j)=atan2(vtop(j),sign(abs(utop(j))+1.e-05,utop(j)))
        nss(j)=speed/fm(j)
        dudz(j)=nss(j)*phim(ztop(j),ilmo(j),h(j))*cos(angtop(j)) &
                              /(ztop(j)+z0(j)) &
                + speed*sin(angtop(j))*angmax*sin(lat(j))/h(j)
        dvdz(j)=nss(j)*phim(ztop(j),ilmo(j),h(j))*sin(angtop(j)) &
                              /(ztop(j)+z0(j)) &
                - speed*cos(angtop(j))*angmax*sin(lat(j))/h(j)
      END DO

      return
      CONTAINS
!  Derivatives of the stability functions

        REAL FUNCTION PHIM(Z,ILMO,H)

        REAL Z,ILMO,H,HH

        HH=MAX(1-Z/H,FACTN-1.0)
        IF(ILMO.GT.0.) THEN
           PHIM=MIN(1.+ASX*BETA*Z*ILMO,0.5*(HH+SQRT(HH**2+ &
                4.*AS*BETA*Z*ILMO*HH)))
        ELSE
           PHIM=(1.-CI*BETA*Z*ILMO)**(-.1666666)
        END IF
        RETURN
        END FUNCTION



        REAL FUNCTION PHIH(Z,ILMO,H)

        REAL Z,ILMO,H,HH

        HH=MAX(1-Z/H,FACTN-1.0)
        IF(ILMO.GT.0.) THEN
           PHIH=BETA*MIN(1.+ASX*BETA*Z*ILMO,0.5*(HH+SQRT(HH**2+ &
                4.*AS*BETA*Z*ILMO*HH)))
        ELSE
           PHIH=(1.-CI*BETA*Z*ILMO)**(-.333333333)
        END IF
        RETURN
        END FUNCTION
      end subroutine sltop_uvStub








!!s/r slfun_uvStub - Testing stub for the slfun_uvStub routine
subroutine slfun_uvStub(fm,z,z0,ilmo,h,n)
!
!AUTHOR
!     J.W. Blezius OCT 2003
!
!REVISION
! v1_3    Blezius J.W.          - new code
!
!OBJECT
!        To fill the role, in part, of the physics library during testing.
!
!ARGUMENTS
  implicit none
      REAL FM(N),Z(N),Z0(N),ILMO(N),H(N)
      INTEGER N
!
!          - Output -
! FM      normalised wind speed at desired output height
!
!          - Input -
! Z       height of desired output
! Z0      roughness length for wind
! ILMO    inverse of MONIN-OBUKHOV lenth
! H       height of boundary layer (for stable case only)
! N       number of horizontal points to process
!
!NOTES
!        This routine was supplied by Yves Delage, who called it slfun_uv.
!
!!
      INTEGER J
      REAL AS,ASX,CI,BS,BETA,FACTN,HMIN,ANGMAX

      REAL RAC3,X,X0,Y,Y0,Z0T(1),HI,LZZ0T(1)
      
      REAL      LZZ0(N)

! Initilisation of constants in the common SURFCON

      AS    = 12.
      ASX   = 13.
      CI    = 40.
      BS    = 1.0
      BETA  = 1.0
      FACTN = 1.2
      HMIN  = 40.
      ANGMAX= 0.85

      RAC3=SQRT(3.)
      DO J=1,N
      LZZ0(J)=LOG(Z(J)/Z0(J)+1)
      IF(ILMO(J).LE.0.) THEN
!---------------------------------------------------------------------
!                      UNSTABLE CASE
           fm(j)= fmi(z(j)+z0(j),j)
      ELSE
!---------------------------------------------------------------------
!                        STABLE CASE
        hi=1/MAX(H(J),hmin,factn/(4*AS*BETA*ilmo(j)), &
             (Z(J)+10*Z0(J))*factn)
        fm(j)=LZZ0(J)+min(psi(Z(J)+Z0(J),j)-psi(Z0(J),j), &
                           ASX*BETA*ILMO(J)*Z(J))
      ENDIF
!---------------------------------------------------------------------
      END DO

      return

      CONTAINS
!   Internal function FMI
!   Stability function for momentum in the unstable regime (ilmo<0)
!   Reference: Delage Y. and Girard C. BLM 58 (19-31) Eq. 19

      REAL FUNCTION FMI(Z,I)

      REAL Z
      integer i

      X=(1-CI*Z*BETA*ILMO(I))**(0.1666666)
      X0=(1-CI*Z0(I)*BETA*ILMO(I))**(0.1666666)
      FMI=LZZ0(I)+LOG((X0+1)**2*SQRT(X0**2-X0+1)*(X0**2+X0+1)**1.5 &
                     /((X+1)**2*SQRT(X**2-X+1)*(X**2+X+1)**1.5)) &
                    +RAC3*ATAN(RAC3*((X**2-1)*X0-(X0**2-1)*X)/ &
                    ((X0**2-1)*(X**2-1)+3*X*X0))

      RETURN
      END FUNCTION FMI

!   Internal function FHI
!   Stability function for heat and moisture in the unstable regime (ilmo<0)
!   Reference: Delage Y. and Girard C. BLM 58 (19-31) Eq. 17

      REAL FUNCTION FHI(Z,I)

      REAL Z
      integer i

      Y=(1-CI*Z*BETA*ILMO(I))**(0.33333333)
      Y0=(1-CI*Z0T(I)*BETA*ILMO(I))**(0.3333333)
      FHI=BETA*(LZZ0T(I)+1.5*LOG((Y0**2+Y0+1)/(Y**2+Y+1))+RAC3* &
              ATAN(RAC3*2*(Y-Y0)/((2*Y0+1)*(2*Y+1)+3)))

      RETURN
      END FUNCTION FHI

!   Internal function psi
!   Stability function for momentum in the stable regime (unsl>0)
!   Reference :  Y. Delage, BLM, 82 (p23-48) (Eqs.33-37)


      REAL FUNCTION PSI(Z,I)

      REAL Z,a,b,c,d
      integer i

      d = 4*AS*BETA*ilmo(i)
      c = d*hi - hi**2
      b = d - 2*hi
      a = sqrt(1 + b*z - c*z**2)
      psi = 0.5 * (a-z*hi-log(1+b*z*0.5+a)- &
                  b/(2*sqrt(c))*asin((b-2*c*z)/d))

      RETURN
      END FUNCTION PSI
      end subroutine slfun_uvStub
