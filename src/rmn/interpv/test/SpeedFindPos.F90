program SpeedFindPos
  implicit none 

  integer, parameter :: NNI=500         ! scales the arrays involved
  real,    dimension(NNI,50) :: vLevelSource, &
                                vLevelDestn
  integer, dimension(NNI,50) :: posnDestInSrc

  integer :: h                          ! horizontal loop index
  integer :: v                          ! vertical loop index


  ! Fill in some source and destination levels
  do v=1,50
    do h=1,NNI
      vLevelSource(h,v)=v
      vLevelDestn (h,v)=v + 0.2
    end do ! h
  end do ! v

  ! Process 200,000 horizontal points
  do h=1,200000/NNI
    call Interp1D_FindPos(NNI, 50, 50, &
                          NNI, NNI, &

                          vLevelSource, &
                          posnDestInSrc, &
                          vLevelDestn &
                         )
  end do

  write(6,*) 'Done'

end program SpeedFindPos
