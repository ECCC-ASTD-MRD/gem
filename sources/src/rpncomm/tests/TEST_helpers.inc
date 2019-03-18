        subroutine TestUserInit(NX,NY) ! try to get NX,NY from file TEST.cfg if it exists
        external :: get_a_free_unit
        integer :: get_a_free_unit
        integer :: iun,ier,i
        iun=get_a_free_unit()
        if(iun<0)return
        print *,'attempting to read TEST.cfg'
        print *,'nx , ny =',nx,ny
        open(UNIT=iun,FILE='TEST.cfg',STATUS='OLD',ACTION='READ',IOSTAT=ier)
        print *,'open iostat=',ier
        if(ier .ne. 0) then
          print *,'attempting auto distribution, nx , ny =',nx,ny
          do i=5,2,-1
            print *,'nx i mod(nx,i) =',nx, i, mod(nx,i)
            if(mod(NX,i) == 0 .and. nx .ne. i)then
              NY = i
              NX = NX / i
              return
            endif
          enddo
          return
        endif
        read(UNIT=iun,IOSTAT=ier,FMT=*)NX,NY
        print *,'read iostat=',ier
        close(UNIT=iun)
        return
        end
        integer function get_a_free_unit()
        implicit none
        integer :: i
        character (len=16) :: access_mode
          get_a_free_unit=-1
          do i = 1 , 99  ! find an available unit number
            inquire(UNIT=i,ACCESS=access_mode)
            if(trim(access_mode) == 'UNDEFINED')then ! found
              get_a_free_unit = i
              exit
            endif
          enddo
        return
        end
