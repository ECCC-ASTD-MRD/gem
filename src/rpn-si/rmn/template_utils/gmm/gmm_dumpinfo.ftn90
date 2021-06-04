#ifdef INTERFACEONLY
       interface
         subroutine gmm_dumpinfo(fldstat)
         logical, intent(in), optional :: fldstat
         end subroutine gmm_dumpinfo
       end interface
#endif

#ifndef INTERFACEONLY
    subroutine gmm_dumpinfo(fldstat)
    use gmm_internals
    implicit none
    logical, intent(in), optional :: fldstat
    integer :: i, l_page, l_entry, nelm, crc, f_calc_crc
    real :: xx
    pointer(px, xx(*))
    type(gmm_layout), dimension(4) :: dims
    external f_calc_crc

    l_page = 1
    l_entry = 1
    print *,'GMM dumpinfo, number of variables in use is',used
    do i=1,used

      dims = directory(l_page)%entry(l_entry)%l
      nelm = ( (dims(1)%high - dims(1)%low +1) * &
              &(dims(2)%high - dims(2)%low +1) * &
              &(dims(3)%high - dims(3)%low +1) * &
              &(dims(4)%high - dims(4)%low +1) )
      if (present(fldstat)) then
        print *,'Appel a statfld a ecrire, fldstat=',fldstat
        write (6,77) 'Name=',directory(l_page)%entry(l_entry)%name,' addr=',&
                      &directory(l_page)%entry(l_entry)%array_addr
      else
        call make_cray_pointer(px,directory(l_page)%entry(l_entry)%array_addr)
        crc = f_calc_crc(xx,nelm,0,1)
        write (6,88) 'Name=',directory(l_page)%entry(l_entry)%name,' addr=',&
                      &directory(l_page)%entry(l_entry)%array_addr,' checksum=',crc 
      endif
      l_entry = l_entry + 1
      if (l_entry .gt. PAGE_SIZE) then
        l_page = l_page + 1
        l_entry = 1
      endif

    enddo
77  format(a,a,a,i10)
88  format(a,a,a,i10,a,i10)
    return
    end subroutine gmm_dumpinfo
#endif