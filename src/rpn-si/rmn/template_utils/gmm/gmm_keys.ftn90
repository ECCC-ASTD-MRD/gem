#ifdef INTERFACEONLY
      interface
      integer function gmm_nkeys()
      end function gmm_nkeys
      end interface 
      interface
      integer function gmm_keys(taglist,pattern)
      character(len=*), intent(out) :: taglist(:)
      character(len=*), intent(in), optional :: pattern
      end function gmm_keys
      end interface
#endif

#ifndef INTERFACEONLY
  integer function gmm_nkeys()
  use gmm_internals
  implicit none
  gmm_nkeys = used
  return
  end function gmm_nkeys

  integer function gmm_keys(taglist,pattern)
  use gmm_internals
  implicit none
  character(len=*), intent(out) :: taglist(:)
  character(len=*), intent(in), optional :: pattern
  integer :: i,strlen_pattern,nkeys,maxkeys

  nkeys = 0
  maxkeys = nkeys
  gmm_keys = -1
!
  maxkeys = size(taglist)
  if (used > maxkeys) then
    return
  endif
  cur_page = 1
  cur_entry = 1
  if (present(pattern)) then
    strlen_pattern = len_trim(pattern)
    gmm_keys = 0
    do i=1,used
      if (directory(cur_page)%entry(cur_entry)%name(1:strlen_pattern) == pattern(1:strlen_pattern)) then
        taglist(nkeys+1) = directory(cur_page)%entry(cur_entry)%name
        nkeys = nkeys + 1
      endif
      cur_entry = cur_entry + 1
      if (cur_entry .gt. PAGE_SIZE) then
        cur_page = cur_page + 1
        cur_entry = 1
      endif
    enddo
    gmm_keys = nkeys
  else
    do i = 1,used
      taglist(i) = directory(cur_page)%entry(cur_entry)%name
      cur_entry = cur_entry + 1
      if (cur_entry .gt. PAGE_SIZE) then
        cur_page = cur_page + 1
        cur_entry = 1
      endif
    enddo
    gmm_keys = used
  endif
  return
#ifdef DEBUG_MODE
  print *,'(GMMKEYS) looking for pattern=',pattern,'='
#endif
  end  function gmm_keys
#endif

!