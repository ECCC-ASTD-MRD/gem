 
      integer,parameter :: PHYINREAD_MAX = 512
      integer(INT64) :: phyinread_jdateo
      integer :: phyinread_n, phyinread_dt, phyinread_step
      integer :: phyinread_list_nk(PHYINREAD_MAX)
      character(len=32) :: phyinread_list_S(PHYINREAD_MAX)

      common /phyinread8/ phyinread_jdateo
      common /phyinread/  phyinread_n, phyinread_dt, phyinread_step
      common /phyinreada/ phyinread_list_nk
      common /phyinreads/ phyinread_list_S
