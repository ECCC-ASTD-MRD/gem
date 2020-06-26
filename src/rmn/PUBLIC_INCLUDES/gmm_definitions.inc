
      type gmm_layout                              ! dimensioning elements
         SEQUENCE
         integer :: low,high,halo,halomax,n      ! dimension is low:high, useful contents 1..n
      end type

      type gmm_attributes
        SEQUENCE
        integer*8 :: key          ! unique key, used mainly for vmmget interface
        integer*8 :: uuid1, uuid2 ! extra keys used for lookup in external tables or additional flags
        integer   :: initmode                   ! how the field must be initialized (default=-1, no init)
        integer   :: flags                      ! flags (valid, restart, .... )
      end type

      type gmm_metadata
        SEQUENCE
        type(gmm_layout), dimension(4) :: l
        type(gmm_attributes) :: a
      end type
