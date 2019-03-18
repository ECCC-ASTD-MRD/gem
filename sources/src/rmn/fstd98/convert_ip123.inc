
      type, BIND(C)  :: FLOAT_IP
      real(C_FLOAT)  :: lo         ! low bound
      real(C_FLOAT)  :: hi         ! high bound
      integer(C_INT) :: kind       ! type code (see table below)
      end type

      integer, parameter :: TO_IP=1
      integer, parameter :: TO_RP=-1
      integer, parameter :: CONVERT_OK=0
      integer, parameter :: CONVERT_GUESS=14
      integer, parameter :: CONVERT_GOOD_GUESS=2
      integer, parameter :: CONVERT_BAD_GUESS=4
      integer, parameter :: CONVERT_TERRIBLE_GUESS=8
      integer, parameter :: CONVERT_WARNING=32
      integer, parameter :: CONVERT_ERROR=64

      integer, parameter :: KIND_ABOVE_SEA=0
      integer, parameter :: KIND_SIGMA=1
      integer, parameter :: KIND_PRESSURE=2
      integer, parameter :: KIND_ARBITRARY=3
      integer, parameter :: KIND_ABOVE_GND=4
      integer, parameter :: KIND_HYBRID=5
      integer, parameter :: KIND_THETA=6
      integer, parameter :: KIND_HOURS=10
      integer, parameter :: KIND_SAMPLES=15
      integer, parameter :: KIND_MTX_IND=17
      integer, parameter :: KIND_M_PRES=21

      interface encode_ip

      function encode_ip_0(IP1,IP2,IP3,RP1,RP2,RP3) RESULT(status)      &
     &                    BIND (C,name='EncodeIp')
        use ISO_C_BINDING
        type, BIND(C)  :: FLOAT_IP
        real(C_FLOAT)  :: lo         ! low bound
        real(C_FLOAT)  :: hi         ! high bound
        integer(C_INT) :: kind       ! type code (see table below)
        end type
        integer(C_INT) :: status
        integer(C_INT), intent(OUT) :: IP1,IP2,IP3
        type(FLOAT_IP), intent(IN)  :: RP1,RP2,RP3
      end FUNCTION encode_ip_0

      function encode_ip_1(IP,RP) result(status)                        &
     &                    BIND (C,name='EncodeIp_v')
        use ISO_C_BINDING
        type, BIND(C)  :: FLOAT_IP
        real(C_FLOAT)  :: lo         ! low bound
        real(C_FLOAT)  :: hi         ! high bound
        integer(C_INT) :: kind       ! type code (see table below)
        end type
        integer(C_INT) :: status
        integer(C_INT), dimension(3), intent(OUT) :: IP
        type(FLOAT_IP), dimension(3), intent(IN)  :: RP
      end function encode_ip_1

      function encode_ip_2(IP1,IP2,IP3,P1,kkind1,P2,kkind2,P3,kkind3)   &
     &                    result(status) BIND(C,name='ConvertPKtoIP')
        use ISO_C_BINDING
        integer(C_INT) :: status
        integer(C_INT),        intent(OUT) :: IP1,IP2,IP3
        real(C_FLOAT), value, intent(IN)   :: P1,P2,P3
        integer(C_INT), value, intent(IN)  :: kkind1,kkind2,kkind3
      end function encode_ip_2

      function encode_ip_3(IP,RP,kind) result(status)                   &
     &                    BIND(C,name='ConvertPKtoIP_v')
        use ISO_C_BINDING
        integer(C_INT) :: status
        integer(C_INT), dimension(3), intent(OUT) :: IP
        real(C_FLOAT), dimension(3), intent(IN)   :: RP
        integer(C_INT), dimension(3), intent(IN)  :: kind
      end function encode_ip_3

      end interface

      interface decode_ip

      function decode_ip_0(RP1,RP2,RP3,IP1V,IP2V,IP3V) result(status)   &
     &                    BIND (C,name='DecodeIp')
        use ISO_C_BINDING
        type, BIND(C)  :: FLOAT_IP
        real(C_FLOAT)  :: lo         ! low bound
        real(C_FLOAT)  :: hi         ! high bound
        integer(C_INT) :: kind       ! type code (see table below)
        end type
        integer(C_INT) :: status
        integer(C_INT), value, intent(IN)  :: IP1V,IP2V,IP3V
        type(FLOAT_IP), intent(OUT) :: RP1,RP2,RP3
      end function decode_ip_0

      function decode_ip_1(RP,IP) result(status)                        &
     &                    BIND (C,name='DecodeIp_v')
        use ISO_C_BINDING
        type, BIND(C)  :: FLOAT_IP
        real(C_FLOAT)  :: lo         ! low bound
        real(C_FLOAT)  :: hi         ! high bound
        integer(C_INT) :: kind       ! type code (see table below)
        end type
        integer(C_INT) :: status
        integer(C_INT), dimension(3), intent(IN)  :: IP
        type(FLOAT_IP), dimension(3), intent(OUT) :: RP
      end function decode_ip_1

      function decode_ip_2(RP1,kind1,RP2,kind2,RP3,kind3,IP1V,IP2V,IP3V)&
     &                    result(status) BIND(C,name='ConvertIPtoPK')
        use ISO_C_BINDING
        integer(C_INT) :: status
        real(C_FLOAT),        intent(OUT)  :: RP1,RP2,RP3
        integer(C_INT),        intent(OUT) :: kind1,kind2,kind3
        integer(C_INT), value, intent(IN)  :: IP1V,IP2V,IP3V
      end function decode_ip_2

      function decode_ip_3(RP,kind,IP) result(status)                   &
     &                    BIND(C,name='ConvertIPtoPK_v')
        use ISO_C_BINDING
        integer(C_INT) :: status
        real(C_FLOAT),  dimension(3), intent(OUT) :: RP
        integer(C_INT), dimension(3), intent(OUT) :: kind
        integer(C_INT), dimension(3), intent(IN)  :: IP
      end function decode_ip_3

      end interface

      include 'convip_plus.inc'

      interface

      subroutine convip_plus(ip, p, kind, mode, string, flagv)
        integer, intent(INOUT) :: ip, kind
        integer, intent(IN) :: mode
        real, intent(INOUT) :: p
        character *(*), intent(OUT) :: string 
        logical, intent(IN) :: flagv
      end subroutine convip_plus

      subroutine c_conv_ip(ip, p, kind, mode) BIND(C,name='ConvertIp') ! reduced inteerface with no string option
        use ISO_C_BINDING
        implicit none
          integer(C_INT), intent(INOUT) :: ip, kind
          integer(C_INT), intent(IN), value :: mode
          real(C_FLOAT), intent(INOUT) :: p
      end subroutine c_conv_ip

      end interface
