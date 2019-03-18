!/@*
!@objective Provide basic set of constants and interfaces to librmn fn
!*@/

!- TYPEDEF
integer, parameter :: RDOUBLE = 8
integer, parameter :: IDOUBLE = 8

!- STATUS (error is <0; use min() to collect, use if <RMN_OK to check)
integer, parameter :: RMN_ERR = -1
integer, parameter :: RMN_OK  = 0

!- I/O
integer, parameter :: RMN_STDOUT = 6
integer, parameter :: RMN_STDERR = 0
integer, parameter :: RMN_PATH_LEN = 512 !adjusted to WB_MAX_STRLEN, should move back to 1024 for longer path
integer, parameter :: RMN_ETK_LEN = 12
integer, parameter :: RMN_VARNAME_LEN = 4
integer, parameter :: RMN_VARTYPE_LEN = 2

integer, parameter :: RMN_DTYPE_UINT = 2 !-unsigned int
integer, parameter :: RMN_DTYPE_SINT = 4 !-signed int
integer, parameter :: RMN_DTYPE_IEEE = 5

logical, parameter :: RMN_REWRITE = .true.
logical, parameter :: RMN_APPEND  = .false.

integer, parameter :: RMN_CONV_31BITS = 0
integer, parameter :: RMN_CONV_IP2P = -1
integer, parameter :: RMN_CONV_P2IP = 1
integer, parameter :: RMN_CONV_P2IPNEW = 2
integer, parameter :: RMN_CONV_P2IPOLD = 3

integer, parameter :: RMN_CONV_MSL = 0
integer, parameter :: RMN_CONV_SIGMA = 1
integer, parameter :: RMN_CONV_PRESS = 2
integer, parameter :: RMN_CONV_ARBITRARY = 3
integer, parameter :: RMN_CONV_MGL = 4
integer, parameter :: RMN_CONV_HY = 5
integer, parameter :: RMN_CONV_TH = 6


logical, parameter :: RMN_CONV_USEFORMAT_L = .true.

logical, parameter :: RMN_OPT_GET = .true.
logical, parameter :: RMN_OPT_SET = .false.

!- newdate mode
integer, parameter :: RMN_DATE_STAMP2TRUE = 1
integer, parameter :: RMN_DATE_TRUE2STAMP = -1
integer, parameter :: RMN_DATE_PRINT2TRUE = 2
integer, parameter :: RMN_DATE_TRUE2PRINT = -2
integer, parameter :: RMN_DATE_PRINT2STAMP = 3
integer, parameter :: RMN_DATE_STAMP2PRINT = -3
integer, parameter :: RMN_DATE_OLD2STAMP = 4
integer, parameter :: RMN_DATE_STAMP2OLD = -4

character(len=RMN_ETK_LEN), parameter :: RMN_ANY_ETK = '            '
character(len=RMN_VARTYPE_LEN), parameter :: RMN_ANY_TYP = '  '
character(len=RMN_VARNAME_LEN), parameter :: RMN_ANY_NAME = '    '
integer, parameter :: RMN_ANY_I  = -1
integer, parameter :: RMN_ANY_DATE  = -1


!-
interface

   integer function fstopc(F_option_S,F_value_S,F_get_L)
      character(len=*),intent(in) :: F_option_S
      character(len=*) :: F_value_S !I/O
      logical,intent(in) :: F_get_L
   end function fstopc

   integer function exdb(F_title,F_version,F_ipcode)
      character(len=*),intent(in) :: F_title,F_version,F_ipcode
   end function exdb

   integer function exfin(F_title,F_version,F_ipcode)
      character(len=*),intent(in) :: F_title,F_version,F_ipcode
   end function exfin

   integer function wkoffit(F_path)
      character(len=*),intent(in) :: F_path
   end function wkoffit

   integer function fnom(F_unit,F_filename_S,F_type_S,F_recl)
      integer,intent(inout) :: F_unit
      character(len=*),intent(in) :: F_filename_S,F_type_S
      integer,intent(in) :: F_recl
   end function fnom

   integer function fstouv(F_unit,F_type_S)
      integer,intent(in) :: F_unit
      character(len=*),intent(in) :: F_type_S
   end function fstouv

   integer function fstlnk(F_unitlist,F_nunits)
      integer,intent(in) :: F_nunits      
      integer,intent(in) :: F_unitlist(F_nunits)
   end function fstlnk

   integer function fstfrm(F_unit)
      integer,intent(in) :: F_unit
   end function fstfrm

   integer function fclos(F_unit)
      integer,intent(in) :: F_unit
   end function fclos

   integer function fstprm(F_key, F_dateo,F_deet,F_npas, F_ni,F_nj,F_nk, &
                   F_nbits, F_datyp, F_ip1, F_ip2, F_ip3, &
                   F_typvar, F_nomvar, F_etiket, &
                   F_grtyp, F_ig1, F_ig2, F_ig3, F_ig4, F_swa, F_lng, F_dltf, &
                   F_ubc, F_extra1, F_extra2, F_extra3)
      integer,intent(in) :: F_key
      integer,intent(out) :: F_dateo,F_deet,F_npas,F_ni,F_nj,F_nk, F_nbits, &
           F_datyp, F_ip1, F_ip2, F_ip3, F_ig1, F_ig2, F_ig3, F_ig4, &
           F_swa, F_lng, F_dltf,F_ubc, F_extra1, F_extra2, F_extra3
      character(len=*),intent(out) :: F_typvar,F_nomvar, F_etiket, F_grtyp
   end function fstprm

   integer function fstinf(F_unit,F_ni,F_nj,F_nk,F_datev,F_etiket,F_ip1,F_ip2,F_ip3,F_typvar,F_nomvar)
      integer,intent(out) :: F_ni,F_nj,F_nk
      integer,intent(in) :: F_unit,F_datev,F_ip1,F_ip2,F_ip3
      character(len=*),intent(in) :: F_etiket,F_typvar,F_nomvar
   end function fstinf

   integer function fstinl(F_unit,F_ni,F_nj,F_nk,F_datev,F_etiket,F_ip1,F_ip2,F_ip3,F_typvar,F_nomvar,F_liste,F_nbfound,F_nmax)
      character(len=*),intent(in) :: F_etiket,F_typvar,F_nomvar
      integer,intent(in) :: F_unit,F_datev,F_ip1,F_ip2,F_ip3,F_nmax
      integer,intent(out) :: F_ni,F_nj,F_nk
      integer,intent(out) :: F_nbfound,F_liste(F_nmax)
   end function fstinl

   integer function ezqkdef(ni, nj, grtyp, ig1, ig2, ig3, ig4, iunit)
      integer,intent(in) :: ni, nj, ig1, ig2, ig3, ig4, iunit
      character(len=*),intent(in) :: grtyp
   end function ezqkdef

!!$   integer function ezgdef_fmem(ni,nj, grtyp, grref, ig1,ig2,ig3,ig4, ax, ay)
!!$      integer,intent(in) :: ni, nj, ig1, ig2, ig3, ig4
!!$      character(len=*),intent(in) :: grtyp,grref
!!$      real,intent(in) :: ax(:), ay(:)
!!$   end function ezgdef_fmem

   integer function ezdefset(gdout, gdin)
      integer,intent(in) :: gdout, gdin
   end function ezdefset

   integer function ezsetopt(option, value)
      character(len=*),intent(in) :: option, value
   end function ezsetopt

!!$   integer function ezsint(zout, zin)
!!$      real,intent(out) :: zout(*)
!!$      real,intent(in) :: zin(*)
!!$   end function ezsint
!!$
!!$   integer function ezuvint(uuout, vvout, uuin, vvin)
!!$      real,intent(out) :: uuout(*), vvout(*)
!!$      real,intent(in) :: uuin(*), vvin(*)
!!$   end function ezuvint
   
   integer function gdrls(gridid)
      integer,intent(in) :: gridid
   end function gdrls

   integer function ezgprm(gridid, grtyp_S, ni, nj, ig1, ig2, ig3, ig4)
      integer,intent(in) :: gridid
      character(len=*),intent(out) :: grtyp_S
      integer,intent(out) :: ni, nj, ig1, ig2, ig3, ig4
   end function ezgprm

   integer function ezgxprm(gridid, ni, nj, grtyp_S, ig1, ig2, ig3, ig4, grref_S, ig1ref, ig2ref, ig3ref, ig4ref)
      integer,intent(in) :: gridid
      character(len=*),intent(out) :: grtyp_S,grref_S
      integer,intent(out) :: ni, nj, ig1, ig2, ig3, ig4, ig1ref, ig2ref, ig3ref, ig4ref
   end function ezgxprm

!!$   integer function gdgaxes(gridid, ax, ay)
!!$      integer,intent(in) :: gridid
!!$      real,intent(out) :: ax(:), ay(:)
!!$   end function gdgaxes

   integer function ezgfstp(gridid, nomvarx_S, typvarx_S, etikx_S, nomvary_S, typvary_S, etiky_S, ip1, ip2, ip3, dateo, deet, npas, nbits)
      integer,intent(in) :: gridid
      character(len=*),intent(out) :: nomvarx_S, typvarx_S, etikx_S, nomvary_S, typvary_S, etiky_S
      integer,intent(out) :: ip1, ip2, ip3, dateo, deet, npas, nbits
   end function ezgfstp

   subroutine convip(ip,p,kind,mode,format_S,useformat_L)
      integer,intent(in) :: mode
      integer,intent(inout) :: ip, kind
      real ,intent(inout) ::  p
      character(len=*),intent(in) :: format_S
      logical,intent(in) :: useformat_L
   end subroutine convip

   integer function ip1_all(level,kind)
      real,intent(in) :: level
      integer,intent(in) :: kind
   end function ip1_all

   subroutine difdatr(idate1, idate2, nhours)
      real(8),intent(out) :: nhours  !date difference expressed in hours
      integer,intent(in)  :: idate1  !CMC date-time stamp integer
      integer,intent(in)  :: idate2  !CMC date-time stamp integer
   end subroutine difdatr

   subroutine incdatr(idate_out, idate_in, nhours)
      real(8),intent(in)  :: nhours
      integer,intent(in)  :: idate_in  !CMC date-time stamp integer
      integer,intent(out) :: idate_out !CMC date-time stamp integer
   end subroutine incdatr

   !TODO: newdate old style is integer,dimension(14)
!!$   integer function newdate(dat1, dat2, dat3, mode)
!!$      integer :: dat1, dat2, dat3
!!$      integer,intent(in) :: mode
!!$   end function newdate

   subroutine cxgaig(grdtyp_S,ig1,ig2,ig3,ig4,xg1,xg2,xg3,xg4)
      character(len=*),intent(in) :: grdtyp_S
      integer,intent(out) :: ig1,ig2,ig3,ig4
      real,intent(in) :: xg1,xg2,xg3,xg4
   end subroutine cxgaig

   subroutine cigaxg(grdtyp_S,xg1,xg2,xg3,xg4,ig1,ig2,ig3,ig4)
      character(len=*),intent(in) :: grdtyp_S
      real,intent(out) :: xg1,xg2,xg3,xg4
      integer,intent(in) :: ig1,ig2,ig3,ig4
   end subroutine cigaxg

end interface

integer,external :: fstecr, fstlir, fstluk
integer,external :: ezsint, ezuvint, gdgaxes,gdxyfll,gdllfxy,ezgdef_fmem
integer,external :: newdate
integer,external :: ezsetval,ezsetival
