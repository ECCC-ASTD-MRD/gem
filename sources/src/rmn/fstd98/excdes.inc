        integer, parameter :: EXCDES_DESIRE = 1
        integer, parameter :: EXCDES_EXCLURE = -1
        integer, parameter :: EXCDES_RANGE = -2
        integer, parameter :: EXCDES_DELTA = -3
        integer, parameter :: MAX_REQUETES = 20

      interface

!       int fst_deactivate_filters()
        function fst_deactivate_filters() result (old_status)
          import :: C_INT
          integer(C_INT) :: old_status
        end function fst_deactivate_filters

!       int fst_reactivate_filters()
        function fst_reactivate_filters() result (old_status)
          import :: C_INT
          integer(C_INT) :: old_status
        end function fst_reactivate_filters

!       int Xc_Select_ip1/2/3(int set_nb, int des_exc, void *iplist, int nelm)
        function Select_ip1(set_nb, des_exc, iplist, nelm)                     &
     &           result (status) BIND(C,name='Xc_Select_ip1')
          import :: C_INT, C_CHAR
          implicit none
          integer(C_INT) :: status
          integer(C_INT), intent(IN), value :: set_nb, des_exc, nelm
          integer(C_INT), dimension(nelm), intent(IN) :: iplist
        end function Select_ip1

        function Select_ip2(set_nb, des_exc, iplist, nelm)                     &
     &           result (status) BIND(C,name='Xc_Select_ip2')
          import :: C_INT, C_CHAR
          implicit none
          integer(C_INT) :: status
          integer(C_INT), intent(IN), value :: set_nb, des_exc, nelm
          integer(C_INT), dimension(nelm), intent(IN) :: iplist
        end function Select_ip2

        function Select_ip3(set_nb, des_exc, iplist, nelm)                     &
     &           result (status) BIND(C,name='Xc_Select_ip3')
          import :: C_INT, C_CHAR
          implicit none
          integer(C_INT) :: status
          integer(C_INT), intent(IN), value :: set_nb, des_exc, nelm
          integer(C_INT), dimension(nelm), intent(IN) :: iplist
        end function Select_ip3

!       int Xc_Select_date(int set_nb, int des_exc, int *date_list, int nelm)
        function Select_date(set_nb, des_exc, date_list, nelm)                 &
     &           result (status) BIND(C,name='Xc_Select_date')
          import :: C_INT, C_CHAR
          implicit none
          integer(C_INT) :: status
          integer(C_INT), intent(IN), value :: set_nb, des_exc, nelm
          integer(C_INT), dimension(nelm), intent(IN) :: date_list
        end function Select_date

!       int Xf_select_etiquette(int set_nb, int des_exc, char *etiq_list, int nelm, int flng)
        function Select_etiquette(set_nb, des_exc, etiq_list,                  &
     &                            nelm, flng)                                  &
     &           result (status) BIND(C,name='Xf_Select_etiquette')
          import :: C_INT, C_CHAR
          implicit none
          integer(C_INT) :: status
          integer(C_INT), intent(IN), value :: set_nb,des_exc,nelm,flng
          character(C_CHAR), dimension(*), intent(IN) :: etiq_list
        end function Select_etiquette

!       int Xc_Select_nomvar(int set_nb, int des_exc, char *nomv_list[], int nelm)
        function Select_nomvar(set_nb, des_exc, nomv_list, nelm, flng)         &
     &           result (status) BIND(C,name='Xf_Select_nomvar')
          import :: C_INT, C_CHAR
          implicit none
          integer(C_INT) :: status
          integer(C_INT), intent(IN), value :: set_nb,des_exc,nelm,flng
          character(C_CHAR), dimension(*), intent(IN) :: nomv_list
        end function Select_nomvar

!       int Xc_Select_typvar(int set_nb, int des_exc, char *typv_list[], int nelm)
        function Select_typvar(set_nb, des_exc, typv_list, nelm, flng)         &
     &           result (status) BIND(C,name='Xf_Select_typvar')
          import :: C_INT, C_CHAR
          implicit none
          integer(C_INT) :: status
          integer(C_INT), intent(IN), value :: set_nb,des_exc,nelm,flng
          character(C_CHAR), dimension(*), intent(IN) :: typv_list
        end function Select_typvar

!        int Xc_Select_suppl(int set_nb, int des_exc, int ni, int nj, int nk, int ig1, int ig2, int ig3, int ig4, char gtyp)
        function Select_suppl(set_nb, des_exc, ni, nj, nk,                     &
     &                        ig1, ig2, ig3, ig4, gtyp)                        &
     &           result (status) BIND(C,name='Xc_Select_suppl')
          import :: C_INT, C_CHAR
          implicit none
          integer(C_INT) :: status
          integer(C_INT), intent(IN), value :: set_nb,des_exc,ni,nj,nk
          integer(C_INT), intent(IN), value :: ig1, ig2, ig3, ig4
          character(C_CHAR), intent(IN), value :: gtyp
        end function Select_suppl

        subroutine Dump_Request_Table() BIND(C,name='DumpRequestTable')
        end subroutine Dump_Request_Table
        
!       int c_fst_match_req(int handle)
        function fst_match_req(handle) result (status)                         &
     &         BIND(C,name='C_fst_match_req')
          import :: C_INT, C_CHAR
          implicit none
          integer(C_INT) :: status
          integer(C_INT), intent(IN), value :: handle
        end function fst_match_req

!       int C_requetes_reset(int set_nb, int nomvars, int typvars, int etikets, int dates, int ip1s, int ip2s, int ip3s)
        function f_requetes_reset(set_nb,nomvars,typvars,etikets,              &
     &         dates,ip1s,ip2s,ip3s) result(status)                            &
     &         BIND(C,name='C_requetes_reset')
          import :: C_INT, C_CHAR
          implicit none
          integer(C_INT) :: status
          integer(C_INT), intent(IN), value :: set_nb,nomvars,typvars,         &
     &                  etikets,dates,ip1s,ip2s,ip3s
        end function f_requetes_reset

      end interface

