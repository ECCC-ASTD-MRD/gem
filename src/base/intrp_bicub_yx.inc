   interface                                                                                      !InTf!
    subroutine set_intrp_bicub_off_yx(qx,qy) bind(C,name='set_intrp_bicub_off_yx')                !InTf!
      integer, intent(IN), value :: qx,qy                                                         !InTf!
    end subroutine set_intrp_bicub_off_yx                                                         !InTf!
    subroutine set_intrp_bicub_off_quv(qx,qy,ux,uy,vx,vy) bind(C,name='set_intrp_bicub_off_quv')  !InTf!
      integer, intent(IN), value :: qx,qy,ux,uy,vx,vy                                             !InTf!
    end subroutine set_intrp_bicub_off_quv                                                        !InTf!
    subroutine set_intrp_bicub_yx(x1,x2,y1,y2) bind(C,name='set_intrp_bicub_yx')                  !InTf!
      integer, intent(IN), value :: x1,x2,y1,y2                                                   !InTf!
    end subroutine set_intrp_bicub_yx                                                             !InTf!
    subroutine set_intrp_bicub_quv(x1,x2,y1,y2,u1,u2,u3,u4,v1,v2,v3,v4) bind(C,name='set_intrp_bicub_quv')   !InTf!
      integer, intent(IN), value :: x1,x2,y1,y2,u1,u2,u3,u4,v1,v2,v3,v4                           !InTf!
    end subroutine set_intrp_bicub_quv                                                            !InTf!
    subroutine intrp_bicub_yx_s(f, r, ni, ninj, nk, x, y) bind(C,name='intrp_bicub_yx_s')         !InTf!
      import :: C_INT, C_FLOAT, C_DOUBLE                                                          !InTf!
      real(C_FLOAT), dimension(*), intent(IN) :: f                                                !InTf!
      real(C_FLOAT), dimension(*), intent(OUT) :: r                                               !InTf!
      real(C_DOUBLE), intent(IN), value :: x, y                                                   !InTf!
      integer(C_INT), intent(IN), value :: ni, ninj, nk                                           !InTf!
    end subroutine intrp_bicub_yx_s                                                               !InTf!
    subroutine intrp_bicub_yx_s_mono(f, r, ni, ninj, nk, x, y) bind(C,name='intrp_bicub_yx_s_mono') !InTf!
      import :: C_INT, C_FLOAT, C_DOUBLE                                                          !InTf!
      real(C_FLOAT), dimension(*), intent(IN) :: f                                                !InTf!
      real(C_FLOAT), dimension(*), intent(OUT) :: r                                               !InTf!
      real(C_DOUBLE), intent(IN), value :: x, y                                                   !InTf!
      integer(C_INT), intent(IN), value :: ni, ninj, nk                                           !InTf!
    end subroutine intrp_bicub_yx_s_mono                                                          !InTf!
    subroutine intrp_bicub_yx_vec(u,v,ru,rv,ni,ninj,nk,xu,yu,xv,yv,s) bind(C,name='intrp_bicub_yx_vec')  !InTf!
      import :: C_INT, C_FLOAT, C_DOUBLE                                                          !InTf!
      real(C_FLOAT), dimension(*), intent(IN) :: u, v                                             !InTf!
      real(C_FLOAT), dimension(*), intent(OUT) :: ru, rv                                          !InTf!
      real(C_DOUBLE), intent(IN), value :: xu, yu, xv, yv                                         !InTf!
      real(C_DOUBLE), intent(IN), dimension(*) :: s                                               !InTf!
      integer(C_INT), intent(IN), value :: ni, ninj, nk                                           !InTf!
    end subroutine intrp_bicub_yx_vec                                                             !InTf!
    subroutine intrp_bicub_yx_uv(u,v,r1,r2,ni,ninj,nk,x,y,s) bind(C,name='intrp_bicub_yx_uv')     !InTf!
      import :: C_INT, C_FLOAT, C_DOUBLE                                                          !InTf!
      real(C_FLOAT), dimension(*), intent(IN) :: u, v                                             !InTf!
      real(C_FLOAT), dimension(*), intent(OUT) :: r1, r2                                          !InTf!
      real(C_DOUBLE), intent(IN), value :: x, y                                                   !InTf!
      real(C_DOUBLE), intent(IN), dimension(*) :: s                                               !InTf!
      integer(C_INT), intent(IN), value :: ni, ninj, nk                                           !InTf!
    end subroutine intrp_bicub_yx_uv                                                              !InTf!
   end interface                                                                                  !InTf!
