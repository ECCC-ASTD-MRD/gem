#ifdef AIX_POWERPC7

#define VPOW_WRAPPERS

      subroutine vpow1n_wrap(r,x,y,n)
!!!#include <arch_specific.hf>
      real(REAL64) :: x,y(*),r(*)
      real(REAL64) :: w(n)
      w(:)=x
      Call vpow(r,w,y,n)
      end

      subroutine vpown1_wrap(r,x,y,n)
!!!#include <arch_specific.hf>
      real(REAL64) :: x(*),y,r(*)
      real(REAL64) :: w(n)
      w(:)=y
      Call vpow(r,x,w,n)
      end

      subroutine vpownn_wrap(r,x,y,n)
!!!#include <arch_specific.hf>
      real(REAL64) :: x(*),y(*),r(*)
      Call vpow(r,x,y,n)
      end

      subroutine vspow1n_wrap(r,x,y,n)
!!!#include <arch_specific.hf>
      real x,y(*),r(*)
      real w(n)
      w(:)=x
      Call vspow(r,w,y,n)
      end

      subroutine vspown1_wrap(r,x,y,n)
!!!#include <arch_specific.hf>
      real x(*),y,r(*)
      real w(n)
      w(:)=y
      Call vspow(r,x,w,n)
      end

      subroutine vspownn_wrap(r,x,y,n)
!!!#include <arch_specific.hf>
      real x(*),y(*),r(*)
      Call vspow(r,x,y,n)
      return
      end

#else

      subroutine vpow_dummy()
      return
      end

#endif
