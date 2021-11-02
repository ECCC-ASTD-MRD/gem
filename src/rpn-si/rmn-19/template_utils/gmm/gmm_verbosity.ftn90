#ifdef INTERFACEONLY
      interface
      integer function gmm_verbosity(verbose_level)
      integer, intent(in) :: verbose_level
      end function gmm_verbosity
      end interface
#endif

#ifndef INTERFACEONLY
  integer function gmm_verbosity(verbose_level)
  use gmm_internals
  implicit none
  integer, intent(in) :: verbose_level

   select case (verbose_level)
   case (GMM_MSG_DEBUG)
      gmm_verbose_level = GMM_MSG_DEBUG
      print *, '(GMM_VERBOSITY) Setting GMM Verbose level to DEBUG'
   case (GMM_MSG_INFO)
      gmm_verbose_level = GMM_MSG_INFO
      print *, '(GMM_VERBOSITY) Setting GMM Verbose level to INFO'
   case (GMM_MSG_WARN)
      gmm_verbose_level = GMM_MSG_WARN
      print *, '(GMM_VERBOSITY) Setting GMM Verbose level to WARN'
   case (GMM_MSG_ERROR)
      gmm_verbose_level = GMM_MSG_ERROR
      print *, '(GMM_VERBOSITY) Setting GMM Verbose level to ERROR'
   case (GMM_MSG_SEVERE)
      gmm_verbose_level = GMM_MSG_SEVERE
      print *, '(GMM_VERBOSITY) Setting GMM Verbose level to SEVERE'
   case (GMM_MSG_FATAL)
      gmm_verbose_level = GMM_MSG_FATAL
      print *, '(GMM_VERBOSITY) Setting GMM Verbose level to FATAL'
   case default
      print *, '(GMM_VERBOSITY) Unknown GMM Verbose level'
      print *, '(GMM_VERBOSITY) Please check GMM Documentation'
   end select
  gmm_verbosity = 0
!
  return

  end function gmm_verbosity
#endif
