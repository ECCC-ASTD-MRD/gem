      integer, parameter                            :: GMM_NULL_FLAGS=0
      type(gmm_layout), parameter                   :: GMM_NULL_LAYOUT=gmm_layout(0,0,0,0,0)
      type(gmm_layout), parameter, dimension(4)     :: GMM_NULL_LAYOUTS = (/GMM_NULL_LAYOUT, GMM_NULL_LAYOUT, GMM_NULL_LAYOUT, GMM_NULL_LAYOUT/)
      type(gmm_attributes), parameter               :: GMM_NULL_ATTRIB=gmm_attributes(GMM_KEY_NOT_FOUND,0,0,0,0)
      type(gmm_metadata), parameter                 :: GMM_NULL_METADATA=gmm_metadata(GMM_NULL_LAYOUTS , GMM_NULL_ATTRIB)
