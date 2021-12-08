
module phybus
   public
   save

#define PHYPTR
#include "phyvar.hf"

   integer :: sigw=0  !#set to sigm/t in sigmalev

   real, pointer, contiguous :: entbus(:,:) => null()
   real, pointer, contiguous :: dynbus(:,:) => null()
   real, pointer, contiguous :: perbus(:,:) => null()
   real, pointer, contiguous :: volbus(:,:) => null()

end module phybus

