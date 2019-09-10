
module phybus
   public
   save

#define PHYPTR
#include "phyvar.hf"

   integer :: sigw=0  !#set to sigm/t in sigmalev

   real, pointer :: entbus(:,:) => null()
   real, pointer :: dynbus(:,:) => null()
   real, pointer :: perbus(:,:) => null()
   real, pointer :: volbus(:,:) => null()

end module phybus

