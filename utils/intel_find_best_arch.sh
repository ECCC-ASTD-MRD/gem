#!/bin/sh

# Since the Intel compilers do not support "-march=native", we try to find the most
# performant architecture

for arch in avx512 avx2 avx; do
    grep flags /proc/cpuinfo | head -n 1 | grep -q $arch && break
done

if [[ "$arch" = "avx512" ]]; then
    echo -n "skylake-avx512"
elif [[ "$arch" = "avx2" ]]; then
    echo -n "core-avx2"
elif [[ "$arch" = "avx" ]]; then
    echo -n "core-avx-i"
else
    # Bail out at core2
   echo -n "core2"
fi
