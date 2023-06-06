#! /bin/bash
## 
## typically the surface and major point emissions are grouped within their respective 
## subdirectories, and there are (12 x 7) files in each of these subdirectories, one for
## each month/day combination
##
## each file carries a name of the form ${prefix}_${month}_${day}${suffix}
##
## the emission task of the GM maestro suite expects the following structure
##
## area_directory/month/area_shortmonth_day.fst 
## and
## major_pt_directory/month/major_shortmonth_day.fst 
##
## where month takes the values {january february march april may june july august september october november december}
##       shortmonth, the values {jan feb mar apr may jun jul aug sep oct nov dec}
## and   day, the values {sun mon tue wed thu fri sat}
##
## this script can help you generated the required structure starting from a set of existing emissions files
## that have name of the form ${prefix}_${month}_${day}${suffix}. Soft links are used.
##
##

#  set -x
## $dir_src is the directory where both the area and major point subdirectories are located
## $dir_target is the directory where the new structure will be located
dir_src=/fs/cmo/data/emissions_can2006_us2012_10km
dir_target=/fs/dev/mrb/armn/armnsyg/emissions/emissions_can2006_us2012_10km
month=(january february march april may june july august september october november december)
short_m=(jan feb mar apr may jun jul aug sep oct nov dec)
day=(sun mon tue wed thu fri sat)

## in $dir_src, the area emissions files are in the subdirectory ${area_src}
## within that directory the files carry names = ${area_prefix}_${area_m}_${area_d}${area_suffix}
area_src=area
area_prefix=area
area_suffix=.fst_int10km
area_m=(jan fev mar avr mai jun jul aou sep oct nov dec)
area_d=(dim lun mar mer jeu ven sam)
## in $dir_src, the major point emissions files are in the subdirectory ${mjr_src}
## within that directory the files carry names = ${mjr_prefix}_${major_m}_${major_d}${mjr_suffix}
mjr_src=major
mjr_prefix=2011
mjr_suffix=_major_GEMMACH.fst
major_m=(Jan Feb March April May June July Aug Sept Oct Nov Dec)
major_d=(Sunday Monday Tuesday Wednesday Thursday Friday Saturday Sunday)

n=0
   while (test "${n}" -le 11); do
    the_month=${month[n]}
    abv_month=${short_m[n]}
    mkdir -p ${dir_target}/area/${the_month}
    mkdir -p ${dir_target}/major/${the_month}

nday=0
    while (test "${nday}" -le 6); do
    AREA=${dir_src}/${area_src}/${area_prefix}_${area_m[n]}_${area_d[nday]}${area_suffix}
    my_area=area_${abv_month}_${day[nday]}.fst
    ln -s -f ${AREA} ${dir_target}/area/${the_month}/${my_area}

    MAJOR=${dir_src}/${mjr_src}/${mjr_prefix}_${major_m[n]}_${major_d[nday]}${mjr_suffix}
    my_major=major_${abv_month}_${day[nday]}.fst
    ln -s -f ${MAJOR} ${dir_target}/major/${the_month}/${my_major}

    let nday=nday+1
    done 

   let n=n+1
   done

exit

