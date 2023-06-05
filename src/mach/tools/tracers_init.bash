#!/bin/bash
################################### LICENCE BEGIN ###############################
# GEM-MACH - Atmospheric chemistry library for the GEM numerical atmospheric model
# Copyright (C) 2007-2013 - Air Quality Research Division &
#                           National Prediction Operations division
#                           Environnement Canada
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
################################### LICENCE END #################################

# les analyses sont ici: /data/aqli05/afsudev/gem-mach/input/trunk/100x100/analyses
# ./tracers_init.bash /cnfs/dev/aq/aq01/afsuhul/2011052000_000 . 2010011100 2010011112 12 8

. r.ssmuse.dot profile
. s.ssmuse.dot SCM/1.0.6

set -x
last_run_file=$1  # full path and file name
inrep=$2          # where the met files are located
start_date=$3     # yyyymmddhh (assuming that delta >= 1h)
end_date=$4       # yyyymmddhh
pilot_delta=$5    # in hour
max_cpu=${6:-8}

tracers_list="TA2 TA3 TALD TARO TBZO TC38 TCO TCR1 TCR2 TCRE TDIA TETH TH22 THCH THN3 THN4 THO2 THON TISO TMC3 TMEK TMGL TN25 TNH3 TNO2 TNO3 TNO TO3 TO TOH TOSD TPAN TR22 TR2N TR2R TRN3 TRO2 TROO TSO2 TSO4 TTOL TSU1 TSU2 TSS1 TSS2 TOC1 TOC2 TNI1 TNI2 TAM1 TAM2 TCM1 TCM2 TEC1 TEC2 TPC1 TPC2"

EDITFST=editfst_6.10
fichier_directives="directives"
traceurs_de_00=les_traceurs_00.fst

date_fst=$(r.date ${start_date})

nb_heures=0
nb_heures=$(r.date -nV ${end_date} ${start_date})
(( modulo = nb_heures % pilot_delta ))
if [[ ${modulo} -ne 0 ]]; then
   echo "ERREUR: L'increment d'heure doit etre un facteur de la duree de l'integration" >&2
   echo "nb_heures = ${nb_heures}, pilot_delta = ${pilot_delta}, modulo = ${modulo}" >&2
   exit 1
fi

echo "Nombre d'heure a traiter: ${nb_heures}"
if [[ $nb_heures -lt 1 ]]; then
   echo "ERREUR: le nombre d'heure a integrer doit etre plus grand que 0.  nb_heures = ${nb_heures}" >&2
   exit 1
fi

# creation du fichier qui contiendra les traceurs a copier dans tous les fichiers de pilotages
# Premierement, voyon si les traceurs sont deja dans le fichier de pilotage 000 (ce sera le cas suivant gmpilot ou son equivalent dans les series)
traceurs_present=$(r.fstliste -izfst ${inrep}/${start_date}_000|cut -b 1-4|sort -u)
va_voir_ailleur=0
le_fichier_source=${inrep}/${start_date}_000
for traceur in ${tracers_list}
do
   est_il_la=$(echo ${traceurs_present}|grep ${traceur})
   if [[ $? -ne 0 ]]; then
      echo "Le champ ${traceur} est introuvable dans le fichier d'analyse meteo ${inrep}/${start_date}_000"
      echo "On va verifier si les traceurs sont dans le fichier de sorties (argument 1: ${last_run_file})"
      va_voir_ailleur=1
      break
   fi
done
# S'il n'est pas present dans le fichier 000, on regarde dans la sortie de la derniere run
if [[ ${va_voir_ailleur} -eq 1 ]]; then
   traceurs_present=$(r.fstliste -izfst ${last_run_file}|cut -b 1-4)
   le_fichier_source=${last_run_file}
   for traceur in ${tracers_list}
   do
      est_il_la=$(echo ${traceurs_present}|grep ${traceur})
      if [[ $? -ne 0 ]]; then
         echo "ERREUR: le champ ${traceur} est introuvable dans le fichier de sortie ${last_run_file}" >&2
         exit 1
      fi
   done
fi

# On genere maintenant les fichiers de directives
# Comme il y a un maximum de 20 directives par fichier, on doit en creer plusieurs
max_directive=20
compteur_directive=1
index_fichier=0
echo "" > ${fichier_directives}_${index_fichier}
for traceur in ${tracers_list}
do   
   if [[ ${compteur_directive} -gt 20 ]]; then
      compteur_directive=1
      (( index_fichier = index_fichier + 1 ))
      echo "" > ${fichier_directives}_${index_fichier}
   fi
   echo "desire (-1,['${traceur}'])" >> ${fichier_directives}_${index_fichier}
   (( compteur_directive = compteur_directive + 1 ))
done

while [[ ${index_fichier} -ge 0 ]];
do
   $EDITFST -s ${le_fichier_source} -d ${traceurs_de_00} -i ${fichier_directives}_${index_fichier}
   (( index_fichier = index_fichier - 1 ))
done

mv ${traceurs_de_00}.fd ${traceurs_de_00}

#boucle sur toutes les dates ou le pilotage est requis
if [[ ${va_voir_ailleur} -eq 1 ]]; then
#    le fichier 000 ne contient pas la chimie, on commence la
   heure_courante=0
else
#    le fichier 000 contenait la chimie, on commence a l'increment suivant
   heure_courante=${pilot_delta}
fi
cpu_count=0
while [[ ${heure_courante} -le ${nb_heures} ]]
do
(
   # On change la date de validite des champs meteo avec forcedate
   (( npas_courant = 8 * heure_courante ))
   for traceur in ${tracers_list}
   do
   # forcedate INPUT OUTPUT FIELD_NAME DATE DEET NPAS
      forcedate ${traceurs_de_00} ${traceurs_de_00}.${heure_courante} ${traceur} ${date_fst} 450 ${npas_courant}
   done
   suffixe_fichier=$(printf %03i ${heure_courante})
   echo "Ajout de ${traceurs_de_00}.${heure_courante} dans ${start_date}_${suffixe_fichier}"
   $EDITFST -s ${traceurs_de_00}.${heure_courante} -d ${start_date}_${suffixe_fichier} -i 0
# << EOZ
#    desire(-1)
#    zap (-1, -1, -1, ${date_fst}, -1, ${heure_courante}, -1)   
# EOZ
# Et maintenant on copie les tic-tac du fichier original
   $EDITFST -s ${le_fichier_source} -d ${start_date}_${suffixe_fichier} -i << EOD
   desire(-1, ['>>','^^'])
EOD
) &
   (( heure_courante = heure_courante + pilot_delta ))
   (( cpu_count = cpu_count + 1 ))
   if [ "$cpu_count" = "${max_cpu}" ]; then
      echo "On attends que les ${cpu_count} threads aient termine" 
      wait
      echo "Threads termine, on reprend"
      cpu_count=0
   fi
done



