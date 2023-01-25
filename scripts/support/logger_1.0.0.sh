#===============================================================================
# Environnement Canada - Service meteorologique du Canada
# Centre meteorologique canadien
# 2121 Route Trans-canadienne
# Dorval, Quebec
# H9P 1J3
#
# Projet     : Logger to generalize the logging mechanism among the jobs.
# Nom        : <Logger.sh>
# Creation   : Juin 2009 - J.P. Gauhier - CMC/CMOE
#
# Description: Functions to generalize the loggin mechanism amongs various jobs:
#
# Global Settings
#   LOG_LEVEL      Logging level (MUST,ERROR,WARNING,INFO,DEBUG,EXTRA)
#   LOG_MAILTO     EMail address for mail reports
#   LOG_MAILTITLE  Title to include in mail reports
#   LOG_FILE       File to use for log (stdout if undefined)
#   LOG_TIME       Include dates in log
#   LOG_JOBID      Job Identificator
#   LOG_JOBCLASS   Job class (SCRIPT,DAEMON,ORJI,HCRON,INTERACTIVE,REPORT)
#   LOG_JOBREPORT  Job class (ALL,WARNING,ERROR)
#
# Functions:
#   log_start     $Job $Version [$Paramfile]
#   log_mail      $Subject $File
#   log_end       $ExitStatus
#   log_print     $Level $Message $Seconds
#   log_time      $Seconds
#
# Retour     :
#
# Remarques  :
#   Aucune.
#===============================================================================

#----- Check for environment settings
LOG_LEVEL=${LOG_LEVEL:=INFO}
LOG_COLOR=${LOG_COLOR:=0}
LOG_TIME=${LOG_TIME:=0}
LOG_PROC=${LOG_PROC:=0}
LOG_FILE=${LOG_FILE:=""}

LOG_MAILTO=${LOG_MAILTO:=""}
LOG_MAILTITLE=${LOG_MAILTITLE:="Job Info"}

LOG_JOB=${LOG_JOB:="Unknown"}
LOG_JOBVERSION=${LOG_JOBVERSION:="Unknown"}
LOG_JOBID=${LOG_JOBID:=""}
LOG_JOBDATE=$(date +%Y%m%d_%H%M%S)
LOG_JOBPATH=${LOG_JOBPATH:=""}
LOG_JOBCLASS=${LOG_JOBCLASS:=SCRIPT}
LOG_JOBREPORT=${LOG_JOBREPORT:=ALL}
LOG_JOBCOMMAND=$0
LOG_JOBARGS=$*

#----- Logger internal variables
LOG_SECTIME=$(date +%s)
LOG_SECSTART=$(date +%s)
LOG_SECEND=$(date +%s)
LOG_WARNINGS=0
LOG_ERRORS=0
LOG_RETRY=0

typeset -A LOG_LEVELS
LOG_LEVELS=( [MUST]=-1 [ERROR]=0 [WARNING]=1 [INFO]=2 [DEBUG]=3 [EXTRA]=4 )

typeset -A APP_COLORS
APP_COLORS=( [BLINK]="\x1b[5m" [BLACK]="\x1b[30m" [RED]="\x1b[31m" [GREEN]="\x1b[32m" [ORANGE]="\x1b[33m" [YELLOW]="\x1b[1m\x1b[33m" [BLUE]="\x1b[34m" [MAGENTA]="\x1b[35m" [CYAN]="\x1b[36m" [LIGHTCYAN]="\x1b[1m\x1b[36m" [GRAY]="\x1b[37m" [RESET]="\x1b[0m" )

typeset -A LOG_COLORS
LOG_COLORS=( [MUST]="" [ERROR]="${APP_COLORS[RED]}" [WARNING]="${APP_COLORS[YELLOW]}" [INFO]="" [DEBUG]="${APP_COLORS[LIGHTCYAN]}" [EXTRA]="${APP_COLORS[CYAN]}" [PROGRESS]="${APP_COLORS[MAGENTA]}" )

#----------------------------------------------------------------------------
# Nom      : <log_start>
# Creation : Octobre 2009 - J.P. Gauthier - CMC/CMOE
#
# But      : Afficher une message de demarrage standard.
#
# Parametres  :
#    <Job>    : Nom de la job
#    <Version>: Version de la job
#    <Input>  : fichier d'entre (afin de recupere le temps d'attente en queue)
#
# Retour:
#
# Remarques :
#----------------------------------------------------------------------------
log_start() {
   LOG_SECSTART=$(date +%s)
   LOG_JOB=${1}
   LOG_JOBVERSION=${2}
   in=${3}

   if [[ ${LOG_JOBID} = "" ]] ; then
      LOG_JOBID=$(echo ${LOG_JOB} | tr '[a-z]' '[A-Z]')
   fi
   
   #----- Simulation run time ID.
   LOG_JOBID="${LOG_JOBID}-${LOG_JOBDATE}"

   if [[ -n "$LOG_FILE" ]] ; then
      rm -f ${LOG_FILE}
   fi

   log_print MUST "-------------------------------------------------------------------------------"
   log_print MUST "Script              : ${LOG_JOB}"
   log_print MUST "Version             : ${LOG_JOBVERSION}"
   log_print MUST "Hostname            : $(hostname)"
   log_print MUST "Architecture        : ${ORDENV_PLAT}"
   log_print MUST "Run ID              : ${LOG_JOBID}"

   if [[ ${LOG_MAILTO} != "" ]] ; then
      log_print MUST "E-mail Address      : ${LOG_MAILTO}"
   fi

   #----- Queue stuff
   if [[ -n "$LOADL_STEP_ID" ]]; then
      log_print MUST "Queue Method        : llv"
      log_print MUST "   Queue            : $LOADL_STEP_CLASS"
      log_print MUST "   Job ID           : $LOADL_STEP_ID"

      if [[ -r ${in} ]]; then
         secs0=$(date +%s)
         eval secs=\`perl -e \'\$mtime=\(stat\(\"${in}\"\)\)\[9\]\; print ${secs0}-\$mtime\'\`
         log_print MUST "   Waiting time     : $(log_time ${secs})"
      fi
    elif [[ -n "$SGE_CELL" ]]; then
      log_print MUST "Queue Method        : sge"
      log_print MUST "   Queue               : $QUEUE"
      log_print MUST "   Job ID              : $JOB_ID"
      if [[ -r ${in} ]]; then
         secs0=$(date +%s)
         eval secs=\`perl -e \'\$mtime=\(stat\(\"${in}\"\)\)\[9\]\; print ${secs0}-\$mtime\'\`
         log_print MUST "   Waiting time     : $(log_time ${secs})"
      fi
    elif [[ -n "$PBS_JOBID" ]]; then
      log_print MUST "Queue Method        : pbs"
      log_print MUST "   Queue               : $PBS_O_QUEUE"
      log_print MUST "   Job ID              : $PBS_JOBID"
      if [[ -r ${in} ]]; then
         secs0=$(date +%s)
         eval secs=\`perl -e \'\$mtime=\(stat\(\"${in}\"\)\)\[9\]\; print ${secs0}-\$mtime\'\`
         log_print MUST "   Waiting time     : $(log_time ${secs})"
      fi
   fi

   log_print MUST "Start time          : $(date +%c)"
   log_print MUST "-------------------------------------------------------------------------------\n"

   if [[ ${LOG_JOBCLASS} = "INTERACTIVE" ]]; then
      log_mail "Job started" ${LOG_FILE}
   fi

   #----- Starting hook
   log_hookstart
}

#----------------------------------------------------------------------------
# Nom      : <log_end>
# Creation : Octobre 2009 - J.P. Gauthier - CMC/CMOE
#
# But      : Afficher une message de fin standard.
#
# Parametres  :
#    <Status> : Code de retour de la job (0=ok, sinon erreur)
#    <Exit>   : Sortie du prorgamme (Default=True)
#
# Retour:
#
# Remarques :
#----------------------------------------------------------------------------
log_end() {
   status=${1}
   exitnow=${2}

   if [[ ${status} -eq -1 ]]; then
      if [[ ${LOG_ERRORS} -gt 0 ]]; then
         status=1
      else
         status=0
      fi
   fi

   LOG_SECEND=$(date +%s)

   log_print MUST "\n-------------------------------------------------------------------------------"
   if [[ ${status} -eq 0 ]]; then
      log_print MUST "Status              : Job has terminated successfully (${LOG_WARNINGS} Warning(s))."
   else
      log_print MUST "Status              : Job has encountered some errors (${LOG_ERRORS} Error(s))."
   fi
   log_print MUST "End time            : $(date +%c)"
   log_print MUST "Total running time  : $(log_time $((LOG_SECEND-LOG_SECSTART)))"
   log_print MUST "-------------------------------------------------------------------------------"

   if [[ ${LOG_JOBCLASS} = "REPORT" ]]; then
      if [[ ${LOG_ERRORS} -gt 0 ]]; then
         log_mail "Job finished (ERROR (${LOG_ERRORS}))" ${LOG_FILE}
      elif [[ ${LOG_WARNINGS} -gt 0 ]]; then
         log_mail "Job finished (WARNING (${LOG_WARNINGS}))" ${LOG_FILE}
      elif [[ ${LOG_JOBREPORT} == "ALL" ]]; then
          log_mail "Job finished (NORMAL)" ${LOG_FILE}
      fi
   elif [[ ${status} -eq 0 ]]; then
      if [[ ${LOG_JOBCLASS} = "INTERACTIVE" ]]; then
         log_mail "Job finished (NORMAL)" ${LOG_FILE}
      fi
   else
      log_mail "Job finished (ERROR (${LOG_ERRORS}))" ${LOG_FILE}
   fi
      
   #----- Ending hook
   log_hookend $status

   if [[ ${exitnow} = "" ]]; then
      exit $status
   fi
}

#----------------------------------------------------------------------------
# Nom      : <log_print>
# Creation : Octobre 2009 - J.P. Gauthier - CMC/CMOE
#
# But      : Afficher une message standard.
#
# Parametres  :
#    <Type>   : Type de mesage (MUST,ERROR,WARNING,INFO,DEBUG,EXTRA)
#    <Message>: Message a afficher
#    <Time>   : Temps specifique
#
# Retour:
#
# Remarques :
#----------------------------------------------------------------------------
log_print() {
   typeset level=${1}
   typeset msg=${2}
   typeset time=${3}

   #----- Adjust count for errors/warnings
   case "$level" in
       WARNING) LOG_WARNINGS=$((LOG_WARNINGS+1));;
       ERROR)   LOG_ERRORS=$((LOG_ERRORS+1));;
   esac

   #----- No need to continue if the level is over the threshold
   [[ ${LOG_LEVELS[$level]} -gt ${LOG_LEVELS[$LOG_LEVEL]} ]] && return 0

   #----- Format the time if it is specified
   [[ -n $time ]] && time=" $(log_time $time)"

   #----- Check for event time
   typeset datetime=""
   [[ ${LOG_TIME} -eq 1 && $level != "MUST" ]] && datetime="($(date)) "

   typeset levels="($level) "
   [[ $level = "MUST" ]] && levels=""

   #----- Build the full message
   msg="${datetime}${levels}${msg}${time}"

   #----- Print on stderr in case of error
   [[ $level = "ERROR" ]] && printf -- "$msg\n" 1>&2

   #----- Add the color if need be
   [[ $LOG_COLOR -eq 1 ]] && msg="${APP_COLORS[RESET]}${LOG_COLORS[$level]}${msg}${APP_COLORS[RESET]}"

   #----- Print the message to the log file
   if [[ -z $LOG_FILE ]]; then
      printf -- "$msg\n"
   else
      printf -- "$msg\n" >> "$LOG_FILE"
   fi
}

#----------------------------------------------------------------------------
# Nom      : <log_progress>
# Creation : Septembre 2014 - E. Legault-Ouellet - CMC/CMOE
#
# But      : Afficher le progret d'une tâche dans un format standard et
#            plus facile à parser automatiquement.
#
# Parametres  :
#    <Percent>: Le pourcentage d'avancement de la job
#    <Msg>    : Message supplémentaire
#
# Retour:
#
# Remarques :
#----------------------------------------------------------------------------
log_progress() {
   typeset percent="$1"
   typeset msg="$2"

   #----- If we need to add the time
   typeset datetime=""
   [[ ${LOG_TIME} -eq 1 ]] && datetime="($(date)) "

   #----- Build message
   msg="${datetime}(PROGRESS) [%6.2f %%] ${msg}"

   #----- Add the color if need be
   [[ $LOG_COLOR -eq 1 ]] && msg="${APP_COLORS[RESET]}${LOG_COLORS[PROGRESS]}${msg}${APP_COLORS[RESET]}"

   #----- Print the progress
   if [[ -z $LOG_FILE ]]; then
      printf -- "$msg\n" "$percent"
   else
      printf -- "$msg\n" "$percent" >> "$LOG_FILE"
   fi
}

#----------------------------------------------------------------------------
# Nom      : <log_mail>
# Creation : Octobre 2009 - J.P. Gauthier - CMC/CMOE
#
# But      : Envoyer un message par courriel.
#
# Parametres  :
#    <Subject>: Sujet du message
#    <File>   : Fichier a envoyer
#    <Address>: Adresse destinataire
#
# Retour:
#
# Remarques :
#----------------------------------------------------------------------------
log_mail() {
   typeset subject=${1}
   typeset file=${2}
   typeset address=${3}

   if [[ $address = "" ]]; then
      address=${LOG_MAILTO}
   fi

   if [[ ${address} = "" ]]; then
      return 0
   fi

   if [[ -r ${file} ]] ;then
      mail -s "${LOG_MAILTITLE} - ${subject} (${LOG_JOBID})" ${address} < ${file}
   else
      printf "$file" | mail -s "${LOG_MAILTITLE} - ${subject} (${LOG_JOBID})" ${address}
   fi
}

#----------------------------------------------------------------------------
# Nom      : <log_time>
# Creation : Mai 2010 - J.P. Gauthier - CMC/CMOE
#
# But      : Formater un temps en secondes.
#
# Parametres  :
#   <Sec>     : Secondes systeme
#
# Retour      :
#   <Date>    : Date formatee
#
# Remarques   :
#----------------------------------------------------------------------------
log_time() {
   typeset secs=$1
   typeset hours mins

   hours=$((secs / 3600))
   secs=$((secs % 3600))
   mins=$((secs / 60))
   secs=$((secs % 60))

   printf "%02d:%02d:%02d\n" $hours $mins $secs
}

#----------------------------------------------------------------------------
# Nom      : <log_hook*>
# Creation : Février 2021 - E. Legault-Ouellet - CMC/CMOE
#
# But      : Hooks called to allow scripts to hook at specific places
#
# Parametres  :
#
# Retour:
#
# Remarques :
#     The main purpose of these hooks is so that other scripts hook onto these
#     functions. It should not contain anything.
#----------------------------------------------------------------------------
log_hookstart() { :; }
log_hookend() { :; }

#----------------------------------------------------------------------------
# Nom      : <log_registererrorhandler>
# Creation : Février 2021 - E. Legault-Ouellet - CMC/CMOE
#
# But      : Register the erorr handler
#
# Parametres  :
#
# Retour:
#
# Remarques :
#----------------------------------------------------------------------------
log_registererrorhandler() {
   #----- Necessary for the trap to be functionnal across bash functions
   set -o errtrace
   trap log_errorhandler ERR
}

#----------------------------------------------------------------------------
# Nom      : <log_errorhandler>
# Creation : Février 2021 - E. Legault-Ouellet - CMC/CMOE
#
# But      : Prints an error message and a stack trace upon a trapped error
#
# Parametres  :
#
# Retour:
#
# Remarques :
#----------------------------------------------------------------------------
log_errorhandler() {
   typeset n=${#FUNCNAME[@]}
   typeset fn
   typeset msg="Trapped error at [${BASH_SOURCE[1]}:${BASH_LINENO[0]}] when executing [$BASH_COMMAND]"

   #----- Add the stacktrace
   for ((i=1; i<$n-1; i++)); do
      case "${FUNCNAME[$i]}" in
         source) fn="sourced from";;
         *)      fn="called from ${FUNCNAME[$i]}(), at";;
      esac
      msg="$msg\n   $fn ${BASH_SOURCE[$((i+1))]}:${BASH_LINENO[$i]}"
   done
   log_print ERROR "$msg"
   return 0
}
