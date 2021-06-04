!**s/p rmnlib_version
!
	subroutine rmnlib_version(F_version_S,F_prnt_L)
        implicit none
!
!Auteur M. Lepine - RPN - Octobre 1998
!
!Objet
!  Retourner un identificateur de version de la programmatheque RMNLIB
!  qui servira de signature pour les differents programmes qui utilisent 
!  la programmatheque.
!
        character(len=*) F_version_S
	logical F_prnt_L
!
        F_version_S =  &

         "  RMNLIB  -  VRelease:"// &
#include "../goas_rmnlib_version.inc"
! include should have something like the following:
!        " 019.6-beta"// 
!        " Linux_x86-64/intel-2019.3.199 Thu Apr 24 2020            "

 	if (F_prnt_L) print *,F_version_S
	return
	end
