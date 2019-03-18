***s/p rmnlib_version
*
	subroutine rmnlib_version(version,prnt)
*
*Auteur M. Lepine - RPN - Octobre 1998
*
*Objet
*  Retourner un identificateur de version de la programmatheque RMNLIB
*  qui servira de signature pour les differents programmes qui utilisent 
*  la programmatheque.
*
	character *(*) version
	logical prnt
*
        version = 
     %  "  RMNLIB  -  VRelease:"//
     %  " 16.2"//
     %  " Linux_x86-64/intel-2016.1.156 Thu May 04 18:00:06 EDT 2018"
	if (prnt) print *,version
	return
	end
