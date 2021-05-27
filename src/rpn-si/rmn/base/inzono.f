*/* RMNLIB - Library of useful routines for C and FORTRAN programming
* * Copyright (C) 1975-2001  Division de Recherche en Prevision Numerique
* *                          Environnement Canada
* *
* * This library is free software; you can redistribute it and/or
* * modify it under the terms of the GNU Lesser General Public
* * License as published by the Free Software Foundation,
* * version 2.1 of the License.
* *
* * This library is distributed in the hope that it will be useful,
* * but WITHOUT ANY WARRANTY; without even the implied warranty of
* * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
* * Lesser General Public License for more details.
* *
* * You should have received a copy of the GNU Lesser General Public
* * License along with this library; if not, write to the
* * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
* * Boston, MA 02111-1307, USA.
* */
*** S/P INZONO - initialisation des tables d'infos pour diag zonaux
*
      SUBROUTINE INZONO (POIDS, RANG, THETA, NDELTAT, DELTAT, MODE,
     X                   DZNSRF, ZSURFAC, DZNPRF, ZPROFIL, LATMIN,
     X                   ROT, IUN, S, ETIKX, IDAYO, NI, NJ, NK)
      IMPLICIT NONE
      Integer       NI,NJ,NK,LATMIN,ROT,IUN
      Integer       NDELTAT,DELTAT,MODE

      Integer       RANG(NI*NJ),IDAYO
      Integer       DZNSRF,DZNPRF
      Integer       ZSURFAC(DZNSRF),ZPROFIL(DZNPRF)
      Real          POIDS(NI*NJ),THETA(NI*NJ),S(NK)

      Character*(*) ETIKX
*
*Auteur: G.Pellerin, CMC - Septembre 1992: Version 1.0
*
*Objet(INZONO) 
*     Information de controle sont stokees pour les extractions des 
*     diagnostics zonaux
*
*REVISION 001 - G. PELLERIN - AOU93 - AJOUT DES NIVEAUX SIGMA
*         002 - G. PELLERIN - DEC93 - ADAPTATION AU MODELE SEF
*         003 - G. PELLERIN - FEV94 - CONVERSION DE CHARACTERES 
*                              A ENTIERS VIA R4ASTRG ET STRGR4A
*         004 - G. PELLERIN - FEV94 - AUGMENTER NombreC a 14
*         005 - G. PELLERIN - FEV94 - L'ANGLE DE ROTATION (ROT)
*         006 - G. PELLERIN - FEV94 - ENLEVER DEPENDANCE ZONETAB.CDK
*
*Arguments
*   IN  POIDS   - poids relatifs des points de grille pour extraction
*   IN  RANG    - le numero de la bande pour position des accumulateurs
*   IN  THETA   - angles de rotation de la grille p/r Greenwich
*   IN  NDELTAT - nombre de pas de temps d'accumulations
*   IN  DELTAT  - nombre de seconde entre chaque pas de temps
*   IN  MODE    - sauve les moyennes, la somme des moyennes, les deux
*   IN  DZNSRF  - nombre de variables de surface
*   IN  ZSURFAC - variables de surface demandees
*   IN  DZNPRF  - nombre de variables de surface
*   IN  ZPROFIL - variables de profil demandees
*   IN  LATMIN  - plus grand cercle de latitude inscrit dans la grille
*   IN  ROT     - l'angle de rotation de l'axe des X de la grille
*   IN  IUN     - le numero de fichier standard ou l'on ecrit  
*   IN  S       - les niveaux du modele
*   IN  ETIKX   - l'etiket de l'experience 
*   IN  IDAYO   - le date stamp lu de analev
*   IN  NI      - dimension horizontale de la grille du modele
*   IN  NJ      - deuxieme dimension horizontale de la grille 
*   IN  NK      - nombre de niveaux du modele
*
**    Declaration des parametres d'appel.
  
      Integer       COMPLET
      Integer       NBIN,SOMNK

*--    Section DECLARATIONS

      Integer       MaxVar,       MaxVarP2
      Parameter   ( MaxVar = 256, MaxVarP2 = MaxVar +2 )


*--    Nom/propiete/position des variables a traiter.

      Integer       var(MaxVarP2),NVAR

      Character*4   listvar(MaxVar)
      Integer       propvar(MaxVar), posvar(0:MaxVar)

      save          listvar,propvar,posvar,NVAR

*--    Declarations des variables statiques.

      Integer       NombreC
      Parameter   ( NombreC = 14 )

      Integer       tabctl(NombreC)
      Save          tabctl

*--   Parametres servant a l'allocation de memoire dynamique.
*     pWk   - espace de travail pour fichier standard

      Real           WK
      Pointer     ( pWK, WK(NI*NJ) )

*--    Drapeaux logiques.

      Logical       Vrai,Faux,  Control, Tourne
      Save          Vrai,Faux,  Control, Tourne

*--    Declaration des fonctions fstxxx et de leurs parametres.

      Integer       ifrm,fstfrm, iecr,fstecr,
     +              inbr,fstnbr, 
     +              ierr,fstouv, nil, exfin, fnom 

      External      fstnbr, fstouv, fstfrm, fstecr,  exfin, fnom,
     +              qqexit, strgr4a

      Integer       errcod
      External      hpalloc,hpdeallc

      Character     typvar*1,etivar*4,nomvar*2,
     +              etiket*8,grtyp*1

      Integer       dateo,deet,npas,        
     +              ip1,ip2,ip20,ip3,
     +              ig1,ig2,ig3,ig4,         
     +              datyp,nbits,npak 

      Logical       rewrit

      Save          dateo,deet,npas,       
     +              ip1,ip2,ip20,ip3,       
     +              ig1,ig2,ig3,ig4,         
     +              datyp,nbits,npak 

*--    Variables de travail non statiques.

      Integer       NIC,NJC,NKC, II,JJ,KK

*--    Declaration des valeurs initiales 

      Data          Control /   .true.     /
     +               Faux   /   .false.    /
     +               Vrai   /   .true.     /

      Data          COMPLET     / 0 /
      Data          nbits,   npas      / 24,  0 /
      Data          ig1,ig2,ig3,ig4    /  4 * 0 /
      Data          ip1,ip20,ip2,ip3   /  4 * 0 /

*---------------------------------------------------------------------
*=====================================================================
*---------------------------------------------------------------------
*---------------------------------------------------------------------
         Tourne  = Vrai
         If (ROT.eq.0) Tourne = Faux
      
         dateo   = IDAYO
         deet    = deltat
         npak    = -nbits

         tabctl(1) = ndeltat
         tabctl(2) = deltat 
         tabctl(3) = mode
         tabctl(4) = ni
         tabctl(5) = nj
         tabctl(6) = nk

*--    Encode les variables de surface   

         posvar(0)=1
         do ii=1,dznsrf
            WRITE(listvar(ii),'(A4)') zsurfac(ii)
            call strgr4a (listvar(ii),var(ii),0,3)
*           ces calculs seront repris dans mzonxst
            propvar(ii) = 0
            if(listvar(ii)(1:1) .eq. '.') propvar(ii) = 1
            posvar(ii) = ii + 1
         end do

*--    Encode les variables de profil   

         kk = dznsrf
         do ii=1,dznprf
            kk=kk + 1
            WRITE(listvar(kk),'(A4)') zprofil(ii)
            call strgr4a (listvar(kk),var(kk),0,3)
*           ces calculs seront repris dans mzonxst
            propvar(kk) = 0
            if(listvar(kk)(1:1) .eq. '.') propvar(kk) = 1
            posvar(kk) = kk + ii*(nk-1)+1
         end do
     
         NVAR=dznsrf+dznprf
*                        somme des deux
                         print *,'NVAR=',nvar
         do ii=1,NVAR
*        var(II)       = listvar(ii)
         print *,listvar(II),propvar(II),posvar(II)
         end do

         SOMNK=posvar(NVAR)-1
                         print *,'SOMNK= ',SOMNK

*--   Determiner NBIN.
 
         NBIN = 0
          Do JJ=1,nj
           Do II=1,ni
            If (NBIN.lt.rang(II+(JJ-1)*ni)) NBIN = rang(II+(JJ-1)*ni)
           End Do
          End Do
                         print *,'NBIN= ',NBIN

*--    Encode la tiquette

         etivar = ETIKX(1:4)
         call strgr4a (etivar,var(nvar+1),0,3) 
         etivar = ETIKX(5:8)
         call strgr4a (etivar,var(nvar+2),0,3)

                         print *,'ETIKX= ',ETIKX

*--    Ecrire le reste des tables attn. maximum de 14   

         tabctl(7)  = NBIN 
         tabctl(8)  = SOMNK 
         tabctl(9)  = COMPLET 
         tabctl(10) = LATMIN  
         tabctl(11) = ROT
         tabctl(12) = 0
         tabctl(13) = 0
         tabctl(14) = 0
         
*---------------------------------------------------------------------
*=====================================================================
*---------------------------------------------------------------------

*--    Declarer le fichier

              ierr = fnom( IUN, 'noutzon', 'STD+RND',0 )
*
*             If (ierr.ne.0)                                   Then
*                     Write(6,6000)
*                     nil = exfin( 'Zonecri', 'Erreur 1', 'NON' )
*                     Call qqexit( 1 )
*                 Else
*                     Control = Faux
*                 End If


*--   Ouvrir le fichier.

              inbr = fstouv( IUN, ' RND ' )

*--   Allouer la memoire pour les vecteurs de travail

              Call hpalloc( pWK,      NI*NJ,    errcod,1 )


*--
*--   Definir les parametres de stockage du fichier standard 
*------------------------------------------------------------------
*--
                  datyp  =   2
                  rewrit = Faux
                  typvar = '+' 
                  etiket = 'CONTROLE'
                  grtyp  = 'X'
                  If (.not. Tourne) grtyp  = 'G'
                   
*--
*--   Ecrire l'information de controle dans IUN.
*--
                  nomvar = 'T/'
                  nic    =   nombreC
                  njc    =   1
                  nkc    =   1 
                  iecr   = fstecr( tabctl, WK,
     +                             -32 ,IUN,dateo,deet,npas,
     +                             nic,njc,nkc, ip1,ip20,ip3,
     +                             typvar,nomvar,etiket,grtyp,
     +                             ig1,ig2,ig3,ig4,datyp,
     +                             rewrit
     +                           )

                  If (iecr.lt.0)                               Then

                      Write(6,6001)
                      nil = exfin( 'Zonecri','Erreur 2','NON' )
                      Call qqexit( 2 )

                  End If
*--
*--   Ecrire LISTVAR. NV+2 est contenu dans NIC.
*---------------------------------------------------------------------
*--
                  datyp = 3
                  nomvar = 'V/'
                  nic    =  ( nvar+2 )*4
                  njc    =   1
                  nkc    =   1
                  iecr   = fstecr( var, WK,
     +                             -32 ,IUN,dateo,deet,npas,
     +                             nic,njc,nkc, ip1,ip2,ip3,
     +                             typvar,nomvar,etiket,grtyp,
     +                             ig1,ig2,ig3,ig4,datyp, 
     +                             rewrit
     +                            )

                      If (iecr.lt.0)                           Then

                          Write(6,6002) 
                          nil = exfin( 'Zonecri','Erreur 3','NON' )
                          Call qqexit( 3 )

                      End If

*--
*--   Ecrire POSVAR.    posvar(0) = 1 par definition. 
*---------------------------------------------------------------------
*--
                  datyp = 2
                  nomvar = 'P/'
                  nic    =   nvar+1
                  njc    =   1
                  nkc    =   1
                  iecr   = fstecr( posvar, WK,
     +                             -32 ,IUN,dateo,deet,npas,
     +                             nic,njc,nkc, ip1,ip2,ip3,
     +                             typvar,nomvar,etiket,grtyp,
     +                             ig1,ig2,ig3,ig4,datyp,
     +                             rewrit
     +                            )

                      If (iecr.lt.0)                           Then

                          Write(6,6003)
                          nil = exfin( 'Zonecri','Erreur 4','NON' )
                          Call qqexit( 4 )
  
                      End If

*--
*--   Ecrire les champs grilles. 
*---------------------------------------------------------------------
*--        Ecrire PDS.

                  datyp = 1
                  nomvar = 'W/'
                  nic = ni
                  njc = nj
                  nkc =  1

                  iecr   = fstecr( poids, WK,
     +                             npak,IUN,dateo,deet,npas,
     +                             nic,njc,nkc, ip1,ip20,ip3,
     +                             typvar,nomvar,etiket,grtyp,
     +                             ig1,ig2,ig3,ig4,datyp,
     +                             rewrit
     +                           )

                  If (iecr.lt.0)                               Then

                      Write(6,6007)
                      nil = exfin( 'Zonecri','Erreur 7','NON' )
                      Call qqexit( 7 )

                  End If
*--*
           If(tourne)                                          Then

               Do JJ=1,nj
                Do II=1,ni
                 poids(II+(JJ-1)*ni) = SIN(theta(II+(JJ-1)*ni))
                End Do
               End Do

           Else
   
                   Do JJ=1,nj
                      Do II=1,ni
                      poids(II+(JJ-1)*ni) = -1.
                   End Do
                   End Do
  
           Endif
*--
*--        Ecrire SINT.
*--
                  datyp = 1
                  nomvar = 'S/'
                  nic = ni
                  njc = nj
                  nkc =  1
                  iecr   = fstecr( poids, WK,
     +                             npak,IUN,dateo,deet,npas,
     +                             nic,njc,nkc, ip1,ip20,ip3,
     +                             typvar,nomvar,etiket,grtyp,
     +                             ig1,ig2,ig3,ig4,datyp,
     +                             rewrit
     +                           )

                  If (iecr.lt.0)                               Then

                      Write(6,6004)
                      nil = exfin( 'Zonecri','Erreur 4','NON' )
                      Call qqexit( 4 )

                  End If

*--*
           If(Tourne)                                          Then

                 Do JJ=1,nj
                  Do II=1,ni
                   poids(II+(JJ-1)*ni) = COS(theta(II+(JJ-1)*ni))
                  End Do
                 End Do

           Else
 
                 Do JJ=1,nj
                  Do II=1,ni
                   poids(II+(JJ-1)*ni) = 0.
                  End Do
                 End Do

           Endif
*--  
*--        Ecrire COST.
*--  

                  nomvar = 'C/'
                  nic = ni
                  njc = nj
                  nkc =  1

                  iecr   = fstecr( poids, WK,
     +                             npak,IUN,dateo,deet,npas,
     +                             nic,njc,nkc, ip1,ip20,ip3,
     +                             typvar,nomvar,etiket,grtyp,
     +                             ig1,ig2,ig3,ig4,datyp,
     +                             rewrit
     +                           )

                  If (iecr.lt.0)                               Then

                      Write(6,6005)
                      nil = exfin( 'Zonecri','Erreur 5','NON' )
                      Call qqexit( 5 )

                  End If

*--
*--        Ecrire BIN.
*--
                  nomvar = 'B/'
                  datyp  = 2
                  nic = ni
                  njc = nj
                  nkc =  1

                  iecr   = fstecr( RANG, WK,
     +                             npak,IUN,dateo,deet,npas,
     +                             nic,njc,nkc, ip1,ip20,ip3,
     +                             typvar,nomvar,etiket,grtyp,
     +                             ig1,ig2,ig3,ig4,datyp,
     +                             rewrit
     +                           )

                  If (iecr.lt.0)                               Then

                      Write(6,6006)
                      nil = exfin( 'Zonecri','Erreur 6','NON' )
                      Call qqexit( 6 )

                  End If
*--
*--         Ecrire les niveaux sigma dans IUN.
*--

                  nomvar = 'S^'
                  datyp  =   1
                  nic    =   nk
                  njc    =   1
                  nkc    =   1

                  iecr   = fstecr( S , WK,
     +                             npak,IUN,dateo,deet,npas,
     +                             nic,njc,nkc, ip1,ip20,ip3,
     +                             typvar,nomvar,etikx ,grtyp,
     +                             ig1,ig2,ig3,ig4,datyp,
     +                             rewrit
     +                           )
  
                  If (iecr.lt.0)                               Then

                      Write(6,6010)
                      nil = exfin( 'Zonecri','Erreur 10','NON' )
                      Call qqexit( 10 )

                  End If
*--
*--         Ecrire les latitudes dans IUN.
*--
            if(.not.Tourne)                               Then

                  nomvar = 'L^'
                  nic    =   nbin
                  njc    =   1
                  nkc    =   1

                  iecr   = fstecr( theta , WK,
     +                             npak,IUN,dateo,deet,npas,
     +                             nic,njc,nkc, ip1,ip20,ip3,
     +                             typvar,nomvar,etikx ,grtyp,
     +                             ig1,ig2,ig3,ig4,datyp,
     +                             rewrit
     +                           )

                  If (iecr.lt.0)                               Then

                      Write(6,6011)
                      nil = exfin( 'Zonecri','Erreur 11','NON' )
                      Call qqexit( 11 )

                  End If
           End If


*--    Liberer la memoire dynamique.

              Call hpdeallc( pWK,   errcod,1 )

*--     Fermer le fichier IUN.

              ifrm = fstfrm( IUN )

        Return    

*--*********************************************************************

*--*********************************************************************

 6000 Format(' Fnom error on file xxxxxx.')
 6001 Format(' Unable to write "T/" control table before closing down.')
 6002 Format(' Unable to write "V/" variable list before closing down.')
 6003 Format(' Unable to write "P/" position list before closing down.')
 6004 Format(' Unable to write "S/" sin array before closing down.')
 6005 Format(' Unable to write "C/" cos array before closing down.')
 6006 Format(' Unable to write "B/" bin array before closing down.')
 6007 Format(' Unable to write "W/" weights array before closing down.')
 6010 Format(' Unable to write "S^" sigma array before closing down.')
 6011 Format(' Unable to write "L^" theta array before closing down.')

      End
*** S/P INZONO2 - initialisation des tables d'infos pour diag zonaux
*
      SUBROUTINE INZONO2(POIDS, RANG, THETA, NDELTAT, DELTAT, MODE,
     X                   DZNSRF, ZSURFAC, DZNPRF, ZPROFIL, LATMIN,
     X                   ROT, IUN, S, ETIKX, IDAYO, NI, NJ, NK,
     X                   Lun_out,noutzon)
*
      IMPLICIT NONE
      Integer       NI,NJ,NK,LATMIN,ROT,IUN
      Integer       NDELTAT,DELTAT,MODE,Lun_out

      Integer       RANG(NI*NJ),IDAYO
      Integer       DZNSRF,DZNPRF
      Integer       ZSURFAC(DZNSRF),ZPROFIL(DZNPRF)
      Real          POIDS(NI*NJ),THETA(NI*NJ),S(NK)

      Character*(*) ETIKX,noutzon
*
*Auteur: G.Pellerin, CMC - Septembre 1992: Version 1.0
*
*Objet(INZONO) 
*     Information de controle sont stokees pour les extractions des 
*     diagnostics zonaux
*
*REVISION 001 - G. PELLERIN - AOU93 - AJOUT DES NIVEAUX SIGMA
*         002 - G. PELLERIN - DEC93 - ADAPTATION AU MODELE SEF
*         003 - G. PELLERIN - FEV94 - CONVERSION DE CHARACTERES 
*                              A ENTIERS VIA R4ASTRG ET STRGR4A
*         004 - G. PELLERIN - FEV94 - AUGMENTER NombreC a 14
*         005 - G. PELLERIN - FEV94 - L'ANGLE DE ROTATION (ROT)
*         006 - G. PELLERIN - FEV94 - ENLEVER DEPENDANCE ZONETAB.CDK
*         007 - B. Dugas    - sep01 - Apdapter aux fichiers std2000
*
*Arguments
*   IN  POIDS   - poids relatifs des points de grille pour extraction
*   IN  RANG    - le numero de la bande pour position des accumulateurs
*   IN  THETA   - angles de rotation de la grille p/r Greenwich
*   IN  NDELTAT - nombre de pas de temps d'accumulations
*   IN  DELTAT  - nombre de seconde entre chaque pas de temps
*   IN  MODE    - sauve les moyennes, la somme des moyennes, les deux
*   IN  DZNSRF  - nombre de variables de surface
*   IN  ZSURFAC - variables de surface demandees
*   IN  DZNPRF  - nombre de variables de surface
*   IN  ZPROFIL - variables de profil demandees
*   IN  LATMIN  - plus grand cercle de latitude inscrit dans la grille
*   IN  ROT     - l'angle de rotation de l'axe des X de la grille
*   OUT IUN     - le numero de fichier standard ou l'on ecrit  
*   IN  S       - les niveaux du modele
*   IN  ETIKX   - l'etiket de l'experience 
*   IN  IDAYO   - le date stamp lu de analev
*   IN  NI      - dimension horizontale de la grille du modele
*   IN  NJ      - deuxieme dimension horizontale de la grille 
*   IN  NK      - nombre de niveaux du modele
*   IN  LUN_OUT - standard out
*   IN  NOUTZON - full pathname of the current output file
*
**    Declaration des parametres d'appel.
  
      Integer       COMPLET
      Integer       NBIN,SOMNK

*--    Section DECLARATIONS

      Integer       MaxVar,       MaxVarP3
      Parameter   ( MaxVar = 256, MaxVarP3 = MaxVar +3 )


*--    Nom/propiete/position des variables a traiter.

      Integer       var(MaxVarP3),NVAR

      Character*4   listvar(MaxVar)
      Integer       propvar(MaxVar), posvar(0:MaxVar)

      save          listvar,propvar,posvar,NVAR

*--    Declarations des variables statiques.

      Integer       NombreC
      Parameter   ( NombreC = 14 )

      Integer       tabctl(NombreC)
      Save          tabctl

*--   Parametres servant a l'allocation de memoire dynamique.
*     pWk   - espace de travail pour fichier standard

      Real           WK
      Pointer     ( pWK, WK(NI*NJ) )

*--    Drapeaux logiques.

      Logical       Vrai,Faux,  Control, Tourne
      Save          Vrai,Faux,  Control, Tourne

*--    Declaration des fonctions fstxxx et de leurs parametres.

      Integer       ifrm,fstfrm, iecr,fstecr,
     +              inbr,fstnbr, 
     +              ierr,fstouv, nil, exfin,
     +              fnom, fclos

      External      fstnbr, fstouv, fstfrm, fstecr,  exfin,
     +              fnom, fclos, qqexit, strgr4a

      Integer       errcod
      External      hpalloc,hpdeallc

      Character     typvar*1,etivar*4,nomvar*2,
     +              etiket*12,grtyp*1

      Integer       dateo,deet,npas,        
     +              ip1,ip2,ip20,ip3,
     +              ig1,ig2,ig3,ig4,         
     +              datyp,nbits,npak 

      Logical       rewrit

      Save          dateo,deet,npas,       
     +              ip1,ip2,ip20,ip3,       
     +              ig1,ig2,ig3,ig4,         
     +              datyp,nbits,npak 

*--    Variables de travail non statiques.

      Integer       NIC,NJC,NKC, II,JJ,KK

*--    Declaration des valeurs initiales 

      Data          Control /   .true.     /
     +               Faux   /   .false.    /
     +               Vrai   /   .true.     /

      Data          COMPLET     / 0 /
      Data          nbits,   npas      / 24,  0 /
      Data          ig1,ig2,ig3,ig4    /  4 * 0 /
      Data          ip1,ip20,ip2,ip3   /  4 * 0 /

*---------------------------------------------------------------------
*=====================================================================
*---------------------------------------------------------------------
*---------------------------------------------------------------------
         Tourne  = Vrai
         If (ROT.eq.0) Tourne = Faux
      
         dateo   = IDAYO
         deet    = deltat
         npak    = -nbits

         tabctl(1) = ndeltat
         tabctl(2) = deltat 
         tabctl(3) = mode
         tabctl(4) = ni
         tabctl(5) = nj
         tabctl(6) = nk

*--    Encode les variables de surface   

         posvar(0)=1
         do ii=1,dznsrf
            WRITE(listvar(ii),'(A4)') zsurfac(ii)
            call strgr4a (listvar(ii),var(ii),0,3)
*           ces calculs seront repris dans mzonxst
            propvar(ii) = 0
            if(listvar(ii)(1:1) .eq. '.') propvar(ii) = 1
            posvar(ii) = ii + 1
         end do

*--    Encode les variables de profil   

         kk = dznsrf
         do ii=1,dznprf
            kk=kk + 1
            WRITE(listvar(kk),'(A4)') zprofil(ii)
            call strgr4a (listvar(kk),var(kk),0,3)
*           ces calculs seront repris dans mzonxst
            propvar(kk) = 0
            if(listvar(kk)(1:1) .eq. '.') propvar(kk) = 1
            posvar(kk) = kk + ii*(nk-1)+1
         end do
     
*      Somme des deux

         NVAR=dznsrf+dznprf
         If (Lun_out.gt.0) write(Lun_out,'(a,i6)') 'NVAR=',nvar
         If (Lun_out.gt.0) write(Lun_out,'(a,2i6)')
     +   (listvar(II),propvar(II),posvar(II),II=1,NVAR)


         SOMNK=posvar(NVAR)-1

         If (Lun_out.gt.0) write(Lun_out,'(a,i6)') 'SOMNK= ',SOMNK

*--   Determiner NBIN.
 
         NBIN = 0
          Do JJ=1,nj
           Do II=1,ni
            If (NBIN.lt.rang(II+(JJ-1)*ni)) NBIN = rang(II+(JJ-1)*ni)
           End Do
          End Do

         If (Lun_out.gt.0) write(Lun_out,'(a,i6)') 'NBIN= ',NBIN

*--    Encode la tiquette

         etivar = ETIKX(1:4)
         call strgr4a (etivar,var(nvar+1),0,3) 
         etivar = ETIKX(5:8)
         call strgr4a (etivar,var(nvar+2),0,3)
         etivar = ETIKX(9:12)
         call strgr4a (etivar,var(nvar+3),0,3)

         If (Lun_out.gt.0) write(Lun_out,'(a,a)') 'ETIKX= ',ETIKX

*--    Ecrire le reste des tables attn. maximum de 14   

         tabctl(7)  = NBIN 
         tabctl(8)  = SOMNK 
         tabctl(9)  = COMPLET 
         tabctl(10) = LATMIN  
         tabctl(11) = ROT
         tabctl(12) = 0
         tabctl(13) = 0
         tabctl(14) = 0
         
*---------------------------------------------------------------------
*=====================================================================
*---------------------------------------------------------------------

*--    Declarer le fichier

              IUN  = 0
              ierr = fnom( IUN, noutzon, 'STD+RND',0 )

*--   Ouvrir le fichier.

              if (ierr.ge.0)                                   Then

                  inbr = fstouv( IUN, ' RND ' )

                  if (inbr.lt.0)                               Then

                      If (Lun_out.gt.0) Write(Lun_out,6002) inbr
                      nil = exfin( 'Inzono2', 'Erreur 2', 'NON' )
                      Call qqexit( 2 )

                  End If

              else

                  If (Lun_out.gt.0) Write(Lun_out,6001) ierr,noutzon
                  nil = exfin( 'Inzono2', 'Erreur 1', 'NON' )
                  Call qqexit( 1 )

              End If


*--   Allouer la memoire pour les vecteurs de travail

              Call hpalloc( pWK,      NI*NJ,    errcod,1 )


*--
*--   Definir les parametres de stockage du fichier standard 
*------------------------------------------------------------------
*--
                  datyp  = 2
                  rewrit = Faux
                  typvar = '+' 
                  etiket = 'CONTROLE'
                  grtyp  = 'X'
                  If (.not. Tourne) grtyp  = 'G'
                   
*--
*--   Ecrire l'information de controle dans IUN.
*--
                  nomvar = 'T/'
                  nic    = nombreC
                  njc    = 1
                  nkc    = 1 
                  iecr   = fstecr( tabctl, WK,
     +                             -32 ,IUN,dateo,deet,npas,
     +                             nic,njc,nkc, ip1,ip20,ip3,
     +                             typvar,nomvar,etiket,grtyp,
     +                             ig1,ig2,ig3,ig4,datyp,
     +                             rewrit
     +                           )

                  If (iecr.lt.0)                               Then

                      If (Lun_out.gt.0) Write(Lun_out,6003) iecr
                      nil = exfin( 'Inzono2','Erreur 3','NON' )
                      Call qqexit( 3 )

                  End If
*--
*--   Ecrire LISTVAR. NV+3 est contenu dans NIC.
*---------------------------------------------------------------------
*--
                  datyp  = 3
                  nomvar = 'V/'
                  nic    = (nvar+3)*4
                  njc    = 1
                  nkc    = 1
                  iecr   = fstecr( var, WK,
     +                             -32 ,IUN,dateo,deet,npas,
     +                             nic,njc,nkc, ip1,ip2,ip3,
     +                             typvar,nomvar,etiket,grtyp,
     +                             ig1,ig2,ig3,ig4,datyp, 
     +                             rewrit
     +                            )

                      If (iecr.lt.0)                           Then

                          If (Lun_out.gt.0) Write(Lun_out,6004) iecr
                          nil = exfin( 'Inzono2','Erreur 4','NON' )
                          Call qqexit( 4 )

                      End If

*--
*--   Ecrire POSVAR.    posvar(0) = 1 par definition. 
*---------------------------------------------------------------------
*--
                  datyp  = 2
                  nomvar = 'P/'
                  nic    = nvar+1
                  njc    = 1
                  nkc    = 1
                  iecr   = fstecr( posvar, WK,
     +                             -32 ,IUN,dateo,deet,npas,
     +                             nic,njc,nkc, ip1,ip2,ip3,
     +                             typvar,nomvar,etiket,grtyp,
     +                             ig1,ig2,ig3,ig4,datyp,
     +                             rewrit
     +                            )

                      If (iecr.lt.0)                           Then

                          If (Lun_out.gt.0) Write(Lun_out,6005) iecr
                          nil = exfin( 'Inzono2','Erreur 5','NON' )
                          Call qqexit( 5 )
  
                      End If

*--
*--   Ecrire les champs grilles. 
*---------------------------------------------------------------------
*--        Ecrire PDS.

                  datyp  = 1
                  nomvar = 'W/'
                  nic    = ni
                  njc    = nj
                  nkc    = 1

                  iecr   = fstecr( poids, WK,
     +                             npak,IUN,dateo,deet,npas,
     +                             nic,njc,nkc, ip1,ip20,ip3,
     +                             typvar,nomvar,etiket,grtyp,
     +                             ig1,ig2,ig3,ig4,datyp,
     +                             rewrit
     +                           )

                  If (iecr.lt.0)                               Then

                      If (Lun_out.gt.0) Write(Lun_out,6009) iecr
                      nil = exfin( 'Inzono2','Erreur 9','NON' )
                      Call qqexit( 9 )

                  End If
*--*
           If(tourne)                                          Then

               Do JJ=1,nj
                Do II=1,ni
                 poids(II+(JJ-1)*ni) = SIN(theta(II+(JJ-1)*ni))
                End Do
               End Do

           Else
   
                   Do JJ=1,nj
                      Do II=1,ni
                      poids(II+(JJ-1)*ni) = -1.
                   End Do
                   End Do
  
           Endif
*--
*--        Ecrire SINT.
*--
                  datyp  = 1
                  nomvar = 'S/'
                  nic    = ni
                  njc    = nj
                  nkc    = 1
                  iecr   = fstecr( poids, WK,
     +                             npak,IUN,dateo,deet,npas,
     +                             nic,njc,nkc, ip1,ip20,ip3,
     +                             typvar,nomvar,etiket,grtyp,
     +                             ig1,ig2,ig3,ig4,datyp,
     +                             rewrit
     +                           )

                  If (iecr.lt.0)                               Then

                      If (Lun_out.gt.0) Write(Lun_out,6006) iecr
                      nil = exfin( 'Inzono2','Erreur 6','NON' )
                      Call qqexit( 6 )

                  End If

*--*
           If(Tourne)                                          Then

                 Do JJ=1,nj
                  Do II=1,ni
                   poids(II+(JJ-1)*ni) = COS(theta(II+(JJ-1)*ni))
                  End Do
                 End Do

           Else
 
                 Do JJ=1,nj
                  Do II=1,ni
                   poids(II+(JJ-1)*ni) = 0.
                  End Do
                 End Do

           Endif
*--  
*--        Ecrire COST.
*--  

                  nomvar = 'C/'
                  nic    = ni
                  njc    = nj
                  nkc    = 1

                  iecr   = fstecr( poids, WK,
     +                             npak,IUN,dateo,deet,npas,
     +                             nic,njc,nkc, ip1,ip20,ip3,
     +                             typvar,nomvar,etiket,grtyp,
     +                             ig1,ig2,ig3,ig4,datyp,
     +                             rewrit
     +                           )

                  If (iecr.lt.0)                               Then

                      If (Lun_out.gt.0) Write(Lun_out,6007) iecr
                      nil = exfin( 'Inzono2','Erreur 7','NON' )
                      Call qqexit( 7 )

                  End If

*--
*--        Ecrire BIN.
*--
                  nomvar = 'B/'
                  datyp  = 2
                  nic    = ni
                  njc    = nj
                  nkc    = 1

                  iecr   = fstecr( RANG, WK,
     +                             npak,IUN,dateo,deet,npas,
     +                             nic,njc,nkc, ip1,ip20,ip3,
     +                             typvar,nomvar,etiket,grtyp,
     +                             ig1,ig2,ig3,ig4,datyp,
     +                             rewrit
     +                           )

                  If (iecr.lt.0)                               Then

                      If (Lun_out.gt.0) Write(Lun_out,6008) iecr
                      nil = exfin( 'Inzono2','Erreur 8','NON' )
                      Call qqexit( 8 )

                  End If
*--
*--         Ecrire les niveaux sigma dans IUN.
*--

                  nomvar = 'S^'
                  datyp  = 1
                  nic    = nk
                  njc    = 1
                  nkc    = 1

                  iecr   = fstecr( S , WK,
     +                             npak,IUN,dateo,deet,npas,
     +                             nic,njc,nkc, ip1,ip20,ip3,
     +                             typvar,nomvar,etikx ,grtyp,
     +                             ig1,ig2,ig3,ig4,datyp,
     +                             rewrit
     +                           )
  
                  If (iecr.lt.0)                               Then

                      If (Lun_out.gt.0) Write(Lun_out,6010) iecr
                      nil = exfin( 'Inzono2','Erreur 10','NON' )
                      Call qqexit( 10 )

                  End If
*--
*--         Ecrire les latitudes dans IUN.
*--
            if(.not.Tourne)                               Then

                  nomvar = 'L^'
                  nic    = nbin
                  njc    = 1
                  nkc    = 1

                  iecr   = fstecr( theta , WK,
     +                             npak,IUN,dateo,deet,npas,
     +                             nic,njc,nkc, ip1,ip20,ip3,
     +                             typvar,nomvar,etikx ,grtyp,
     +                             ig1,ig2,ig3,ig4,datyp,
     +                             rewrit
     +                           )

                  If (iecr.lt.0)                               Then

                      If (Lun_out.gt.0) Write(Lun_out,6011) iecr
                      nil = exfin( 'Inzono2','Erreur 11','NON' )
                      Call qqexit( 11 )

                  End If
           End If


*--    Liberer la memoire dynamique.

              Call hpdeallc( pWK,   errcod,1 )

*--     Fermer le fichier IUN.

              ifrm = fstfrm( IUN )
              ierr = fclos( IUN )

        Return    

*--*********************************************************************

*--*********************************************************************

 6001 Format(' Fnom error ',I5,' on file ',A)
 6002 Format(' Fstouv error ',I5)
 6003 Format(' Unable to write "T/" control table, fstecr error =',I5)
 6004 Format(' Unable to write "V/" variable list, fstecr error =',I5)
 6005 Format(' Unable to write "P/" position list, fstecr error =',I5)
 6006 Format(' Unable to write "S/" sin array, fstecr error =',I5)
 6007 Format(' Unable to write "C/" cos array, fstecr error =',I5)
 6008 Format(' Unable to write "B/" bin array, fstecr error =',I5)
 6009 Format(' Unable to write "W/" weights array, fstecr error =',I5)
 6010 Format(' Unable to write "S^" sigma array, fstecr error =',I5)
 6011 Format(' Unable to write "L^" theta array, fstecr error =',I5)

      End

*** S/P INZONO3 - initialisation des tables d'infos pour diag zonaux
*
      SUBROUTINE INZONO3(POIDS, RANG, THETA, NDELTAT, DELTAT, MODE,
     X                   DZNSRF, SURFAC, DZNPRF, PROFIL, LATMIN,
     X                   ROT, IUN, S, ETIKX, IDAYO, NI, NJ, NK, NBIN,
     X                   Lun_out,noutzon)
*
      IMPLICIT NONE
      Integer       NI,NJ,NK,NBIN,LATMIN,ROT,IUN
      Integer       NDELTAT,DELTAT,MODE,Lun_out

      Integer       RANG(NI*NJ),IDAYO
      Integer       DZNSRF,DZNPRF
      Character*6   SURFAC(DZNSRF),PROFIL(DZNPRF)
      Real          POIDS(NI*NJ),THETA(NI*NJ),S(NK)

      Character*(*) ETIKX,noutzon
*
*Auteur: G.Pellerin, CMC - Septembre 1992: Version 1.0
*
*Objet(INZONO) 
*     Information de controle sont stokees pour les extractions des 
*     diagnostics zonaux
*
*REVISION 001 - G. PELLERIN - AOU93 - AJOUT DES NIVEAUX SIGMA
*         002 - G. PELLERIN - DEC93 - ADAPTATION AU MODELE SEF
*         003 - G. PELLERIN - FEV94 - CONVERSION DE CHARACTERES 
*                              A ENTIERS VIA R4ASTRG ET STRGR4A
*         004 - G. PELLERIN - FEV94 - AUGMENTER NombreC a 14
*         005 - G. PELLERIN - FEV94 - L'ANGLE DE ROTATION (ROT)
*         006 - G. PELLERIN - FEV94 - ENLEVER DEPENDANCE ZONETAB.CDK
*         007 - B. Dugas    - sep01 - Apdapter aux fichiers std2000
*         008 - K. Winger   - nov06 - new routine name because of changes
*                                     in the list of input parameters
*                                   - parameter NBIN added for LAM mode
*                                   - Write variable names as characters
*                                     instead of r4a
*                                   - Allow 4 character variable names 
*
*Arguments
*   IN  POIDS   - poids relatifs des points de grille pour extraction
*   IN  RANG    - le numero de la bande pour position des accumulateurs
*   IN  THETA   - angles de rotation de la grille p/r Greenwich
*   IN  NDELTAT - nombre de pas de temps d'accumulations
*   IN  DELTAT  - nombre de seconde entre chaque pas de temps
*   IN  MODE    - sauve les moyennes, la somme des moyennes, les deux
*   IN  DZNSRF  - nombre de variables de surface
*   IN  SURFAC  - variables de surface demandees
*   IN  DZNPRF  - nombre de variables de surface
*   IN  PROFIL  - variables de profil demandees
*   IN  LATMIN  - plus grand cercle de latitude inscrit dans la grille
*   IN  ROT     - l'angle de rotation de l'axe des X de la grille
*   OUT IUN     - le numero de fichier standard ou l'on ecrit  
*   IN  S       - les niveaux du modele
*   IN  ETIKX   - l'etiket de l'experience 
*   IN  IDAYO   - le date stamp lu de analev
*   IN  NI      - dimension horizontale de la grille du modele
*   IN  NJ      - deuxieme dimension horizontale de la grille 
*   IN  NK      - nombre de niveaux du modele
*   IN  NBIN    - nombre de bandes zonales d'un pole a l'autre
*   IN  LUN_OUT - standard out
*   IN  NOUTZON - full pathname of the current output file
*
**    Declaration des parametres d'appel.
  
      Integer       COMPLET
      Integer       SOMNK

*--    Section DECLARATIONS

      Integer       MaxVar,       MaxVarP2
      Parameter   ( MaxVar = 256, MaxVarP2 = MaxVar +2 )


*--    Nom/propiete/position des variables a traiter.

      Character     Chaine*10000
      Character*8   listvar(MaxVarP2)
      Integer       propvar(MaxVarP2), posvar(0:MaxVarP2), NVAR

      save          listvar,propvar,posvar,NVAR

*--    Declarations des variables statiques.

      Integer       NombreC
      Parameter   ( NombreC = 14 )

      Integer       tabctl(NombreC)
      Save          tabctl

*--   Parametres servant a l'allocation de memoire dynamique.
*     Wk   - espace de travail pour fichier standard

      Real, dimension (:), Allocatable :: WK

*--    Drapeaux logiques.

      Logical       Vrai,Faux,  Control, Tourne
      Save          Vrai,Faux,  Control, Tourne

*--    Declaration des fonctions fstxxx et de leurs parametres.

      Integer       ifrm,fstfrm, iecr,fstecr,fstecr_s,
     +              inbr,fstnbr, 
     +              ierr,fstouv, nil, exfin,
     +              fnom, fclos

      External      fstnbr, fstouv, fstfrm, fstecr,fstecr_s,
     +              exfin, fnom, fclos, qqexit

      Character     typvar*1,etivar*4,nomvar*4,
     +              etiket*12,grtyp*1

      Integer       dateo,deet,npas,        
     +              ip1,ip2,ip20,ip3,
     +              ig1,ig2,ig3,ig4,         
     +              datyp,nbits,npak 

      Logical       rewrit

      Save          dateo,deet,npas,       
     +              ip1,ip2,ip20,ip3,       
     +              ig1,ig2,ig3,ig4,         
     +              datyp,nbits,npak 

*--    Variables de travail non statiques.

      Integer       NIC,NJC,NKC, II,JJ,KK, I,I0,I1

*--    Declaration des valeurs initiales 

      Data          Control /   .true.     /
     +               Faux   /   .false.    /
     +               Vrai   /   .true.     /

      Data          COMPLET     / 0 /
      Data          nbits,   npas      / 24,  0 /
      Data          ig1,ig2,ig3,ig4    /  4 * 0 /
      Data          ip1,ip20,ip2,ip3   /  4 * 0 /

*---------------------------------------------------------------------
*=====================================================================
*---------------------------------------------------------------------
*---------------------------------------------------------------------
         Tourne  = Vrai
         If (ROT.eq.0) Tourne = Faux
      
         dateo   = IDAYO
         deet    = deltat
         npak    = -nbits

         tabctl(1) = ndeltat
         tabctl(2) = deltat 
         tabctl(3) = mode
         tabctl(4) = ni
         tabctl(5) = nj
         tabctl(6) = nk

         NVAR = DZNSRF+DZNPRF
         If (NVAR.gt.MaxVar)                                   Then

             If (Lun_out.gt.0) Write(Lun_out,6012) NVAR,MaxVar
             nil = exfin( 'InZono3', 'Erreur 12', 'NON' )
             Call qqexit( 12 )

         End If

*--    Encode les variables de surface   

         posvar(0)=1
         do ii=1,dznsrf
            listvar(ii) = '        '
            listvar(ii) = surfac(ii)
*           ces calculs seront repris dans mzonxst
            propvar(ii) = 0
            if(listvar(ii)(1:1) .eq. '.') propvar(ii) = 1
            posvar(ii) = ii + 1
         end do

*--    Encode les variables de profil   

         kk = dznsrf
         do ii=1,dznprf
            kk=kk + 1
            listvar(kk) = '        '
            listvar(kk) = profil(ii)
*           ces calculs seront repris dans mzonxst
            propvar(kk) = 0
            if(listvar(kk)(1:1) .eq. '.') propvar(kk) = 1
            posvar(kk) = kk + ii*(nk-1)+1
         end do
     
*      Somme des deux

         If (Lun_out.gt.0) write(Lun_out,'(a,i6)') 'NVAR=',nvar
         If (Lun_out.gt.0) write(Lun_out,'(a,2i6)')
     +   (listvar(II),propvar(II),posvar(II),II=1,NVAR)

         SOMNK=posvar(NVAR)-1

         If (Lun_out.gt.0) write(Lun_out,'(a,i6)') 'SOMNK= ',SOMNK

         If (Lun_out.gt.0) write(Lun_out,'(a,i6)') 'NBIN= ',NBIN

*--    Encode la etiquette

         listvar(nvar+1) = ETIKX(1:8)
         listvar(nvar+2) = ETIKX(9:12)

         If (Lun_out.gt.0) write(Lun_out,'(a,a)') 'ETIKX= ',ETIKX

*--    Ecrire le reste des tables attn. maximum de 14   

         tabctl(7)  = NBIN 
         tabctl(8)  = SOMNK 
         tabctl(9)  = COMPLET 
         tabctl(10) = LATMIN  
         tabctl(11) = ROT
         tabctl(12) = 0
         tabctl(13) = 0
         tabctl(14) = 0
         
*---------------------------------------------------------------------
*=====================================================================
*---------------------------------------------------------------------

*--    Declarer le fichier

              IUN  = 0
              ierr = fnom( IUN, noutzon, 'STD+RND',0 )

*--   Ouvrir le fichier.

              if (ierr.ge.0)                                   Then

                  inbr = fstouv( IUN, ' RND ' )

                  if (inbr.lt.0)                               Then

                      If (Lun_out.gt.0) Write(Lun_out,6002) inbr
                      nil = exfin( 'InZono3', 'Erreur 2', 'NON' )
                      Call qqexit( 2 )

                  End If

              else

                  If (Lun_out.gt.0) Write(Lun_out,6001) ierr,noutzon
                  nil = exfin( 'InZono3', 'Erreur 1', 'NON' )
                  Call qqexit( 1 )

              End If


*--   Allouer la memoire pour les vecteurs de travail

              allocate( WK(NI*NJ) )

*--
*--   Definir les parametres de stockage du fichier standard 
*------------------------------------------------------------------
*--
                  datyp  = 2
                  rewrit = Faux
                  typvar = '+' 
                  etiket = 'CONTROLE'
                  grtyp  = 'X'
                  If (.not. Tourne) grtyp  = 'G'
                   
*--
*--   Ecrire l'information de controle dans IUN.
*--
                  nomvar = 'T/'
                  nic    = nombreC
                  njc    = 1
                  nkc    = 1 
                  iecr   = fstecr( tabctl, WK,
     +                             -32 ,IUN,dateo,deet,npas,
     +                             nic,njc,nkc, ip1,ip20,ip3,
     +                             typvar,nomvar,etiket,grtyp,
     +                             ig1,ig2,ig3,ig4,datyp,
     +                             rewrit
     +                           )

                  If (iecr.lt.0)                               Then

                      If (Lun_out.gt.0) Write(Lun_out,6003) iecr
                      nil = exfin( 'InZono3','Erreur 3','NON' )
                      Call qqexit( 3 )

                  End If
*--
*--   Ecrire LISTVAR. NV+2 est contenu dans NIC.
*---------------------------------------------------------------------
*--
                  datyp  = 7
                  nomvar = 'VC/'
                  nic    = (nvar+2)*8
                  njc    = 1
                  nkc    = 1

                  I0 = 1
                  I1 = 8
                  do  II=1,nvar+2
                     Chaine(I0:I1) = listvar(II)
                     I0 = I1+1
                     I1 = i0+7
                  end do

                  iecr   = fstecr_s( trim( chaine ), WK,
     +                             -8 ,IUN,dateo,deet,npas,
     +                             nic,njc,nkc, ip1,ip2,ip3,
     +                             typvar,nomvar,etiket,grtyp,
     +                             ig1,ig2,ig3,ig4,datyp, 
     +                             rewrit
     +                            )

                  If (iecr.lt.0)                               Then

                      If (Lun_out.gt.0) Write(Lun_out,6004) iecr
                      nil = exfin( 'InZono3','Erreur 4','NON' )
                      Call qqexit( 4 )

                  End If

*--
*--   Ecrire POSVAR.    posvar(0) = 1 par definition. 
*---------------------------------------------------------------------
*--
                  datyp  = 2
                  nomvar = 'P/'
                  nic    = nvar+1
                  njc    = 1
                  nkc    = 1
                  iecr   = fstecr( posvar, WK,
     +                             -32 ,IUN,dateo,deet,npas,
     +                             nic,njc,nkc, ip1,ip2,ip3,
     +                             typvar,nomvar,etiket,grtyp,
     +                             ig1,ig2,ig3,ig4,datyp,
     +                             rewrit
     +                            )

                      If (iecr.lt.0)                           Then

                          If (Lun_out.gt.0) Write(Lun_out,6005) iecr
                          nil = exfin( 'InZono3','Erreur 5','NON' )
                          Call qqexit( 5 )
  
                      End If

*--
*--   Ecrire les champs grilles. 
*---------------------------------------------------------------------
*--        Ecrire PDS.

                  datyp  = 1
                  nomvar = 'W/'
                  nic    = ni
                  njc    = nj
                  nkc    = 1

                  iecr   = fstecr( poids, WK,
     +                             npak,IUN,dateo,deet,npas,
     +                             nic,njc,nkc, ip1,ip20,ip3,
     +                             typvar,nomvar,etiket,grtyp,
     +                             ig1,ig2,ig3,ig4,datyp,
     +                             rewrit
     +                           )

                  If (iecr.lt.0)                               Then

                      If (Lun_out.gt.0) Write(Lun_out,6009) iecr
                      nil = exfin( 'InZono3','Erreur 9','NON' )
                      Call qqexit( 9 )

                  End If
*--*
           If(tourne)                                          Then

               Do JJ=1,nj
                Do II=1,ni
                 poids(II+(JJ-1)*ni) = SIN(theta(II+(JJ-1)*ni))
                End Do
               End Do

           Else
   
                   Do JJ=1,nj
                      Do II=1,ni
                      poids(II+(JJ-1)*ni) = -1.
                   End Do
                   End Do
  
           Endif
*--
*--        Ecrire SINT.
*--
                  datyp  = 1
                  nomvar = 'S/'
                  nic    = ni
                  njc    = nj
                  nkc    = 1
                  iecr   = fstecr( poids, WK,
     +                             npak,IUN,dateo,deet,npas,
     +                             nic,njc,nkc, ip1,ip20,ip3,
     +                             typvar,nomvar,etiket,grtyp,
     +                             ig1,ig2,ig3,ig4,datyp,
     +                             rewrit
     +                           )

                  If (iecr.lt.0)                               Then

                      If (Lun_out.gt.0) Write(Lun_out,6006) iecr
                      nil = exfin( 'InZono3','Erreur 6','NON' )
                      Call qqexit( 6 )

                  End If

*--*
           If(Tourne)                                          Then

                 Do JJ=1,nj
                  Do II=1,ni
                   poids(II+(JJ-1)*ni) = COS(theta(II+(JJ-1)*ni))
                  End Do
                 End Do

           Else
 
                 Do JJ=1,nj
                  Do II=1,ni
                   poids(II+(JJ-1)*ni) = 0.
                  End Do
                 End Do

           Endif
*--  
*--        Ecrire COST.
*--  

                  nomvar = 'C/'
                  nic    = ni
                  njc    = nj
                  nkc    = 1

                  iecr   = fstecr( poids, WK,
     +                             npak,IUN,dateo,deet,npas,
     +                             nic,njc,nkc, ip1,ip20,ip3,
     +                             typvar,nomvar,etiket,grtyp,
     +                             ig1,ig2,ig3,ig4,datyp,
     +                             rewrit
     +                           )

                  If (iecr.lt.0)                               Then

                      If (Lun_out.gt.0) Write(Lun_out,6007) iecr
                      nil = exfin( 'InZono3','Erreur 7','NON' )
                      Call qqexit( 7 )

                  End If

*--
*--        Ecrire BIN.
*--
                  nomvar = 'B/'
                  datyp  = 2
                  nic    = ni
                  njc    = nj
                  nkc    = 1

                  iecr   = fstecr( RANG, WK,
     +                             npak,IUN,dateo,deet,npas,
     +                             nic,njc,nkc, ip1,ip20,ip3,
     +                             typvar,nomvar,etiket,grtyp,
     +                             ig1,ig2,ig3,ig4,datyp,
     +                             rewrit
     +                           )

                  If (iecr.lt.0)                               Then

                      If (Lun_out.gt.0) Write(Lun_out,6008) iecr
                      nil = exfin( 'InZono3','Erreur 8','NON' )
                      Call qqexit( 8 )

                  End If
*--
*--         Ecrire les niveaux sigma dans IUN.
*--

                  nomvar = 'S^'
                  datyp  = 1
                  nic    = nk
                  njc    = 1
                  nkc    = 1

                  iecr   = fstecr( S , WK,
     +                             npak,IUN,dateo,deet,npas,
     +                             nic,njc,nkc, ip1,ip20,ip3,
     +                             typvar,nomvar,etikx ,grtyp,
     +                             ig1,ig2,ig3,ig4,datyp,
     +                             rewrit
     +                           )
  
                  If (iecr.lt.0)                               Then

                      If (Lun_out.gt.0) Write(Lun_out,6010) iecr
                      nil = exfin( 'InZono3','Erreur 10','NON' )
                      Call qqexit( 10 )

                  End If
*--
*--         Ecrire les latitudes dans IUN.
*--
            if(.not.Tourne)                               Then

                  nomvar = 'L^'
                  nic    = nbin
                  njc    = 1
                  nkc    = 1

                  iecr   = fstecr( theta , WK,
     +                             npak,IUN,dateo,deet,npas,
     +                             nic,njc,nkc, ip1,ip20,ip3,
     +                             typvar,nomvar,etikx ,grtyp,
     +                             ig1,ig2,ig3,ig4,datyp,
     +                             rewrit
     +                           )

                  If (iecr.lt.0)                               Then

                      If (Lun_out.gt.0) Write(Lun_out,6011) iecr
                      nil = exfin( 'InZono3','Erreur 11','NON' )
                      Call qqexit( 11 )

                  End If
           End If


*--    Liberer la memoire dynamique.

              deallocate( WK )

*--     Fermer le fichier IUN.

              ifrm = fstfrm( IUN )
              ierr = fclos( IUN )

        Return    

*--*********************************************************************

*--*********************************************************************

 6001 Format(' Fnom error ',I5,' on file ',A)
 6002 Format(' Fstouv error ',I5)
 6003 Format(' Unable to write "T/" control table, fstecr error =',I5)
 6004 Format(' Unable to write "VC/" variable list,',
     &       ' fstecr error =',I5)
 6005 Format(' Unable to write "P/" position list, fstecr error =',I5)
 6006 Format(' Unable to write "S/" sin array, fstecr error =',I5)
 6007 Format(' Unable to write "C/" cos array, fstecr error =',I5)
 6008 Format(' Unable to write "B/" bin array, fstecr error =',I5)
 6009 Format(' Unable to write "W/" weights array, fstecr error =',I5)
 6010 Format(' Unable to write "S^" sigma array, fstecr error =',I5)
 6011 Format(' Unable to write "L^" theta array, fstecr error =',I5)
 6012 Format(' NVAR ',I4,' greater than MaxVar ',I4)

      End
