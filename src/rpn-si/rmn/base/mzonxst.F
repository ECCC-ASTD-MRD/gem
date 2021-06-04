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
***S/P  MZONXST
#     if defined(AUTO_DOC)
*
*     Routine names  -  MZONOPR (1), MVZNXST (3), MZONFIL (5)
*                       MZONXST (2), MZONIV  (4)
*
*
*     Purpose        -  Setup/control/operate the zonal extraction
*                       package from within the CMC forecast models
*
*
*     Language       -  Fortran77
*
*
*     Author         -  Bernard Dugas, RPN (July/1992)
*
*
*     Description    - (1) The first routine is the setup and control module
*                       for the package. It (A) does the memory allocation and
*                       de-allocation, (B) initilizes the required arrays and
*                       constants, (C) controls the status transitions within
*                       the package, (D) writes the zonal information to disk.
*                       MZONOPR finds most of its control information is found 
*                       in the standard file used for accumulating the results. 
*                          The routine distinguishes between four operating 
*                       states/modes as well as an information returning mode:
*                               MOD = 0 - No initializing/allocation done
*                               MOD = 1 - Initializing/allocation done
*                               MOD = 2 - Timestep started
*                               MOD = 3 - Vertical slab started
*                       Transitions are only allowed from one line to the next 
*                       or to the previous line, except when MOD = 3, since
*                       multi-tasking may require that more one vertical slab
*                       be active at once. It is important to note that until
*                       the initializing step has been done, all of the other
*                       calls to the package return without doing a thing.
*                       Modes values of 100 or more are information returning
*                       and do not affect the package's operation per se.
*
*                    - (2) The second routine transfers information from model
*                       internal scalar arrays to the zonal extraction package
*                       arrays.
*
*                    - (3) The third routine transfers information from model
*                       internal two-dimensional vector arrays to the zonal
*                       extraction package arrays. The two components have the
*                       same sampling characteristics.
*
*                    - (4) The fourth routine modifies the number of vertical
*                       levels to be extracted in a particular vertical slab.
*                       If the new value is negative, the actual number of
*                       levels is its absolute value, and the lowest
*                       extracted level is then reproduced till the
*                       bottom is reached.
*
*                    - (5) The fifth routine modifies the default name of the
*                       standard file used for accumulating the results. The
*                       default value for this name is 'noutzon'.
*
*     Parameters     - (1) a) Control oprating mode
*                          OPR  :   Requested Transition ( +:Start  -:Finish )
*                                 ( If |EQ|    1 ) Integration
*                                 ( If |EQ|    2 ) Time Step
*                                 ( If |EQ|    3 ) Vertical Slab
*                                 ( If  EQ     4 ) Reset KA      (when MOD=3)
*                                 ( If  EQ     5 ) Reset MaxTask (when MOD=0)
*                                 ( If  EQ     6 ) Reset NDELTAT (when MOD=1)
*                                 ( If  EQ     7 ) Reset StartI  (when MOD=0)
*                                 ( If  EQ     8 ) Reset StartJ  (when MOD=0)
*                                 ( If  EQ     9 ) Toggle on/off messages (when MOD=0)
*                          VAL  :   Input value associated with OPR transition
*                                 ( If |OPR| = 1 ) I/O Unit Number IUN
*                                 ( If |OPR| = 2 ) Time Step KOUNT
*                                 ( If |OPR| = 3 ) Vertical Slab J
*                                 ( If  OPR  = 4 ) Number of Vertical Levels KA
*                                 ( If  OPR  = 5 ) Max number of allowed task
*                                 ( If  OPR  = 6 ) Saving interval in timesteps
*                                 ( If  OPR  = 7 ) Starting longitude index
*                                 ( If  OPR  = 8 ) Starting latitude index
*                                 ( If  OPR  = 9 ) Ignored value
*
*                          b) Informative oprating mode
*                          OPR  :   Requested Information
*                                 ( If  EQ   100 ) Return Active mode
*                                 ( If  EQ   101 ) Return current KA
*                                 ( If  EQ   102 ) Return current MaxTask
*                                 ( If  EQ   103 ) Return current NDELTAT
*                                 ( If  EQ   104 ) Return current StartI
*                                 ( If  EQ   105 ) Return current StartJ
*                          VAL  :   Output Information
*
*                    - (2) Z    :   Source Model Array
*                          NOM  :   Variable Name of Data ( NOTE: only the
*                                   first four characters are significant )
*                          J    :   Latitude Ordinal of Data
*                          NDIM :   First dimension of Z
*                          hhhh :   Ignored parameter (backward compatibility)
*                          FACF :   Multiplicative factor applied at transfer
*                                   can be a scalar, a 1D array or a 2D array
*                                   as specified by th evalue of ORD
*                          ORD  :   Specifies how much and what is Saved
*                                 ( If EQ -3 ) Everything with a 2D FACF
*                                 ( If EQ -2 ) Everything with a 1D FACF
*                                 ( If EQ -1 ) Everything
*                                 ( If EQ  0 ) Scalar broadcast of Z(1,1)
*                                 ( If EQ  K ) Vertical broadcast of Z(.,K)
*                          tttt :   Ignored parameter (backward compatibility)
*
*                    - (3) X    :   First Vector Component Source Model Array
*                          Y    :   Second Vector Component Source Model Array
*                          XNOM :   First Variable Name ( NOTE: only the first
*                                   four characters are significant )
*                          YNOM :   Second Variable Name ( NOTE: only the first
*                                   four characters are significant )
*                          J    :   Latitude Ordinal of Data
*                          NDIM :   First dimension of X and Y
*                          FACF :   Multiplicative factor applied at transfer
*                          ORD  :   Specifies how much and what is Saved
*                                 ( If EQ -1 ) Everything
*                                 ( If EQ  0 ) Scalar broadcast of X,Y(1,1)
*                                 ( If EQ  K ) Vertical broadcast of X,Y(.,K)
*                          tttt :   Ignored parameter (backward compatibility)
*
*                    - (4) J    :   Latitude Ordinal of Data
*                          VAL  :   Number of Vertical Levels
*
*                    - (5) NAME :   Name that supersides 'noutzon'
*
*     Error codes    - (1) 01 illegal transition requested
*                          02 request to close an I/O unit never that was openned
*                          03 unable to write incomplete somx linear accumaltor
*                          04 unable to write incomplete somx2 quadratic accumaltor
*                          05 unable to update control information in state 1
*                          06 file NOUTZON does not exist at onset
*                          07 file NOUTZON is empty 
*                          08 unable to read initial control information
*                          09 unable to read SINT rotation field
*                          10 unable to read COST rotation field
*                          11 unable to read binning information
*                          12 unable to read geographical weights
*                          13 unable to read variable list/number of levels
*                          14 unable to read previous incomplete somx
*                          15 unable to read previous incomplete somx2
*                          16 request to finish a timestep that was not started
*                          17 unable to write somx linear accumaltor
*                          18 unable to write somx2 quadratic accumaltor
*                          19 unable to update control information in state 2
*                          22 request to free a vertical slab that is not allocated
*                          23 no vertical slabs free to allocate
*                          24 attempt to set MaxTask beyond allowed maximum
*                    - (2),(3) and (4)
*                          01 Routine called while no vertical slabs open
*                          02 no vertical slabs correspond to latitude ordinal J
*                          03 no vertical levels correspond to a selected variable
*                    - (2),(3)
*                          04 NDIM,NI and StartI inconsistencies
*                          05 illegal ORD processing option.
*
*
* Last Revision  -   032 B.Dugas (mai00) - redefine NIP and NBIP
*                                        - apply limits to intrants in
*                                          calculations of SOMX and SOMX2
*
#     endif
#     if !defined(NO_SOURCE)
      Subroutine mzonopr ( OPR,VAL )

***
* Historique des revisions...
*
*revision 001 G.Pellerin (oct92) ecrit les tables d'informations
*                                dans fichier standard NOUTZON
*revision 002 G.Pellerin (oct92) common zonetab pour nec
*         003 G.Pellerin (oct92) test iluk lt 0 apres lecture
*         004 G.Pellerin (nov92) test transit pour +3 a -3
*         005 G.Pellerin (nov92) variable COMPLET dans common
*         006 G.Pellerin (nov92) at kount=1 SavPas = vrai modulo()   
*         007 G.Pellerin/B.Dugas (jan93) mvznxst pour les vents
*         008 G.Pellerin (jan93) NIV (1 ou ka) calculer de posvar 
*         009 G.Pellerin (jan93) temp(ii,place+jj,task) transfert Z(*)
*         010 G.Pellerin (jan93) util(place+jj,task) position de temp
*         011 G.Pellerin (jan93) character*2 pour NON,XNOM,YNOM 
*         012 B.Dugas    (fev93) remplacer ",task)" par ",tranche)"
*         013 G.Pellerin (oct93) ecriture a kount=1 pour tendances 
*         014 G.Pellerin (nov93) permettre transit 3 a 3 (multi) 
*         015 G.Pellerin (nov93) activer le clef Control adequatement
*         016 G.Pellerin (fev94) conversion d'entier a characteres via
*                                r4astrg et strgr4a
*         017 G.Pellerin (fev94) augmenter Nombrec a 14 (ROT)
*         018 G.Pellerin (fev94) traitement differents des variables 
*                                de surface et de profil
*         019 B.Dugas    (jun94) - mzonopr devient le point d'entree
*                                  principal du fichier
*                                - deplacer l'initialisation des champs
*                                  somx/somx2 au niveau de controle +2
*         020 B.Dugas    (sep94) le mode de fonctionnement par defaut 
*                                est de ne rien faire
*         021 B.Dugas    (dec94) utiliser la fonction wkoffit afin de
*                                verifier qu'un fichier noutzon de type
*                                rnd existe deja
*         022 B.Dugas    (jan95) implanter un mode informatif
*         023 B.Dugas    (mar95) ajouter le point d'entree MZONIV
*         024 B.Dugas    (dec95) ne pas diviser par NDELTAT si KOUNT=0,1
*         025 M.Valin    (dec95) modifier le traitement du facteur FACF
*         026 M.Valin/B.Dugas (mar96) implanter les KA negatifs
*         027 B.Dugas    (jun96) - corriger la taille de LISTVAR a l'entree
*                                - fclos sur NOUTZON si OPR=-1
*                                - KA negatifs dans MVZNXST (oubli)
*         028 B.Dugas    (jul96) - ajouter le point d'entree MZONFIL
*                                - modifier mzonopr pour changer MaxTask
*         029 B.Dugas    (sep96) modifier mzonopr pour changer NDELTAT
*         030 G.Pellerin (Mai98) augmenter MaxTask et MaxVar (32 et 512)
*         031 ?????????? (?????) augmenter MaxVar a 1024
*         032 B.Dugas    (mai00) - redefinir NIP et NBIP
*                                - mettre des bornes dans les
*                                  calculs de SOMX et SOMX2
*         033 B.Dugas    (sep01) modifier mzonopr pour changer StartI,StartJ.
*                                MZONXST/MVZNXST tiennent compte de ces modifs
*         034 K. Winger  (nov06) - Read variable names as characters
*                                - Allow 4 character variables
*         035 B.Dugas    (jan07) eliminer les message d'erreur (en option)
***

      IMPLICIT NONE

*---------------------------------------------------------------------
***    Section DECLARATIONS

      Integer       MaxTask,      MaxVar
      Parameter   ( MaxTask = 32, MaxVar = 1024 )

      Integer       MaxVarP3
      Parameter   ( MaxVarP3 = MaxVar +3 )

      Real          GrosNb
      Parameter   ( GrosNb = 1.E+18 )

***    Declaration des parametres communs a mzonxst/mvznxst.

      Integer       J,NDIM,hhhh,ORD,tttt
      Real          FACF, FACF3D(NDIM,*)

***    Declaration des parametres specifiques a mzonxst.

      Real          Z(NDIM,*)
      Character *(*) NOM

***    Declaration des parametres specifiques a mvznxst.

      Character *(*) XNOM, YNOM
      Real          X(NDIM,*),Y(NDIM,*)

***    Declaration des parametres d'appel de mzonopr.

      Integer       OPR,VAL,NIV

***    Declaration du parametre d'appel de mzonfil.

      Character*(*) NAME

***    Nom/propiete/position des variables a traiter.

      Integer       var(MaxVarP3),NVAR
      Character*8   listvar(MaxVarP3)
      Integer       propvar(MaxVar), posvar(0:MaxVar)

      Save          listvar,propvar,posvar,NVAR

***    Declarations des variables statiques.

      Integer       NombreC
      Parameter   ( NombreC = 14)

      Logical       IOPrint
      Save          IOPrint
      Character*256 NOUTZON
      Save          NOUTZON
      Integer       IUN,KA,NCPU,StartI,StartJ,NDIMS
      Save          IUN,KA,NCPU,StartI,StartJ,NDIMS
      Integer       KAJ(MaxTask),listej(MaxTask),tabctl(NombreC)
      Save          KAJ,         listej,         tabctl

      Integer       NDELTAT,DELTAT,MODE, NI,NIP,NJ,NK,  
     +              NBIN,NBIP,SOMNK,COMPLET,LATMIN,ROT
     
***   Equivalence
***  +( NDELTAT, tabctl(1) ),( DELTAT, tabctl(2) ),
***  +( MODE,    tabctl(3) ),( NI,     tabctl(4) ),
***  +( NJ,      tabctl(5) ),( NK,     tabctl(6) ),
***  +( NBIN,    tabctl(7) ),( SOMNK,  tabctl(8) ),
***  +( COMPLET, tabctl(9) ),( LATMIN, tabctl(10)),
***  +( ROT,     tabctl(11))

      Common /zontab/ NDELTAT,DELTAT,MODE,NI,NIP,NJ,NK,
     +                NBIN,NBIP,SOMNK,COMPLET,LATMIN,ROT 

***    Declaration de la variable d'etat.

      Integer       ETAT, transit(6,0:3)
      Save          ETAT, transit

***    Parametres servant a l'allocation de memoire dynamique.

      Integer       sizWK
      save          sizWK

      Real          sint,cost,pds,somx,somx2,temp,WK
      Integer       bin,util

      Pointer     ( psint,  sint (NI  ,NJ   ) ),
     +            ( pcost,  cost (NI  ,NJ   ) ),
     +            ( pbin,   bin  (NI  ,NJ   ) ),
     +            ( ppds,   pds  (NI  ,NJ   ) ),
     +            ( psomx,  somx (NBIP,SOMNK) ),
     +            ( psomx2, somx2(NBIP,SOMNK) )

      Pointer     ( ptemp,  temp (NIP ,SOMNK,MaxTask) ), 
     +            ( putil,  util (     SOMNK,MaxTask) )

      Pointer     ( pWK,    WK   (NBIN*SOMNK) )              

      Save          psint,pcost, pbin, ppds,pWK,
     +              psomx,psomx2,ptemp,putil

***    Drapeaux logiques.

      Logical       Vrai,Faux, NOTOK, SavPas, Hard, Debug
      Save          Vrai,Faux,        SavPas, Hard, Debug

***    Declaration des fonctions fstxxx et de leurs parametres.

      Integer       ifrm,fstfrm, iecr,fstecr, iprm,fstprm,
     +              inbr,fstnbr, inf, fstinf, iluk,fstluk,
     +              ierr,fstouv, nil, exfin,  wkoffit,
     +              fnom,fclos

      External      fstnbr, fstouv, fstfrm, fstecr,  exfin,
     +              fstprm, fstluk, fstinf, qqexit,
     +              fnom,   fclos,  wkoffit,r4astrg 

      Integer       errcod
      External      hpalloc,hpdeallc

      Character     typvar*1,                nomvar*4,
     +              etiket*12,               grtyp*1

      Integer       dateo,deet,npas,         kount,
     +              ip1,ip2,ip20,ip3,        datyp,nbits,npak,
     +              ig1,ig2,ig3,ig4,         swa,lng,dltf,ubc,
     +              extra1,extra2,extra3

      Logical       rewrit

      Save          dateo,deet,npas,         typvar,etiket,grtyp,
     +              ip1,ip2,ip20,ip3,        datyp,nbits,
     +              ig1,ig2,ig3,ig4,         swa,lng,dltf,ubc

***    Variables de travail non statiques.

      Character     variab*4
      Real          echant,echantX,echantY,poids
      Integer       NIC,NJC,NKC, bini, deb,fin,lo,DIFP, 
     +              II,JJ,J0,JJ0,KK,numero,numeroX,numeroY,
     +              place,placeX,placeY,tranche,LocVAL,
     +              val11,val11X,val11Y, v

***    Declaration des valeurs initiales 

      Data          Vrai     /     .true.     /,
     +              Faux     /     .false.    /,
     +              Hard     /     .false.    /,
     +              Debug    /     .false.    /
      Data          IOPrint  /     .true.     /
      Data          listej   /  MaxTask * -1  /
      Data          IUN,ETAT / -1  ,  0  /
      Data          NCPU     /  1  /
      Data          NDIMS    /  0  /
      Data          StartI   /  0  /
      Data          StartJ   /  0  /
      Data          transit  /  5  ,  1  ,  5  ,  7  ,  8  ,  9  , 
     +                          3  ,  2  , -1  ,  6  , 999 , 999 ,
     +                          3  ,  3  , -2  ,  4  , 999 , 999 ,
     +                          2  ,  3  , -3  , 999 , 999 , 999 /
      Data          NOUTZON  / 'noutzon' /
*---------------------------------------------------------------------
*=====================================================================
*---------------------------------------------------------------------
      If (OPR.ge.100)                                          Then

***        Section mode informatif de la routine.

          If (OPR.eq.100)                                      Then
              VAL = ETAT
          Else If (OPR.eq.101)                                 Then
              VAL = KA
          Else If (OPR.eq.102)                                 Then
              VAL = NCPU
          Else If (OPR.eq.103)                                 Then
              VAL = NDELTAT
          Else If (OPR.eq.104)                                 Then
              VAL = StartI+1
          Else If (OPR.eq.105)                                 Then
              VAL = StartJ+1
          Else
              VAL = -1
          End If

          Return

      End If

************************************************************************

***    Section initialisation/controle de la routine.

      If ( .not.(OPR.eq.1  .or.
     +           OPR.eq.5  .or.
     +           OPR.eq.7  .or.
     +           OPR.eq.8)
     +     .and.
     +           ETAT.eq.0 ) Return
      
***    Verifier que l'operation demandee est permise.

      If (Debug)                                               Then
          If (abs( OPR ).eq.3)                                 Then
              Write(6,6201) OPR,VAL+StartJ,STartI+1
          Else
              Write(6,6202) OPR,VAL
          End IF
      End IF

      NOTOK = Vrai
      Do  II=2,transit( 1, ETAT ) +1 
          If (transit( II, ETAT ).eq.OPR)                      Then
              NOTOK = Faux
              GoTo 100
          End If
      End Do

  100 If (NOTOK .or. OPR.eq.999)                               Then

          If (IOPrint)                                         Then
              Write(6,6001) OPR,ETAT
              nil = exfin( 'Mzonopr', 'Erreur 1', 'NON' )
          End If

          Call qqexit( 1 )

      End If

************************************************************************
      If (abs( OPR ).eq.1)                                     Then

          If (OPR.gt.0)                                        Then
************************************************************************

***     Debut d'une tranche d'integration.

***    1) On ouvre et lit le fichier de controle.
***    2) On fait l'allocation de memoire dynamique.
***    3) On lit les champs d'accumulation, s'il le faut.

              IUN = VAL
 
              ierr = wkoffit( NOUTZON )

              If (ierr.ne.1 .and. ierr.ne.33)                  Then
                  If (Hard .or. ierr.ge.-1)                    Then
                      If (IOPrint)                             Then
                          Write(6,6006)
                          nil = exfin( 'Mzonopr', 'Erreur 6', 'NON' )
                      End If
                      Call qqexit( 6 )
                  Else
                      Return
                  End If
              End If

              ierr = fnom( IUN, NOUTZON, 'STD+RND+OLD',0 )

***      Verifier que IUN fait reference a un fichier de 
***      type valide. 

              inbr = fstnbr( IUN )
              If (inbr.le. 0)                                  Then

***      Le fichier est vide. Abort ou mode non-controle
***      selon la valeur du drapeau Hard.

                  If (Hard)                                    Then
                      If (IOPrint)                             Then
                          Write(6,6007)
                          nil = exfin( 'Mzonopr','Erreur 7','NON' )
                      End If
                      Call qqexit( 7 )
                  Else
                      return
                  End If

              End If

***      Ouvrir le fichier.

              inbr = fstouv( IUN, 'RND' )

***      Lire l'information de controle se trouvant dans IUN.

              inf  = fstinf( IUN,    NIC,NJC,NKC,
     +                               -1,'CONTROLE',
     +                               -1,-1,-1,
     +                               '+','T/' )

              iluk = fstluk( tabctl, inf, NIC,NJC,NKC )

              NDELTAT = tabctl( 1) 
              DELTAT  = tabctl( 2) 
              MODE    = tabctl( 3) 
              NI      = tabctl( 4)
              NJ      = tabctl( 5)
              NK      = tabctl( 6)
              NBIN    = tabctl( 7)
              SOMNK   = tabctl( 8)
              COMPLET = tabctl( 9)  
              LATMIN  = tabctl(10) 
              ROT     = tabctl(11) 

              NIP     = NI  +mod( NI  +1,2 )
              NBIP    = NBIN+mod( NBIN+1,2 )

***      Initialiser KA avec la valeur de NK.

              KA = NK

***      On doit ecrire l'etiket dans la table (temps initial).

              if (NIC*NJC*NKC .ne. NombreC .or.
     +            iluk        .lt. 0       )                   Then

                  If (IOPrint)                                 Then
                      Write(6,6008) NIC,NJC,NKC
                      nil = exfin( 'Mzonopr','Erreur 8','NON' )
                  End If

                  Call qqexit( 8 )

              End If

***      Lire tous les parametres du fichier standart de controle.

              iprm  = fstprm( inf,
     +                        dateo,deet,npas, NIC,NJC,NKC, 
     +                        nbits,datyp,    ip1,ip20,ip3, 
     +                        typvar,nomvar,etiket, 
     +                        grtyp,ig1,ig2,ig3,ig4, 
     +                        swa,lng,dltf,ubc,
     +                        extra1,extra2,extra3 )


***      Allouer l'espace d'une premiere serie de champs 
***         en fonction de cette information.

              Call hpalloc( psint, NI*NJ, errcod,0 )
              Call hpalloc( pcost, NI*NJ, errcod,0 )
              Call hpalloc( pbin,  NI*NJ, errcod,0 )
              Call hpalloc( ppds,  NI*NJ, errcod,0 )

***      Lire ces champs grilles. Commencer par SINT.

              inf  = fstinf( IUN,  NIC,NJC,NKC, 
     +                             -1,'CONTROLE ',
     +                             ip1,ip20,ip3,
     +                             '+','S/' )
              iluk = fstluk( sint, inf, NIC,NJC,NKC )

              if (NIC.NE.NI .or.
     +            NJC.NE.NJ .or.
     +            NKC.NE. 1 .or. iluk.lt.0)                    Then

                  If (IOPrint)                                 Then
                      Write(6,6009)
                      nil = exfin( 'Mzonopr','Erreur 9','NON' )
                  End If

                  Call qqexit( 9 )

              End If

***      Lire COST.

              inf  = fstinf( IUN,  NIC,NJC,NKC, 
     +                             -1,'CONTROLE ',
     +                             ip1,ip20,ip3,
     +                             '+','C/' )
              iluk = fstluk( cost, inf, NIC,NJC,NKC )

              if (NIC.NE.NI .or.
     +            NJC.NE.NJ .or.
     +            NKC.NE. 1 .or. iluk.lt.0)                    Then

                  If (IOPrint)                                 Then
                      Write(6,6010)
                      nil = exfin( 'Mzonopr','Erreur 10','NON' )
                  End If

                  Call qqexit( 10 )

              End If

***      Lire BIN.

              inf  = fstinf( IUN, NIC,NJC,NKC,
     +                            -1,'CONTROLE ',
     +                            ip1,ip20,ip3,
     +                            '+','B/' )
              iluk = fstluk( bin, inf, NIC,NJC,NKC )

              if (NIC.NE.NI .or.
     +            NJC.NE.NJ .or.
     +            NKC.NE. 1 .or. iluk.lt.0)                    Then

                  If (IOPrint)                                 Then
                      Write(6,6011)
                      nil = exfin( 'Mzonopr','Erreur 11','NON' )
                  End If

                  Call qqexit( 11 )

              End If

***      Lire PDS.

              inf  = fstinf( IUN, NIC,NJC,NKC,
     +                            -1,'CONTROLE ',
     +                            ip1,ip20,ip3,
     +                            '+','W/' )
              iluk = fstluk( pds, inf, NIC,NJC,NKC )

              if (NIC.NE.NI .or.
     +            NJC.NE.NJ .or.
     +            NKC.NE. 1 .or. iluk.lt.0)                    Then

                  If (IOPrint)                                 Then
                      Write(6,6012)
                      nil = exfin( 'Mzonopr','Erreur 12','NON' )
                  End If

                  Call qqexit( 12 )

              End If

***      Lire LISTVAR.

              nomvar = 'VC/'

              inf  = fstinf( IUN, NIC,NJC,NKC,
     +                            -1,'CONTROLE ',
     +                            ip1,ip20,ip3,
     +                            '+',nomvar )
              If (inf.lt.0)                                    Then

                  nomvar = 'V/'

                  inf  = fstinf( IUN, NIC,NJC,NKC,
     +                                -1,'CONTROLE ',
     +                                ip1,ip20,ip3,
     +                                '+',nomvar )

              End If

              If (inf.lt.0 
     +       .or.(nomvar.eq.'VC/' .and. NIC.gt.MaxVarP3*8)
     +       .or.(nomvar.eq.'V/'  .and. NIC.gt.MaxVarP3*4))    Then

                  If (IOPrint)                                 Then
                      Write(6,6013)
                      nil = exfin( 'Mzonopr','Erreur 13','NON' )
                  End If

                  Call qqexit( 13 )

              End If

              If (nomvar.eq.'VC/')                             Then

                  NVAR = (NIC/8)-2

                  iluk = fstluk( listvar, inf, NIC,NJC,NKC )

              Else If (nomvar.eq.'V/')                         Then

                  NVAR = (NIC/4)-3

                  iluk = fstluk( var, inf, NIC,NJC,NKC )

              End If

              If (iluk.lt.0)                                   Then

                  If (IOPrint)                                 Then
                      Write(6,6013)
                      nil = exfin( 'Mzonopr','Erreur 13','NON' )
                  End If

                  Call qqexit( 13 )

              End If

              If (nomvar.eq.'V/')                              Then

                  do  II=1,NVAR
                      call r4astrg (listvar(ii),var(ii),0,3)
                  enddo

              End If

***      Definir PROPVAR avec LISTVAR.

              Do  II=1,NVAR
                  propvar(ii) = 0
                  If (listvar(ii)(1:1) .eq. '.') propvar(ii) = 1
              End Do

***      Lire POSVAR.

              inf  = fstinf( IUN, NIC,NJC,NKC,
     +                            -1,'CONTROLE ',
     +                            ip1,ip20,ip3,
     +                            '+','P/' )
              iluk = fstluk( posvar, inf, NIC,NJC,NKC )
       
              If (NIC.gt.MaxVar .or. iluk.lt.0)                Then

                  If (IOPrint)                                 Then
                      Write(6,6033)
                      nil = exfin( 'Mzonopr','Erreur 33','NON' )
                  End If

                  Call qqexit( 33 )

              End If

***       posvar(0) = 1 par definition.

              posvar(0) = 1
       
              SOMNK = posvar(NVAR)-1
                  
***       Allouer la deuxieme serie de grilles. 

              If (mod(MODE,2).eq.1)
     +        Call hpalloc( psomx,    SOMNK*NBIP, errcod,0 )

              If (MODE.ge.2)
     +        Call hpalloc( psomx2,   SOMNK*NBIP, errcod,0 )

              Call hpalloc( ptemp,NIP*SOMNK*NCPU, errcod,0 )
              Call hpalloc( putil,    SOMNK*NCPU, errcod,0 )


***       Allouer l'espace de travail des fichier standart WK.

              sizWK = max( NBIN*SOMNK, NombreC )
              Call hpalloc( pWK, sizWK, errcod,0 )

              If (COMPLET.ne.0)                                Then

***       Aller lire lire le dernier jeux de champs.
***       accumules. Il contient alors le dernier
***       IP2 sauve sur l'unite I/O IUN.

                  place = mod( NBIN+1,2 )

                  If (mod(MODE,2).eq.1)                        Then

                      inf  = fstinf( IUN,  NIC,NJC,NKC,
     +                                     -1,'MZONXST',
     +                                     ip1,COMPLET,ip3, 
     +                                     '+','1/' )
                      iluk = fstluk( WK, inf, NIC,NJC,NKC )

                      If (inf .lt.0 .or. 
     +                    iluk.lt.0    )                       Then

                          If (IOPrint)                         Then
                              Write(6,6014) COMPLET
                              nil = exfin( 'Mzonopr','Erreur 14','NON' )
                          End If

                          Call qqexit( 14 )

                      End If

***                    S'assurer que la premiere dimension
***                    de somx est bien impaire.

                      placeX = 1
                      placeY = 1
                      Do  JJ=1,SOMNK
                          Do  II=1,NBIN
                              somx(placeX,1) = WK(placeY)
                              placeX = placeX+1
                              placeY = placeY+1
                          End Do
                          placeX = placeX+place
                      End Do

                  End If
                  If (MODE.ge.2)                               Then

                      inf  = fstinf( IUN,   NIC,NJC,NKC,
     +                                      -1,'MZONXST',
     +                                      ip1,COMPLET,ip3, 
     +                                      '+','2/' )
                      iluk = fstluk( WK, inf, NIC,NJC,NKC )

                      If (inf .lt.0 .or. 
     +                    iluk.lt.0    )                       Then

                          If (IOPrint)                         Then
                              Write(6,6015) COMPLET
                              nil = exfin( 'Mzonopr','Erreur 15','NON' )
                          End If

                          Call qqexit( 15 )

                      End If

***                    S'assurer que la premiere dimension
***                    de somx2 est bien impaire.

                      placeX = 1
                      placeY = 1
                      Do  JJ=1,SOMNK
                          Do  II=1,NBIN
                              somx2(placeX,1) = WK(placeY)
                              placeX = placeX+1
                              placeY = placeY+1
                          End Do
                          placeX = placeX+place
                      End Do

                  End If

              End If

***       Avant de proceder, mettre ETAT a jour.

              ETAT    = 1
                      
************************************************************************
          Else
************************************************************************

***       Fin d'une tranche d'integration.
***
***        1) On ecrit sur le fichier IUN ce que l'on a 
***             accumule/echantillonne sur psomx/psomx2.
***        2) On de-alloue la memoire dynamique.

              If (IUN.ne.VAL)                                  Then

                  If (IOPrint)                                 Then
                      Write(6,6002) IUN,VAL
                      nil = exfin( 'Mzonopr', 'Erreur 2', 'NON' )
                  End If

                  Call qqexit( 2 )

              End If

              If (COMPLET.ne.0)                                Then


***       Sauver le contenu des SOMX/SOMX2.

                  datyp  =   5
                  rewrit = Vrai
                  npak   = -32   
                  etiket = 'MZONXST'

                  If (mod(MODE,2).eq.1)                        Then

***                    Boucher les trous dans somx.

                      place = 1
                      Do  JJ=1,SOMNK
                          Do  II=1,NBIN
                              WK(place) = somx(II,JJ)
                              place     = place+1
                          End Do
                      End Do

                      nomvar = '1/'
                      iecr   = fstecr( WK, somx, 
     +                                 npak,IUN,dateo,deet,npas,
     +                                 NBIN,SOMNK,1, ip1,ip2,ip3,
     +                                 typvar,nomvar,etiket,grtyp,
     +                                 ig1,ig2,ig3,ig4,datyp, 
     +                                 rewrit )

                      If (iecr.lt.0)                           Then

                          If (IOPrint)                         Then
                              Write(6,6003) 
                              nil = exfin( 'Mzonopr','Erreur 3','NON' )
                          End If

                          Call qqexit( 3 )

                      End If

                  End If

                  If (MODE.gt.1)                               Then

***                    Boucher les trous dans somx2.

                      place = 1
                      Do  JJ=1,SOMNK
                          Do  II=1,NBIN
                              WK(place) = somx2(II,JJ)
                              place     = place+1
                          End Do
                      End Do

                      nomvar = '2/'
                      iecr   = fstecr( WK, somx2,
     +                                 npak,IUN,dateo,deet,npas,
     +                                 NBIN,SOMNK,1, ip1,ip2,ip3,
     +                                 typvar,nomvar,etiket,grtyp,
     +                                 ig1,ig2,ig3,ig4,datyp, 
     +                                 rewrit )

                      If (iecr.lt.0)                           Then

                          If (IOPrint)                         Then
                              Write(6,6004) 
                              nil = exfin( 'Mzonopr','Erreur 4','NON' )
                          End If

                          Call qqexit( 4 )

                      End If

                  End If

***       Mettre a jour le tableau 'T/'

                  datyp  =   2
                  npak   = -32
                  etiket = 'CONTROLE'
                  nomvar = 'T/'

                  nic    =   nombreC
                  njc    =   1

                  tabctl(1) = NDELTAT
                  tabctl(9) = COMPLET
 
                  iecr   = fstecr( tabctl, WK, 
     +                             npak,IUN,dateo,deet,npas,
     +                             nic,njc,1, ip1,ip20,ip3,
     +                             typvar,nomvar,etiket,grtyp,
     +                             ig1,ig2,ig3,ig4,datyp, 
     +                             rewrit )

                  If (iecr.lt.0)                               Then

                      If (IOPrint)                             Then
                          Write(6,6005) 
                          nil = exfin( 'Mzonopr','Erreur 5','NON' )
                      End If

                      Call qqexit( 5 )

                  End If

              End If

***       Liberer la memoire dynamique.

              Call hpdeallc( psint, errcod,1 )
              Call hpdeallc( pcost, errcod,1 )
              Call hpdeallc( pbin,  errcod,1 )
              Call hpdeallc( ppds,  errcod,1 )

              If (mod(MODE,2).eq.1)
     +        Call hpdeallc( psomx, errcod,1 )

              If (MODE.ge.2)
     +        Call hpdeallc( psomx2,errcod,1 )

              Call hpdeallc( ptemp, errcod,1 )
              Call hpdeallc( putil, errcod,1 )

              Call hpdeallc( pWK,   errcod,1 )

***       Fermer le fichier IUN.

              ifrm = fstfrm( IUN )
              ierr = fclos( IUN )

***            Avant de proceder, remettre a jour la variable ETAT.

              ETAT = 0

************************************************************************
          End If

      Else If (abs( OPR ).eq.2)                                Then

          If (OPR.gt.0)                                        Then
************************************************************************

***      Debut d'un pas de temps.

***        1) On met KOUNT a jour. 
***        2) Initialiser somx/somx2 au besoin.
***        3) Decider si on doit sauver quoi que
***               ce soit durant ce pas de temps.

              KOUNT = VAL
            
***      Mettre IP2/NPAS a jour.

              deet = DELTAT
              npas = KOUNT

              If (COMPLET.EQ.0)                                Then

***       Initialiser les accumulateurs somx et/ou somx2.

                  If (mod(MODE,2).eq.1)                        Then

                      Do  II=1,NBIP*SOMNK
                          somx(II,1) = 0.0
                      End Do

                  End If
                  If (MODE.gt.1)                               Then

                      Do  II=1,NBIP*SOMNK
                          somx2(II,1) = 0.0
                      End Do

                  End If

              End If

***      On devra sauver les variables pour ce pas de temps.

              If (mod(KOUNT,NDELTAT).eq.0) SavPas = Vrai
              If (kount .eq.1)             SavPas = Vrai

              If (COMPLET.ne.   0   .and.
     +            COMPLET.gt. KOUNT )                          Then

***      La valeur du pas de temps est plus petite que celle 
***      lue sur le fichier IUN.  Il faudra alors tenter de
***      lire les donnees du pas precedent.

                  ip2   = KOUNT-1
                  place = mod( NBIN+1,2 )

                  If (mod(MODE,2).eq.1)                        Then


                      inf  = fstinf( IUN,  NIC,NJC,NKC,
     +                                     -1,' ',
     +                                     ip1,ip2,ip3,
     +                                     '+','1/' )
                      iluk = fstluk( WK, inf, NIC,NJC,NKC )

                      If (inf .lt.0 .or. 
     +                    iluk.lt.0    )                       Then

                          If (IOPrint)                         Then
                              Write(6,6020) ip2
                              nil = exfin( 'Mzonopr','Erreur 20','NON' )
                          End If

                          Call qqexit( 20 )

                      End If

***                    S'assurer que la premiere dimension
***                    de somx est bien impaire.

                      placeX = 1
                      placeY = 1
                      Do  JJ=1,SOMNK
                          Do  II=1,NBIN
                              somx(placeX,1) = WK(placeY)
                              placeX = placeX+1
                              placeY = placeY+1
                          End Do
                          placeX = placeX+place
                      End Do

                  End If
                  If (MODE.ge.2)                               Then

                      inf  = fstinf( IUN,   NIC,NJC,NKC,
     +                                      -1,' ',
     +                                      ip1,ip2,ip3,
     +                                      '+','2/' )
                      iluk = fstluk( WK, inf, NIC,NJC,NKC )

                      If (inf .lt.0 .or. 
     +                    iluk.lt.0    )                       Then

                          If (IOPrint)                         Then
                              Write(6,6021) ip2
                              nil = exfin( 'Mzonopr','Erreur 21','NON' )
                          End If

                          Call qqexit( 21 )

                      End If

***                    S'assurer que la premiere dimension
***                    de somx2 est bien impaire.

                      placeX = 1
                      placeY = 1
                      Do  JJ=1,SOMNK
                          Do  II=1,NBIN
                              somx2(placeX,1) = WK(placeY)
                              placeX = placeX+1
                              placeY = placeY+1
                          End Do
                          placeX = placeX+place
                      End Do

                  End If

              End If

***     Mettre ip2 et COMPLET a jour.

              ip2     = npas
              COMPLET = ip2

***     Avant de proceder, remettre a jour la variable ETAT.

              ETAT = 2

************************************************************************
          Else
************************************************************************

***     Fin d'un pas de temps. Sommes-nous en mode SavPas ?

              KOUNT = VAL
*          
              If (KOUNT.ne.npas)                               Then

                  If (IOPrint)                                 Then
                      Write(6,6016) KOUNT,npas
                      nil = exfin( 'Mzonopr','Erreur 16','NON' )
                  End If

                  Call qqexit( 16 )
 
              End If   
         
              If (SavPas)                                      Then

***      Sauver ce que l'on a sur somx/somx2.

                  datyp  =   5
                  rewrit =  Faux
                  nkc    =   1
                  npak   =  -32 
                  etiket = 'MZONXST'

                  If (mod(MODE,2).eq.1)                        Then

***                    Boucher les trous dans somx.

                      place = 1
                      Do  JJ=1,SOMNK
                          Do  II=1,NBIN
                              WK(place) = somx(II,JJ)
                              place     = place+1
                          End Do
                      End Do

                      nomvar = '1/'
                      iecr   = fstecr( WK, somx,
     +                                 npak,IUN,dateo,deet,npas,
     +                                 NBIN,SOMNK,nkc, ip1,ip2,ip3,
     +                                 typvar,nomvar,etiket,
     +                                 grtyp,ig1,ig2,ig3,ig4,
     +                                 datyp,rewrit )

                      If (iecr.lt.0)                           Then

                          If (IOPrint)                         Then
                              Write(6,6017) ip2
                              nil = exfin( 'Mzonopr','Erreur 17','NON' )
                          End If

                          Call qqexit( 17 )

                      End If

                  End If
                  If (MODE.gt.1)                               Then

***                    Boucher les trous dans somx2.

                      place = 1
                      Do  JJ=1,SOMNK
                          Do  II=1,NBIN
                              WK(place) = somx2(II,JJ)
                              place     = place+1
                          End Do
                      End Do

                      nomvar = '2/'
                      iecr   = fstecr( WK, somx2,
     +                                 npak,IUN,dateo,deet,npas,
     +                                 NBIN,SOMNK,nkc, ip1,ip2,ip3,
     +                                 typvar,nomvar,etiket,
     +                                 grtyp,ig1,ig2,ig3,ig4,
     +                                 datyp,rewrit )

                      If (iecr.lt.0)                           Then

                          If (IOPrint)                         Then
                              Write(6,6018) ip2
                              nil = exfin( 'Mzonopr','Erreur 18','NON' )
                          End If

                          Call qqexit( 18 )

                      End If

                  End If

***       Mettre a jour le tableau 'T/'

                  COMPLET  = 0
                  datyp  =   2
                  rewrit = Vrai
                  npak   = -32
                  etiket = 'CONTROLE'
                  nomvar = 'T/'
                  nic    =   nombreC
                  njc    =   1
                 
                  tabctl(1) = NDELTAT
                  tabctl(9) = COMPLET 

                  iecr   = fstecr( tabctl, WK,
     +                             npak,IUN,dateo,deet,npas,
     +                             nic,njc,nkc, ip1,ip20,ip3,
     +                             typvar,nomvar,etiket,
     +                             grtyp,ig1,ig2,ig3,ig4,
     +                             datyp,rewrit )

                  If (iecr.lt.0)                           Then

                      If (IOPrint)                         Then
                          Write(6,6019) ip2
                          nil = exfin( 'Mzonopr','Erreur 19','NON' )
                      End If

                      Call qqexit( 19 )

                  End If

***       Fermer/Ouvrir le fichier IUN.

                  inbr = fstfrm( IUN )
                  inbr = fstouv( IUN , 'RND' )

***       ReMettre le drapeau SavPas a faux.

                  SavPas = Faux

              End If

***       Avant de proceder, remettre a jour la variable ETAT.

              ETAT = 1

************************************************************************
          End If

      Else If (abs( OPR ).eq.3)                                Then

          LocVAL = VAL+StartJ

          If (OPR.gt.0)                                        Then
************************************************************************

***        Debut d'une  tranche verticale.
***
***        1) Allouer une tranche temporaire dans TEMP.
***        2) Initialiser temp et util pour cette tranche.

              Do KK=1,NCPU

                  If (listej(KK).eq.-1)                        Then

                      listej(KK) = LocVAL
                      KAJ(KK)    = KA
#if defined(NEC)
*vdir novector
#endif
                      Do  JJ=1,SOMNK
                          util(JJ,KK) = -1
                      End Do

                      Do  II=1,NIP*SOMNK
                          temp(II,1,KK) = 0.0
                      End Do

***                    Avant de proceder, remettre a jour la
***                    variable ETAT.

                      ETAT = 3

                      Return

                  End If

              End Do

***      Aucune tranche disponible.

              If (IOPrint)                                     Then
                  Write(6,6023)
                  nil = exfin( 'Mzonopr','Erreur 23','NON' )
              End If

              Call qqexit( 23 )

************************************************************************
          Else
************************************************************************

***      Fin d'une tranche verticale.
***
***            1) Transferer l'information de cette tranche dans les
***               accumulateurs (i.e. temp ===> [somx/somx2]).
***            2) Liberer la tranche.

              Do  KK=1,NCPU
     
                  If (listej(KK).eq.LocVAL)                    Then

                      listej(KK) = -1

***      Determiner le debut et la fin des
***      niveaux contenant des donnees.

                      deb = 1
                      fin = SOMNK

                      Do  JJ=1,SOMNK
                          If (util(JJ,KK).ne.-1) GoTo 150
                          deb = deb+1
                      End Do
  150                 Continue

                      Do  JJ=SOMNK,1,-1
                          If (util(JJ,KK).ne.-1) GoTo 200
                          fin = fin-1
                      End Do
  200                 Continue

                      If (deb.le.fin)                          Then

***     Effectuer les transferts vers somx/somx2.

                          lo = StartI+1

                          Do  II=StartI+1,StartI+NDIMS

                              bini  = bin(lo,LocVAL)
                              poids = pds(lo,LocVAL)

                              If (mod(MODE,2).eq.1)            Then

                                  Do  JJ=deb,fin
                                      somx (bini,JJ) = somx(bini,JJ) +
     +                                temp(lo,JJ,KK) * poids
                                  End Do

                              End If

                              If (MODE.gt.1)                   Then

                                  Do  JJ=deb,fin
                                      somx2(bini,JJ) = somx2(bini,JJ) +
     +                                temp(lo,JJ,KK) * temp(lo,JJ,KK) *
     +                                                 poids
                                  End Do

                              End If

                              lo = lo+1

                          End Do

                      End If


***      Avant de proceder, remettre a jour la variable
***      ETAT (mais seulement si toutes les tranches
***             verticales sont libres).

                      Do  JJ=1,NCPU
                          If (listej(JJ).ne.-1) return
                      End Do

                      ETAT = 2
                      Return

                  End If

              End Do

***      Cette tranche ne fait pas partie de listej.

              If (IOPrint)                                     Then
                  Write(6,6022) LocVAL
                  nil = exfin( 'Mzonopr','Erreur 22','NON' )
              End If

              Call qqexit( 22 )

************************************************************************
          End If
************************************************************************
      Else If (OPR.eq.4)                                       Then

***      Redefinir KA.

          KA = VAL

************************************************************************
      Else If (OPR.eq.5)                                       Then

***      Redefinir NCPU.

          NCPU = VAL

          If (NCPU.gt.MaxTask)                                 Then
              If (IOPrint)                                     Then
                  Write(6,6024) NCPU,MaxTask
                  nil = exfin( 'Mzonopr','Erreur 24','NON' )
              End If
              Call qqexit( 24 )
          End If

************************************************************************
      Else If (OPR.eq.6)                                       Then


***      Redefinir NDELTAT.

          NDELTAT   = VAL

************************************************************************
      Else If (OPR.eq.7)                                       Then


***      Redefinir StartI.

          StartI    = VAL-1

************************************************************************
      Else If (OPR.eq.8)                                       Then


***      Redefinir StartJ.

          StartJ    = VAL-1

************************************************************************
      Else If (OPR.eq.9)                                       Then


***      Basculer la valeur de IOPrint.

          IOPrint   = .NOT.IOPrint

************************************************************************
      End If

      Return

*---------------------------------------------------------------------
*=====================================================================
*---------------------------------------------------------------------

      Entry mzonxst ( Z,NOM,J,NDIM,hhhh,FACF3D,ORD,tttt )

***    Section mise-a-jour de la routine pour les champs scalaire
***   (Valide seulement dans le mode de controle actif du package).

      If (ETAT.eq.0) Return
      
      If (ETAT.ne.3)                                           Then

          If (IOPrint)                                         Then
              Write(6,6101) ETAT
              nil = exfin( 'Mzonxst', 'Erreur 1', 'NON' )
          End If

          Call qqexit( 1 )

      End If

***    Que doit-on faire de cette variable.

      JJ = 1
      Do  II=1,NVAR
          variab = '    '
          variab(1:2) = listvar(II)(2:3)
          if ((listvar(II)(4:4) .ge. '0' .and. listvar(II)(4:4) .le. '9'
     % ) .or.
     +        (listvar(II)(4:4) .ge. 'A' .and. listvar(II)(4:4) .le. 'Z'
     % ) .or.
     +        (listvar(II)(4:4) .ge. 'a' .and. listvar(II)(4:4) .le. 'z'
     % )) then
            variab(3:3) = listvar(II)(4:4)
            if ((listvar(II)(5:5) .ge. '0' .and. listvar(II)(5:5) .le. '
     %9' ) .or.
     +          (listvar(II)(5:5) .ge. 'A' .and. listvar(II)(5:5) .le. '
     %Z' ) .or.
     +          (listvar(II)(5:5) .ge. 'a' .and. listvar(II)(5:5) .le. '
     %z' )) then
              variab(4:4) = listvar(II)(5:5)
            endif
          endif

          If (trim(variab) .eq. NOM) GoTo 1000
          JJ     = JJ+1
      End Do

 1000 numero = JJ
      If (numero.gt.NVAR)                                      Then

***        Cette variable n'est pas dans notre liste.

          Return

      Else If (.not.SavPas    .and. 
     +    propvar(numero).eq.1     )                           Then

***        C'est une variable echantillonnee
***        et on ne doit pas encore la sauver.

          Return

      Else

          J0 = J+StartJ

***     Determiner le numero de la tranche.

          JJ = 1
          Do  II=1,MaxTask
              If (J0.eq.listej(II)) GoTo 1050
              JJ = JJ+1
          End Do

 1050     tranche = JJ

          If (tranche.gt.MaxTask)                              Then

***     Mauvaise latitude.

              If (IOPrint)                                     Then
                  Write(6,6102) J0
                  nil = exfin( 'Mzonxst', 'Erreur 2', 'NON' )
              End If

              Call qqexit( 2 )

          End If

***     Determiner la place que prendront ces donnees dans temp.

          place = posvar(numero-1) -1
          DIFP  = posvar(numero)-posvar(numero-1)

          If (DIFP .eq. 0)                                     Then

              If (IOPrint)                                     Then
                  Write(6,6103) NOM
                  nil = exfin( 'Mzonxst', 'Erreur 3', 'NON' )
              End If

              Call qqexit( 3 )

          End If

***    Sauver NDIM dans NDIMS.

          If (NDIMS.ne.NDIM)                                   Then
              NDIMS = NDIM
              If (Debug) write(6,6205) NDIMS
          End If

***    Tenir compte de KA negatifs.

          If (KAJ(tranche).lt.0) NIV = DIFP
          If (KAJ(tranche).gt.0) NIV = Min( DIFP ,KAJ(tranche) )

***    Variable echantillonnee ou moyennee ?

          If (propvar(numero).eq.1)                            Then
              echant = 1.0
          Else
              If (npas.le.1)                                   Then
                  echant = FACF3D(1,1)
              Else If (npas.le.NDELTAT)                        Then
                  echant = 1.0 / Float( NDELTAT-1 )
              Else
                  echant = 1.0 / Float( NDELTAT )
              End If
          End If

***    On transfert le contenu du champs Z dans temp 
***    selon la valeur trouvee dans le parametre ORD.

          NIC = Min( NDIM, NI )

          If (NIC+StartI .gt. NI)                              Then

              If (IOPrint)                                     Then
                  Write(6,6104) StartI+1,NDIM,NI
                  nil = exfin( 'Mzonxst','Erreur 4','NON' )
              End If

              Call qqexit( 4 )

          End If

          If (ORD.eq.-1)                                       Then

***    Copier NIV niveaux. Mode normal.

              echant = echant*FACF3D(1,1)

              Do  JJ0=1,NIV
                  JJ = Min( JJ0,Abs( KAJ(tranche) ) )
                  Do  II=1,NIC

                      temp(II+StartI,place+JJ0,tranche) =
     +                min(  GrosNb,
     +                max( -GrosNb,
     +                Z(II,JJ)*echant
     +                   ) )

                  End Do
              End Do


          Else If (ORD.eq.0)                                   Then

***    Broadcast scalaire de Z(1,1).

              val11 = min(  GrosNb,
     +                max( -GrosNb,
     +                Z(1,1)*echant*FACF3D(1,1)
     +                   ) )

              Do  JJ0=1,NIV
                  Do  II=1,NIC

                      temp(II+startI,place+JJ0,tranche) = val11

                  End Do
              End Do

          Else If (ORD.gt.0 .and. ORD.le.KA)                   Then

***    Broadcast vertical de Z(II,ORD).

              echant = echant*FACF3D(1,1)

              Do  JJ0=1,NIV
                  Do  II=1,NIC

                      temp(II+StartI,place+JJ0,tranche) =
     +                min(  GrosNb,
     +                max( -GrosNb,
     +                Z(II,ORD)*echant
     +                   ) )

                  End Do
              End Do

          Else If (ORD.eq.-2)                                  Then

***    Copier NIV niveaux. Mode avec normalisation 1D 

              Do  JJ0=1,NIV
                  JJ = Min( JJ0,Abs( KAJ(tranche) ) )
                  Do  II=1,NIC

                      temp(II+StartI,place+JJ0,tranche) =
     +                min(  GrosNb,
     +                max( -GrosNb,
     +                Z(II,JJ)*echant*FACF3D(II,1)
     +                   ) )

                  End Do
              End Do

          Else If (ORD.eq.-3)                                  Then

***    Copier NIV niveaux. Mode avec normalisation 2D 

              Do  JJ0=1,NIV
                  JJ = Min( JJ0,Abs( KAJ(tranche) ) )
                  Do  II=1,NIC

                      temp(II+StartI,place+JJ0,tranche) =
     +                min(  GrosNb,
     +                max( -GrosNb,
     +                Z(II,JJ)*echant*FACF3D(II,JJ)
     +                   ) )

                  End Do
              End Do

          Else

***    Valeur non-reconnue pour ORD.

              If (IOPrint)                                     Then
                  Write(6,6105) ORD,NIV
                  nil = exfin( 'Mzonxst', 'Erreur 5', 'NON' )
              End If

              Call qqexit( 5 )

          End If

***    Signaler que ces places sont utilisees.

          Do  JJ=1,NIV
              util(place+JJ,tranche) = 0
          End Do

      End If

      Return

*---------------------------------------------------------------------
*=====================================================================
*---------------------------------------------------------------------

***    Section mise-a-jour de la routine pour les champs vectoriels
***   (Valide seulement dans le mode de controle actif du package).

      Entry MVZNXST ( X,Y,XNOM,YNOM,J,NDIM,FACF3D,ORD,tttt )

      If (ETAT.eq.0) Return
      
      If (ETAT.ne.3)                                           Then

          If (IOPrint)                                         Then
              Write(6,6101) ETAT
              nil = exfin( 'Mvznxst', 'Erreur 1', 'NON' )
          End If

          Call qqexit( 1 )

      End If

***    Que doit-on faire de ces variables.

      JJ = 1
      Do  II=1,NVAR
          variab = '    '
          variab(1:2) = listvar(II)(2:3)
          if ((listvar(II)(4:4) .ge. '0' .and. listvar(II)(4:4) .le. '9'
     % ) .or.
     +        (listvar(II)(4:4) .ge. 'A' .and. listvar(II)(4:4) .le. 'Z'
     % ) .or.
     +        (listvar(II)(4:4) .ge. 'a' .and. listvar(II)(4:4) .le. 'z'
     % )) then
            variab(3:3) = listvar(II)(4:4)
            if ((listvar(II)(5:5) .ge. '0' .and. listvar(II)(5:5) .le. '
     %9' ) .or.
     +          (listvar(II)(5:5) .ge. 'A' .and. listvar(II)(5:5) .le. '
     %Z' ) .or.
     +          (listvar(II)(5:5) .ge. 'a' .and. listvar(II)(5:5) .le. '
     %z' )) then
              variab(4:4) = listvar(II)(5:5)
            endif
          endif   

          If (trim(variab) .eq. XNOM) GoTo 2000
          JJ     = JJ+1
      End Do

 2000 numeroX = JJ

      JJ = 1
      Do  II=1,NVAR
          variab = '    '
          variab(1:2) = listvar(II)(2:3)
          if ((listvar(II)(4:4) .ge. '0' .and. listvar(II)(4:4) .le. '9'
     % ) .or.
     +        (listvar(II)(4:4) .ge. 'A' .and. listvar(II)(4:4) .le. 'Z'
     % ) .or.
     +        (listvar(II)(4:4) .ge. 'a' .and. listvar(II)(4:4) .le. 'z'
     % )) then
            variab(3:3) = listvar(II)(4:4)
            if ((listvar(II)(5:5) .ge. '0' .and. listvar(II)(5:5) .le. '
     %9' ) .or.
     +          (listvar(II)(5:5) .ge. 'A' .and. listvar(II)(5:5) .le. '
     %Z' ) .or.
     +          (listvar(II)(5:5) .ge. 'a' .and. listvar(II)(5:5) .le. '
     %z' )) then
              variab(4:4) = listvar(II)(5:5)
            endif
          endif   

          If (trim(variab) .eq. YNOM) GoTo 2050
          JJ     = JJ+1
      End Do

 2050 numeroY = JJ
      If (numeroX.gt.NVAR .or.
     +    numeroY.gt.NVAR )                                    Then

***    Au moins une de ces variables n'est pas dans notre liste.

*         If (IOPrint)                                         Then
*             Write(6,6122) XNOM,YNOM
*             nil = exfin( 'Mvznxst', 'Erreur 2', 'NON' )
*         End If 
*         Call qqexit( 2 )

          Return

      Else If (.not.SavPas      .and. 
     +    (propvar(numeroX).eq.1 .or. propvar(numeroY).eq.1 )) Then

***     Au moins une des variables est une variable echantillonnee
***     et on ne doit pas encore la sauver.
        
          Return

      Else

          J0 = J+StartJ

***        Determiner le numero de la tranche.

          JJ = 1
          Do  II=1,MaxTask
              If (J0.eq.listej(II)) GoTo 2100
              JJ = JJ+1
          End Do

 2100     tranche = JJ

          If (tranche.gt.MaxTask)                              Then

***     Mauvaise latitude.

              If (IOPrint)                                     Then
                  Write(6,6102) J0
                  nil = exfin( 'Mvznxst', 'Erreur 2', 'NON' )
              End If

              Call qqexit( 2 )

          End If

***    Sauver NDIM dans NDIMS.

          If (NDIMS.ne.NDIM)                                   Then
              NDIMS = NDIM
              If (Debug) write(6,6205) NDIM
          End If

***     Determiner la place que prendront ces donnees dans temp.

          placeX = posvar(numeroX-1) -1
          placeY = posvar(numeroY-1) -1

          DIFP   = posvar(numeroX)-posvar(numeroX-1)

***     Tenir compte de KA negatifs.

          If (KAJ(tranche).lt.0) NIV = DIFP
          If (KAJ(tranche).gt.0) NIV = Min( DIFP ,KAJ(tranche) )

***     Variables echantillonnees ou moyennees ? 

          If (propvar(numeroX).eq.1 .or.
     +        propvar(numeroY).eq.1 )                          Then

              echantX = 1.0
              echantY = 1.0

          Else

***     La moyenne temporelle est requise.

              If (npas.le.1)                                   Then

***     Mais il ne peut etre question de moyennes a ce stade.

                  echantX = 1.0
                  echantY = 1.0

              Else If (npas.le.NDELTAT)                        Then

***     Et puisque que les donnees sont aussi sauvees si npas=1 ...
 
                  echantX = 1.0 / Float( NDELTAT-1 )
                  echantY = 1.0 / Float( NDELTAT-1 )

              Else

***     Ceci est le cas normal. Yeah !

                  echantX = 1.0 / Float( NDELTAT )
                  echantY = 1.0 / Float( NDELTAT )

              End If

          End If

          NIC = Min( NDIM, NI )

          If (NIC+StartI .gt. NI)                              Then

              If (IOPrint)                                     Then
                  Write(6,6104) StartI+1,NDIM,NI
                  nil = exfin( 'Mvznxst','Erreur 4','NON' )
              End If

              Call qqexit( 4 )

          End If

***     Ajuster l'orientation des composantes vectorielles.
***     Ceci est necessaire dans le modele EFR.  Dans le cas
***     du modele SEF, les champs COST et SINT contiennent
***     respectivement des 0s et des -1s. Placer le resultat 
***     dans temp selon le mode correspondant a la valeur du
***     parametre ORD.

          If (ORD.eq.-1)                                       Then

***     Copier NIV niveaux. Mode normal.

              echantX = echantX*FACF3D(1,1)
              echantY = echantY*FACF3D(1,1)

              Do  JJ0 = 1,NIV
                  JJ = Min( JJ0,Abs( KAJ(tranche) ) )
                  Do  II= 1,NIC

                      temp(II+StartI,placeX+JJ,tranche) =
     +                - echantX*X(II,JJ) * SINT(II+StartI,J0)
     +                + echantY*Y(II,JJ) * COST(II+StartI,J0)

                      temp(II+StartI,placeY+JJ,tranche) =
     +                - echantX*X(II,JJ) * COST(II+StartI,J0)
     +                - echantY*Y(II,JJ) * SINT(II+StartI,J0)

                  End Do 
              End Do 

          Else If (ORD.eq.0)                                   Then

***     Broadcast scalaire de X,Y(1,1).

              val11X =
     +        ( - echantX*X(1,1) * SINT(1+StartI,J0)
     +          + echantY*Y(1,1) * COST(1+StartI,J0) ) * FACF3D(1,1)

              val11Y =
     +        ( - echantX*X(1,1) * COST(1+StartI,J0)
     +          - echantY*Y(1,1) * SINT(1+StartI,J0) ) * FACF3D(1,1)

              Do  JJ0=1,NIV
                  Do  II=1,NIC

                      temp(II+StartI,placeX+JJ0,tranche) = val11X
                      temp(II+StartI,placeY+JJ0,tranche) = val11Y

                  End Do
              End Do

          Else If (ORD.gt.0 .and. ORD.le.KA)                   Then

***     Broadcast vertical de X,Y(II,ORD).

              echantX = echantX*FACF3D(1,1)
              echantY = echantY*FACF3D(1,1)

              Do  II=1,NIC

                  temp(II+StartI,placeX+1,tranche) =
     +            - echantX*X(II,ORD) * SINT(II+StartI,J0)
     +            + echantY*Y(II,ORD) * COST(II+StartI,J0)

                  temp(II+StartI,placeY+1,tranche) =
     +            - echantX*X(II,ORD) * COST(II+StartI,J0)
     +            - echantY*Y(II,ORD) * SINT(II+StartI,J0)

              End Do

              Do  JJ0=2,NIV
                  Do  II=1,NIC

                      temp(II+StartI,placeX+JJ0,tranche) =
     +                temp(II+StartI,placeX+1  ,tranche)

                      temp(II+StartI,placeY+JJ0,tranche) =
     +                temp(II+StartI,placeY+1  ,tranche)

                  End Do
              End Do

          Else If (ORD.eq.-2)                                  Then

***     Copier NIV niveaux. Mode avec normalisation 1D 

              Do  JJ0=1,NIV
                  JJ = Min( JJ0,Abs( KAJ(tranche) ) )
                  Do  II=1,NIC

                      temp(II+StartI,placeX+JJ0,tranche) =
     +              ( - X(II,JJ)*echantX * SINT(II+StartI,J0)
     +                + Y(II,JJ)*echantY * COST(II+StartI,J0) )
     +                * FACF3D(II,1)

                      temp(II+StartI,placeY+JJ0,tranche) =
     +              ( - X(II,JJ)*echantX * COST(II+StartI,J0)
     +                - Y(II,JJ)*echantY * SINT(II+StartI,J0) )
     +                * FACF3D(II,1)

                  End Do
              End Do

          Else If (ORD.eq.-3)                                  Then

***     Copier NIV niveaux. Mode avec normalisation 2D 

              Do  JJ0=1,NIV
                  JJ = Min( JJ0,Abs( KAJ(tranche) ) )
                  Do  II=1,NIC

                      temp(II+StartI,placeX+JJ0,tranche) =
     +              ( - X(II,JJ)*echantX * SINT(II+StartI,J0)
     +                + Y(II,JJ)*echantY * COST(II+StartI,J0) )
     +                * FACF3D(II,JJ)

                      temp(II+StartI,placeY+JJ0,tranche) =
     +              ( - X(II,JJ)*echantX * COST(II+StartI,J0)
     +                - Y(II,JJ)*echantY * SINT(II+StartI,J0) )
     +                * FACF3D(II,JJ)

                  End Do
              End Do

          Else

***     Valeur non-reconnue pour ORD.

              If (IOPrint)                                     Then
                  Write(6,6105) ORD,KA
                  nil = exfin( 'Mvznxst', 'Erreur 5', 'NON' )
              End If

              Call qqexit( 5 )

          End If

***     Signaler que ces places sont utilisees.

          Do  JJ=1,NIV
              util(placeX+JJ,tranche) = 0
              util(placeY+JJ,tranche) = 0
          End Do

      End If

      Return

*---------------------------------------------------------------------
*=====================================================================
*---------------------------------------------------------------------

***    Section de modification du KA local KAJ (encore une fois, ceci
***    n'est valide que dans le mode de controle actif du package).

      Entry MZONIV ( J,VAL )

      If (ETAT.eq.0) Return

      If (Debug) write(6,6203) J+StartJ,VAL

      If (ETAT.ne.3)                                           Then

          If (IOPrint)                                         Then
              Write(6,6101) ETAT
              nil = exfin( 'Mzoniv', 'Erreur 1', 'NON' )
          End If

          Call qqexit( 1 )

      Else

***     Determiner le numero de la tranche.

          JJ = 1
          Do  II=1,MaxTask
              If (J+StartJ.eq.listej(II)) GoTo 3000
              JJ = JJ+1
          End Do

 3000     tranche = JJ

          If (tranche.gt.MaxTask)                              Then

***     Mauvaise latitude.

              If (IOPrint)                                     Then
                  Write(6,6102) J+StartJ
                  nil = exfin( 'Mzoniv', 'Erreur 2', 'NON' )
              End If

              Call qqexit( 2 )

          End If

      End If

***     Mettre le KAJ a jour.

      KAJ(tranche) = VAL

      Return

*---------------------------------------------------------------------
*=====================================================================
*---------------------------------------------------------------------

***    Section de modification de NOUTZON.

      Entry MZONFIL ( NAME )

      If (ETAT.ne.0) Return

      If (Debug) write(6,6204) NAME

      NOUTZON = NAME

      Return

*---------------------------------------------------------------------
*=====================================================================
*---------------------------------------------------------------------

    4 Format(A4)

 6001 Format(' Illegal transition ',I3,' from state ',I3)
 6002 Format(' User requests closing of I/O unit ',I3,
     +       ' while only unit ',I3,' should be open.')
 6003 Format(' Unable to write "1/" somx before closing down.')
 6004 Format(' Unable to write "2/" somx2 before closing down.')
 6005 Format(' Unable to write "T/" control table before closing down.')
 6006 Format(" noutzon n'est pas un fichier standard.")
 6007 Format(' File xxxxxx is empty according to fstnbr.')
 6008 Format(' Unable to read "T/" control table or',
     +       ' Incorrect control table sizes :',3I5)
 6009 Format(' Unable to read "S/" sin array or',
     +       ' mismatch between expected and real dimensions.')
 6010 Format(' Unable to read "C/" cos array or',
     +       ' mismatch between expected and real dimensions.')
 6011 Format(' Unable to read "B/" bin array or',
     +       ' mismatch between expected and real dimensions.')
 6012 Format(' Unable to read "W/" weights array or',
     +       ' mismatch between expected and real dimensions.')
 6013 Format(' Unable to read "VC/" variable list or',
     +       ' mismatch between expected and real dimensions.')
 6033 Format(' Unable to read "P/" position list or',
     +       ' mismatch between expected and real dimensions.')
 6014 Format(' Unable to read "1/" somx for KOUNT =',I5,' at startup.')
 6015 Format(' Unable to read "2/" somx2 for KOUNT =',I5,' at startup.')
 6016 Format(' Non-recognized value: KOUNT =',I5,', while npas =',I5)
 6017 Format(' Unable to write "1/" somx at KOUNT =',I5)
 6018 Format(' Unable to write "2/" somx2 at KOUNT =',I5)
 6019 Format(' Unable to write "T/" control table at KOUNT =',I5)
 6020 Format(' Unable to read "1/" somx at KOUNT =',I5)
 6021 Format(' Unable to read "2/" somx2 at KOUNT =',I5)
 6022 Format(' Unknown latitude ',I5,'; unable to desallocate.')
 6023 Format(' No more available slices at this point.')
 6024 Format(' MaxTask set to ',I2,', maximum value is ',I2)

 6101 Format(' Not in a state to execute command, STATE=',I1)
 6102 Format(' Latitude ',I3,' is not being currently processed.')
 6103 Format(' Nothing to save, NOM=',A)
 6104 Format(' StartI+NDIM-1 > NI, as StartI,NDIM,NI =',3I5)
 6105 Format(' Non-recognized value: ORD =',I5,', while KA =',I5)

 6201 Format(' MzonOPR: OPR,VAl =',2I5,',and StartI =',I5)
 6202 Format(' MzonOPR: OPR,VAl =',2I5)
 6203 Format(' MzonNIV: J,VAL =',2I5)
 6204 Format(' MzonFIL: NAME =',A)
 6205 Format(' Mzonxst/Mvznxst: NDIM =',I5)

      End
#     endif
