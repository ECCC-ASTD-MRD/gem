!-------------------------------------- LICENCE BEGIN ------------------------
!Environment Canada - Atmospheric Science and Technology License/Disclaimer,
!                     version 3; Last Modified: May 7, 2008.
!This is free but copyrighted software; you can use/redistribute/modify it under the terms
!of the Environment Canada - Atmospheric Science and Technology License/Disclaimer
!version 3 or (at your option) any later version that should be found at:
!http://collaboration.cmc.ec.gc.ca/science/rpn.comm/license.html
!
!This software is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
!without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!See the above mentioned License/Disclaimer for more details.
!You should have received a copy of the License/Disclaimer along with this software;
!if not, you can write to: EC-RPN COMM Group, 2121 TransCanada, suite 500, Dorval (Quebec),
!CANADA, H9P 1J3; or send e-mail to service.rpn@ec.gc.ca
!-------------------------------------- LICENCE END ---------------------------

subroutine fomichev(taux, t,oz,sh,ps, ni1,ni2,nl,nk)
   use, intrinsic :: iso_fortran_env, only: REAL64
   implicit none
!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>

   integer ni1,ni2,nl,nk
   real taux(ni2,nk),t(ni1,nk),oz(ni2,nk),sh(ni2,nk),ps(nl)

   !@Author
   !          B. Dugas (Jan 2002; Following P.A. Michelangeli)
   !@Object
   !          calculate infra-red radiation rate of cooling
   !          according to (Fomichev &t Blanchet, AO 1995)
   !@Arguments
   !          - input
   ! taux     IR rate of cooling/warming in k/s as provided by radir5
   !          - output -
   ! taux     final IR rate of cooling/warming in k/s. Modified from 40 mb upward
   !          - input -
   ! t        temperature in kelvin
   ! oz       o3 mixing ratio in kg/kg
   ! sh       sigma levels for t, q, oz
   ! ps       surface pressure (n/m**2)
   ! nl       actual number of profiles to process
   ! ni1      1st dimension of t and q
   ! ni2      maximum number of profiles to process
   ! nk       number of layers
   !          - Variables utilisees pour le code radiatif de Fomichev -
   ! hfo      heating rate due to Fomichev in erg/g/s (ou cm**2/s**3)
   ! hf          "     "    "   "    "      " K/s
   ! xo       GEM levels in log(1000/pressure level (mb))  (pressure scale height)
   ! xf       temperature levels for parameterization 'pcool'
   ! xhf      heating rate levels for parameterization 'pcool'
   ! PDR      coefficients for interpolation from xo to xf grids
   ! PCD           "        "       "          "  xhf to xo grids
   ! IDR      indexes of levels for interpolation to xf grid
   ! ICD         "    "    "     "       "         " xo grid
   ! tf       temperature on levels xf
   ! o3f      ozone mixing ratio (ppmv) on levels xf
   ! vsp      ratio between original code and Fomichev code

   !     THESE ARRAYS ARE USED IN THE MATRIX PARAMETERIZATION

   real(REAL64) :: tf(69),o3f(45),hfo(61),xf(69),xhf(61), xo(120), &
        PDR(3,69),PCD(3,120),hf,vsp
   integer IDR(69),ICD(120),ik,il,i

   equivalence (xf(9),xhf(1))

   data xf / 0.00, 0.25, 0.50, 0.75, 1.00, 1.25, 1.50, 1.75, 2.00, &
        2.25, 2.50, 2.75, 3.00, 3.25, 3.50, 3.75, 4.00, 4.25, &
        4.50, 4.75, 5.00, 5.25, 5.50, 5.75, 6.00, 6.25, 6.50, &
        6.75, 7.00, 7.25, 7.50, 7.75, 8.00, 8.25, 8.50, 8.75, &
        9.00, 9.25, 9.50, 9.75,10.00,10.25,10.50,10.75,11.00, &
        11.25,11.50,11.75,12.00,12.25,12.50,12.75,13.00,13.25, &
        13.50,13.75,14.00,14.25,14.50,14.75,15.00,15.25,15.50, &
        15.75,16.00,16.25,16.50,16.75,17.00/

   save xf

   external precl,pcool

   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   ! Insertion du code radiatif de Fomichev (Tire de MAM3)
   ! (V.I. Fomichev et J.-P. Blanchet, 1995, Atmosphere-Ocean, pp.513-531)

   ! Ce code n'est pris en compte qu'a partir de 39mb. De 39mb a 7mb il remplace
   ! graduellement le code original pour rester le seul utilise a partir de 7 mb

   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   !     MATRIX PARAMETERIZATION FOR COOLING RATES IN 15 UM CO2 AND
   !     9.6 UM O3 BAND IS USED:

   do il = 1,nl

      !     Attention, dans Fomichev, on calcule du sol vers le toit alors que
      !     dans GEM c'est le contraire. On inverse donc les index sur les niveaux.
      !     La grille thermodynamique de MAM3 est sur les 1/2 niveaux alors que pour GEM
      !     il n'y a pas de difference (xo(i) n'est pas calcule avec la moyenne des niveaux
      !     i et i+1 contrairement a MAM3).

      do i = 1, nk
         xo(i) = alog(1000./(sh(il,nk-i+1)*ps(il)*0.01))  ! ps est en Pa
      enddo

      call PRECL(xo, xf, xhf, PDR, PCD, IDR, ICD, nk)

      !     temperature&ozone interpolation to parameterizations grid
      !     convert ozone mixing ratio in ppmv

      !     oz est en kg/kg. [ppmv] = 28.964/48 *1.e6 [kg/kg] = .6034e6 [kg/kg]

      !     De 1000mb a .02mb

      do i = 1,45
         if(xf(i).le.xo(nk)) then
            ik = nk - IDR(i)
            o3f(i) = (PDR(1,i)*oz(il,ik+1)+PDR(2,i)*oz(il,ik)+ &
                 PDR(3,i)*oz(il,ik-1))*.6034e6
            tf(i) = PDR(1,i)*t(il,ik+1)+PDR(2,i)*t(il,ik)+ &
                 PDR(3,i)*t(il,ik-1)
         else
            o3f(i) = oz(il,1)*.6034e6
            tf(i) = t(il,1)
         end if
      enddo

      !     Au-dessus de 0.02mb

      do i = 46, 69
         if(xf(i).le.xo(nk)) then
            ik = nk - IDR(i)
            tf(i) = PDR(1,i)*t(il,ik+1)+PDR(2,i)*t(il,ik)+ &
                 PDR(3,i)*t(il,ik-1)
         else
            tf(i) = t(il,1)
         end if
      enddo

      !     Calculate cooling rates for input grid:

      !     Convertit le taux de refroidissement de Fomichev, qui est
      !     en cm**2/s**3 (erg/g/s), en K/s. Les valeurs positives sortant
      !     de pcool correspondent a un rechauffement alors que les valeurs
      !     negatives correspondent a un refroidissement parce que pcool
      !     calcule un taux de rechauffement (comme son nom l'indique !).

      call pcool(hfo,xf,tf,o3f)

      do i = 1, nk

         if(xo(i).le.3.25) goto 9998   ! Plus bas que 38.77mb

         if(xo(i).le.5.00) then        ! Entre 38.77mb et 6.74mb

            ik = ICD(i)
            hf = 8.6e-3/86400.*(PCD(1,i)*hfo(ik)+PCD(2,i)*hfo(ik+1)+ &
                 PCD(3,i)*hfo(ik+2))
            vsp = (xo(i)-3.25)/1.75
            taux(il,nk-i+1) = taux(il,nk-i+1)*(1.-vsp) + hf*vsp

         else                          ! Plus haut que 6.74mb

            ik = ICD(i)
            taux(il,nk-i+1)= 8.6e-3/86400.*(PCD(1,i)*hfo(ik)+ &
                 PCD(2,i)*hfo(ik+1)+PCD(3,i)*hfo(ik+2))

         end if

9998     continue

      enddo
   enddo

   return

   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

end subroutine fomichev



subroutine DETINT(X,XN,N,Z,Z1,Z2,K)
   use, intrinsic :: iso_fortran_env, only: REAL64
   implicit none
!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>

   integer N,K,L1,I
   real(REAL64) :: X,XN(N),Z,Z1,Z2 &
        ,A,A1,A2

   !@Author
   !          V. Fomichev (Octobre 1993) - First implemented in MAM3 model
   !@Revision
   ! 001      P.-A. Michelangeli (March 1998) - Adaptation to RPN physics
   !@Object
   !          to calculate the arrays of second order interpolation
   !          coefficients from XN(N) grid to X level.
   !          V. Fomichev, July 17, 1992.
   !          Called by PRECL.
   !@Arguments
   !          - Input -
   ! X        target level
   ! XN(N)    grid levels
   ! N        number of levels of XN grid
   !          - Output -
   ! K        such as XN(K) .LE. X .LT. XN(K+1)
   ! Z        interpolation coefficients (weight) for XN(K)
   ! Z1       interpolation coefficients (weight) for XN(K+1)
   ! Z2       interpolation coefficients (weight) for XN(K+2)
   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   K=1
   L1=N-1
1  I=int((K+L1)/2.)
   if(I.eq.K) goto 4
   if (X-XN(I) < 0.) then
      L1=I
   else
      K=I
   endif
   goto 1
4  continue

   A=(X-XN(K+1))/(XN(K)-XN(K+2))
   A1=((X-XN(K))/(XN(K)-XN(K+1)))*(1.+A)
   A2=((X-XN(K))/(XN(K+1)-XN(K+2)))*A

   Z=1.+A1
   Z1=-A1-A2
   Z2=A2

   return
end subroutine DETINT


subroutine PRECL(X, XR, XC, PDR, PCD, IDR, ICD, N)
   use, intrinsic :: iso_fortran_env, only: REAL64
   implicit none
!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>

   integer N,I,K
   real(REAL64) :: X(N), XR(69), XC(61)
   real(REAL64) :: PDR(3,69),PCD(3,120),Z,Z1,Z2
   integer IDR(69), ICD(120)

   !@Author
   !          V. Fomichev (Octobre 1993) - First implemented in MAM3 model
   !@Revision
   ! 001      P.-A. Michelangeli (March 1998) - Adaptation to RPN physics
   !@Object
   !          To calculate the interpolation coefficients which link
   !          temperature and cooling rate grids. Second order interpolation
   !          is used.
   !@Arguments
   !          - Input -
   ! X(N)       model level set
   ! XR(69)     temperature levels for parameterization 'pcool'.
   ! XC(61)     cooling rate levels for parameterization 'pcool'.
   ! N          number of levels of XN grid
   !          - Output -
   ! PDR(3,69)  coefficients for interpolation from X(N) to XR(69) grid
   ! IDR(69)    the indexes of levels (for X(N) ) for second order
   !            interpolation to XR(69) grid
   ! PCD(3,N)   coefficients of levels for interpolation from XC(63)
   !            to X(N)
   ! ICD(N)     indexes of levels for interpolation from XC(63) to X(N)


   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   !     to form PDR(3,69),IDR(69) arrays

   do I=1,69
      if(XR(I).le.X(1)) then
         PDR(1,I) = 1.
         PDR(2,I) = 0.
         PDR(3,I) = 0.
         IDR(I) = 1
      else
         call DETINT(XR(I),X,N,Z,Z1,Z2,K)
         PDR(1,I)=Z
         PDR(2,I)=Z1
         PDR(3,I)=Z2
         IDR(I)=K
      end if
   enddo

   !     to form the PCD(3,KX),ICD(KX) arrays

   do I=1,N
      call DETINT(X(I),XC,61,Z,Z1,Z2,K)
      PCD(1,I)=Z
      PCD(2,I)=Z1
      PCD(3,I)=Z2
      ICD(I)=K
   enddo

   return
end subroutine PRECL


subroutine PCOOL(H,X,T,O3)
   use, intrinsic :: iso_fortran_env, only: REAL64
   implicit none
!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>

   real(REAL64) :: H(61),X(69),T(69),O3(45)

   !@Author
   !          V. Fomichev (Octobre 1993) - First implemented in MAM3 model
   !@Revision
   ! 001      P.-A. Michelangeli (March 1998) - Adaptation to RPN physics
   ! 002      B. Dugas and P.-A. Michelangeli (Feb 1999) -
   !          Correct multitasking bug caused by variable "FO3"
   ! 003      B. Bilodeau (Jan 2001) - Remove calls to stkmemw
   !@Object
   !          to calculate heating rates due to both the 15 um CO2 and
   !          9.6 um O3 bands
   !@Arguments
   !       - Output -
   ! H(61)    array of heating rate (erg/g/sec) at X = 2-17, 0.25
   !       - Input -
   ! X(69)    pressure scale height array; X from 0 up to 17, step 0.25;
   ! T(69)    temperature (K) at X(69) levels
   ! O3(45)   ozone volume mixing ratio (ppmv) at X = 0-11., 0.25


   !***********************************************************************
   ! Variables utilisees
   ! -------------------

   ! AU(39,6), BU(39,6)   the parameters for atmosphere layer X < 11.5
   !                                       (for 15 um CO2 band)
   ! AO3(35,6) parameters for 9.6 um O3 band
   !           (for atmosphere layer X < 10.5)
   ! KO        parameter for selection of deactivation rate coefficient for
   !           collision CO2(010)-O.
   ! SN2(23)   N2 volume mixing ratio at X = 11.5-17., 0.25
   ! O2(23)    O2 volume mixing ratio at X = 11.5-17., 0.25
   ! O(23)     O volume mixing ratio at X = 11.5-17., 0.25
   ! CO2(23)   CO2 volume mixing ratio at X = 11.5-17., 0.25
   ! XL(23)    parameters for NLTE region (X > 11.5)
   ! IG(6)     pressure scale height distance = 0.25*IG
   ! FU(69)    LTE source functions for 15 um CO2 band at p.s.h. = 0-17, 0.
   ! FO3(45)   the same for 9.6 um O3 band at p.s.h. = 0-11, 0.25
   ! AL(23)    quantum survival probability for p.s.h. = 11.5-17, 0.25
   !                                         (for 15 um CO2 band only)

   !***********************************************************************

   real(REAL64) :: O(23),O2(23),SN2(23), &
        FU(69),AU(39,6),BU(39,6),AU1(117),AU2(117),BU1(117), &
        BU2(117),XL(23),AL(23),CO2(23),FO3(50),AO3(35,6),A31(105),A32(105)
   integer IG(6)

   equivalence (AU(1,1),AU1(1)),(AU(1,4),AU2(1)), &
        (BU(1,1),BU1(1)),(BU(1,4),BU2(1)), &
        (AO3(1,1),A31(1)),(AO3(1,4),A32(1))


   integer KO
   data KO/0/

   ! USSA-1976 model: vmr for N2, O2, O, x= 11.5-17, 0.25.

   data SN2/.7816E+00, .7814E+00, .7812E+00, .7810E+00, .7807E+00, &
        .7804E+00, .7800E+00, .7793E+00, .7785E+00, .7776E+00, &
        .7768E+00, .7763E+00, .7758E+00, .7753E+00, .7748E+00, &
        .7747E+00, .7742E+00, .7735E+00, .7717E+00, .7691E+00, &
        .7647E+00, .7585E+00, .7505E+00/

   data O2/.2101E+00, .2100E+00, .2098E+00, .2096E+00, .2094E+00, &
        .2092E+00, .2086E+00, .2078E+00, .2062E+00, .2040E+00, &
        .2009E+00, .1970E+00, .1922E+00, .1865E+00, .1796E+00, &
        .1712E+00, .1621E+00, .1521E+00, .1418E+00, .1310E+00, &
        .1206E+00, .1101E+00, .1008E+00/

   data O/.7815E-04, .1033E-03, .1982E-03, .2546E-03, .5970E-03, &
        .1060E-02, .2140E-02, .3751E-02, .6155E-02, .9374E-02, &
        .1334E-01, .1793E-01, .2351E-01, .2997E-01, .3765E-01, &
        .4655E-01, .5666E-01, .6786E-01, .8049E-01, .9453E-01, &
        .1100E+00, .1274E+00, .1452E+00/



   data CO2/4*3.3E-4,3.295E-4,3.289E-4,3.284E-4,3.279E-4,3.241E-4, &
        3.206E-4,3.175E-4,3.100E-4,2.973E-4,2.859E-4,2.760E-4, &
        2.397E-4,2.058E-4,1.759E-4,1.386E-4,1.042E-4,.7606E-4, &
        .5696E-4,.4243E-4/

   data XL/.0391,.0482,.0596,.0746,.0936,.1150,.143,.183,.222,.262, &
        .304,.423,.542,.550,.687,.845,.979,.99,1.16,1.70,1.25, &
        1.40,2.05/

   data AU1/88.40,112.3,166.1,183.3,204.0,222.6,256.7,276.3,298.6, &
        396.8,447.2,625.7,691.1,699.6,849.6,856.6,952.0,1020.,1105.,1149. &
        ,1116.,1104.,1135.,1192.,1226.,1251.,1386.,1486.,1631.,1857.,2088. &
        ,2354.,2592.,2957.,3225.,3355.,3627.,4105.,4517.,983.3,1290.,1730. &
        ,1955.,2356.,2801.,3228.,3699.,4180.,4594.,4912.,4897.,4914.,4996. &
        ,4999.,5278.,5312.,5596.,5923.,6374.,6842.,7338.,7751.,8121.,8495. &
        ,8967.,9302.,9809.,10270.,10760.,11430.,12440.,14120.,16010., &
        17560.,18140.,19500.,21940.,24260.,2034.,2398.,2684.,2937.,3089., &
        3189.,3268.,3245.,3015.,3131.,20170.,21110.,21070.,22650.,24320., &
        26020.,27730.,29240., 30670.,31650.,32310.,32700.,33300.,34230., &
        35530.,37020.,41180.,45760.,49800.,52510.,59210.,69940.,83180., &
        101100.,117300.,131100.,153800.,187700.,225000./

   data AU2/-7324.,-9321.,-10320.,-12120.,-13340.,-16170.,-17050., &
        -18680.,-18920.,-18230.,-49320.,-52480.,-54040.,-57710.,-61150., &
        -64660.,-68090.,-71290.,-74420.,-76810.,-78670.,-80010.,-81870., &
        -84750.,-88740.,-93730.,-103300.,-113600.,-122100.,-130800., &
        -145200.,-169400.,-200300.,-242200.,-280300.,-312100.,-362600., &
        -439800.,-525000.,531.8,1266.,801.3,1196.,1243.,3254.,2664.,2823., &
        2280.,582.8,14360.,15480.,16710.,18210.,19430.,20460.,21320., &
        21890.,22400.,22520.,22540.,22360.,22680.,23860.,25790.,28410., &
        33090.,37350.,40150.,44130.,49230.,58980.,71870.,89740.,106500., &
        120800.,142200.,173800.,207800.,1173.,1115.,1356.,1338.,1384., &
        894.0,1092.,1078.,876.6,806.0,860.1,877.7,958.4,996.0,1085.,1184. &
        ,1239.,1336.,1427.,1505.,1585.,1652.,1712.,1742.,1860.,1980.,2104. &
        ,2386.,2704.,3145.,3666.,4407.,5369.,6636.,7987.,9202.,11060., &
        13920.,17150./

   data BU1/1815.,2427.,3437.,4853.,5538.,6175.,7058.,7919.,9322., &
        12610.,24000.,34860.,33390.,35440.,41800.,43220.,51600.,55020., &
        59550.,61970.,63990.,66330.,66910.,67900.,68020.,70470.,75420., &
        76170.,81140.,90770.,103500.,123300.,150400.,173600.,202300., &
        235400., 278800.,335800.,397300.,23540.,31480.,40210.,48680., &
        58520.,67810.,78870.,88770.,99880.,107300.,182000.,168200., &
        161300.,150600.,141200.,137900.,128700.,123500.,118900.,121100., &
        123000.,133200.,144600.,168000.,203000.,249200.,315900., 405800., &
        526700.,672900.,854800.,1059000.,1263000.,1375000.,1514000., &
        1506000.,1505000.,1495000.,1461000.,60000.,68950.,76590.,81310., &
        83430.,84940.,84860.,83560.,79140.,73140.,518800.,507300.,439300. &
        ,434300.,432400.,432500.,441600.,461200.,490800.,541000.,613300., &
        711500.,839900.,992800.,1173000.,1373000.,1578000.,1766000., &
        1917000.,2014000., 2032000.,1966000.,1773000.,1417000.,1090000., &
        755300., 499200.,299500.,122000./

   data BU2/-157100.,-202200.,-228200.,-255100.,-281600.,-335700., &
        -344400.,-361100.,-371900.,-362700.,-1346000.,-1335000.,-1288000., &
        -1294000.,-1280000.,-1273000.,-1288000.,-1322000.,-1395000., &
        -1511000.,-1679000.,-1894000.,-2184000.,-2543000.,-2944000., &
        -3410000.,-3868000.,-4318000.,-4735000.,-5107000.,-5400000., &
        -5618000.,-5661000.,-5261000.,-4985000.,-4459000.,-4110000., &
        -3823000.,-3474000.,-1914.,15870.,-996.9,4348.,5815.,42480., &
        28630.,26970.,21240.,-15490.,342100.,336400.,335900.,367200., &
        370100.,370200.,384700.,400800.,439200.,483000.,551400.,622600., &
        718200.,832100.,920000.,1017000.,1046000.,1027000.,960500., &
        863400.,728500.,582700.,448000.,287900.,158500.,60150.,5460., &
        -31750.,-41100.,27940.,26390.,32850.,33160.,34440.,25390.,30530., &
        30960.,26910.,25100.,45650.,45540.,46610.,37620.,35280.,34540., &
        35320.,33430.,35480.,37370.,38990.,50430.,64730.,82980.,112300., &
        144200.,188100., 231300.,270600.,296200.,313700.,318800.,301300., &
        264000.,235000.,195700.,163500.,136800.,108600./

   data A31/6812.,5810.,4765.,3858.,3187.,2501.,2016.,1672.,1423., &
        1389.,1390.,1384.,1406.,1445.,1474.,1578.,1678.,1816.,2067.,2294. &
        ,2496.,2849.,3366.,3912.,4739.,5730.,6981.,8266.,9803.,11470., &
        13680.,16010.,18700.,14260.,7682.,2485.,3120.,3357.,3118.,2780., &
        2586.,2422., 2299.,2250.,1900.,1593.,1378.,1196.,1139.,1191., &
        1308.,1581.,2082.,2679.,3774.,5268.,6851.,8381.,9778.,10990., &
        11460.,11970.,12030.,11170.,9974.,7343.,4209.,-197.6,-4099., &
        -4180.,14490.,14530.,14300.,13920.,13680.,13630.,13520.,13510., &
        13430.,13630.,13840.,14080.,14390.,14730.,15220.,15820.,16640., &
        17480.,18640.,19330.,19300.,18710.,17450.,15730.,13630.,11740., &
        9204.,7041.,5583.,4396.,4163.,4532.,5945.,5556.,3694./

   data A32/-49280.,-47830.,-45050.,-41980.,-39490.,-37650.,-36240., &
        -35170.,-34390.,-34050.,-33840.,-33950.,-34280.,-35080.,-36250., &
        -38150.,-40390.,-43400.,-46520.,-49360.,-51400.,-53010.,-53570., &
        -53760.,-53490.,-52890.,-52070.,-51200.,-50390.,-49570.,-48910., &
        -48470.,-48150.,-31120.,-14300.,17320.,17120.,16260.,15390.,14620. &
        ,14050.,13590.,13200.,12900.,12780.,12670.,12690.,12660.,12900., &
        13010.,13270.,13060.,12800.,11770.,10170.,8252.,6737.,4899.,3622. &
        ,2473.,1617.,1095.,821.9,543.8,275.8,174.3,113.8,70.18,13.74, &
        3.553,3912.,3298.,2689.,2205.,1829.,1535.,1329.,1120.,960.2,878.0 &
        ,811.2,725.4,785.0,800.5,895.3,1049.,1222.,1435.,1513.,1541.,1374. &
        ,1138.,864.1,633.7,466.8,232.1,206.0,112.7,132.0,11.95,7.399, &
        4.947,2.857,2*0./

   ! grid levels for height integration   (p.s.h. distance = 0.25*IG)

   data IG/-16,-5,-1,0,1,7/

   integer I,J,II,JJ,IM
   real(REAL64) :: TT,ZO,H1,H2,H3,FJ,AA1,AA2,D1,D2

   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   do I=1,45
      FU(I)=exp(-960.24/T(I))
      FO3(I)=exp(-1500./T(I))
   enddo
   do I=46,69
      FU(I)=exp(-960.24/T(I))
   enddo
   do I=46,50
      FO3(I)=0.0
   enddo

   ! if ko=-1 then use Dickinson (1984) value for rate coefficient for
   !          the deactivation CO2(010) by O
   !    ko=0  then use Shved et al. (1991) value
   !    ko=1  then use Sharma&Wintersteiner (1990) value.

   do I=1,23
      TT=T(I+46)
      if (KO.eq.-1) then
         ZO = 0.5E7
      else if (KO.eq.0) then
         ZO = 3.5E7
      else
         ZO = (2.57E9*sqrt(TT)+1.70E13*exp(-76.75*TT**(-1./3.)))/TT
      endif
      AL(I)=1.5638/(1.5638+exp(-X(I+46))*(SN2(I)*(2.9*TT*TT-1060.*TT+ &
           145000.)+O2(I)*(4.23*TT*TT-1490.*TT+180000.)+O(I)*ZO))
   enddo

   ! *********************************************************************
   ! calculate the heating rates for layer below s.h.p. = 11.5
   !   15 um CO2 + 9.6 um O3:

   do I=1,9
      H1=FU(1)*AU(I,1)
      H2=FU(1)*BU(I,1)
      H3=FO3(1)*AO3(I,1)
      II=I+8

      do J=2,6
         JJ=II+IG(J)
         H1=H1+AU(I,J)*FU(JJ)
         H2=H2+BU(I,J)*FU(JJ)
         H3=H3+AO3(I,J)*FO3(JJ)
      enddo

      H(I)=H1+H2*FU(II)+H3*O3(II)
   enddo

   do I=10,35
      H1=0.
      H2=0.
      H3=0.
      II=I+8

      do J=1,6
         JJ=II+IG(J)
         H1=H1+AU(I,J)*FU(JJ)
         H2=H2+BU(I,J)*FU(JJ)
         H3=H3+AO3(I,J)*FO3(JJ)
      enddo

      H(I)=H1+H2*FU(II)+H3*O3(II)
   enddo

   !   above p.s.h. = 10.5 only 15 um CO2 band is considered

   do I=36,39
      H1=0.
      H2=0.
      II=I+8

      do J=1,6
         FJ=FU(II+IG(J))
         H1=H1+AU(I,J)*FJ
         H2=H2+BU(I,J)*FJ
      enddo

      H(I)=H1+H2*FU(II)
   enddo

   ! calculate the heating rate above s.p.h. = 11.5 (the 15 um CO2 band only)
   !  using the reccurence formula and boundary condition at p.s.h. = 11.5

   H1=H(39)/(CO2(1)*(1.-AL(1)))*1.1008E-10

   do I=2,23
      IM=I-1
      AA1=1.-AL(IM)*(1.-.25*XL(I)-.75*XL(IM))
      AA2=1.-AL(I)*(1.-.75*XL(I)-.25*XL(IM))
      D1=-.25*(XL(I)+3.*XL(IM))
      D2=.25*(3.*XL(I)+XL(IM))
      H2=(AA1*H1-D1*FU(IM+46)-D2*FU(I+46))/AA2
      H(I+38)=H2*CO2(I)*(1.-AL(I))*8.6301E9
      H1=H2
   enddo

   return
end subroutine PCOOL
