!
!   definition des fonctions thermodynamiques de base
!   pour les constantes, utiliser le common /consphy/
!     note: toutes les fonctions travaillent avec les unites s.i.
!           i.e. ttt en deg k, prs en pa, qqq en kg/kg
!          *** n. brunet - mai 90 ***
!          * revision 01 - mai 94 - n. brunet
!                          nouvelle version pour faibles pressions
!          * revision 02 - aout 2000 - j-p toviessi
!                          calcul en real(REAL64) ::
!          * revision 03 - sept 2000 - n. brunet
!                          ajout de nouvelles fonctions
!          * revision 04 - janv 2000 - j. mailhot
!                          fonctions en phase mixte
!          * revision 05 - dec 2001 - g. lemay
!                          double precision pour phase mixte
!          * revision 06 - avr 2002 - a. plante
!                          ajout des nouvelles fonctions fottvh et fotvht
!
!     fonction de tension de vapeur saturante (tetens) - ew ou ei selon tt
      foewf(ttt) = dmin1(dsign(ttns3w,                                     &
       Dble(ttt)-Dble(trpl)),dsign                                         &
       (ttns3i,Dble(ttt)-Dble(trpl)))*dabs(Dble(ttt)-Dble(trpl))/          &
       (Dble(ttt)-ttns4w+dmax1(0.d0,dsign                                  &
       (ttns4w-ttns4i,Dble(trpl)-Dble(ttt))))
      fomult(ddd)=ttns1*ddd
      foew(ttt) = fomult(dexp(foewf(ttt)))
!
!     fonction calculant la derivee selon t de  ln ew (ou ln ei)
      fodle(ttt)=(4097.93d0+dmax1(0.d0,dsign(1709.88d0,                    &
       Dble(trpl)-Dble(ttt))))                                             &
       /((Dble(ttt)-ttns4w+dmax1(0.d0,dsign(ttns4w-ttns4i,                 &
       Dble(trpl)-Dble(ttt))))*(Dble(ttt)-ttns4w+dmax1(0.d0                &
       ,dsign(ttns4w-ttns4i,Dble(trpl)-Dble(ttt)))))
!
!     fonction calculant l'humidite specifique saturante (qsat)
      foqst(ttt,prs) = Dble(eps1)/(dmax1(1.d0,Dble(prs)/foew(ttt))-        &
       Dble(eps2))
      foqstx(prs,ddd) = Dble(eps1)/(dmax1(1.d0,Dble(prs)/ ddd )- &
       Dble(eps2))
!
!     fonction calculant la derivee de qsat selon t
      fodqs(qst,ttt)=Dble(qst)*(1.d0+Dble(delta)*Dble(qst))*fodle(ttt)
!     qst est la sortie de foqst
!
!     fonction calculant tension vap (eee) fn de hum sp (qqq) et prs
      foefq(qqq,prs) = dmin1(dble(prs),(dble(qqq)*dble(prs)) /             &
       (dble(eps1) + dble(eps2)*dble(qqq)))
!
!      fonction calculant hum sp (qqq) de tens. vap (eee) et pres (prs)
      foqfe(eee,prs) = dmin1(1.d0,Dble(eps1)*Dble(eee)/(Dble(prs)-         &
       Dble(eps2)*Dble(eee)))
!
!      fonction calculant temp virt. (tvi) de temp (ttt) et hum sp (qqq)
      fotvt(ttt,qqq) = Dble(ttt) * (1.0d0 + Dble(delta)*Dble(qqq))

!      fonction calculant temp virt. (tvi) de temp (ttt), hum sp (qqq) et
!      masse sp des hydrometeores.
      fotvht(ttt,qqq,qqh) = Dble(ttt) *                                    &
           (1.0d0 + Dble(delta)*Dble(qqq) - Dble(qqh))
!
!      fonction calculant ttt de temp virt. (tvi) et hum sp (qqq)
      fottv(tvi,qqq) = Dble(tvi) / (1.0d0 + Dble(delta)*Dble(qqq))

!      fonction calculant ttt de temp virt. (tvi), hum sp (qqq) et
!      masse sp des hydrometeores (qqh)
      fottvh(tvi,qqq,qqh) = Dble(tvi) /                                    &
           (1.0d0 + Dble(delta)*Dble(qqq) - Dble(qqh))
!
!      fonction calculant hum rel de hum sp (qqq), temp (ttt) et pres (prs)
!      hr = e/esat
       fohr(qqq,ttt,prs) = dmin1(Dble(prs),foefq(qqq,prs)) / foew(ttt)
       fohrx(qqq,prs,ddd) = dmin1(Dble(prs),foefq(qqq,prs)) / ddd
!
!     fonction calculant la chaleur latente de condensation
      folv(ttt) =Dble(chlc) - 2317.d0*(Dble(ttt)-Dble(trpl))
!
!     fonction calculant la chaleur latente de sublimation
      fols(ttt) = dble(chlc)+dble(chlf)+(dble(cpv)-                        &
                  (7.24d0*dble(ttt)+128.4d0))*(dble(ttt)-dble(trpl))
!      fols(ttt) = Dble(chlc)+Dble(chlf) + (Dble(cpv) - Dble(cpi) -          &
!      3.895d0*(Dble(ttt)-3.d0*Dble(trpl)))*(Dble(ttt)-Dble(trpl))
!
!     fonction resolvant l'eqn. de poisson pour la temperature
!     note: si pf=1000*100, "fopoit" donne le theta standard
      fopoit(t00,pr0,pf)=Dble(t00)*(Dble(pr0)/Dble(pf))**                  &
                       (-Dble(cappa))
!
!     fonction resolvant l'eqn. de poisson pour la pression
      fopoip(t00,tf,pr0)=Dble(pr0)*dexp(-(dlog(Dble(t00)/Dble(tf))/        &
                       Dble(cappa)))
!
!     les 5 fonctions suivantes sont valides dans le contexte ou on
!     ne desire pas tenir compte de la phase glace dans les calculs
!     de saturation.
!   fonction de vapeur saturante (tetens)
      foewaf(ttt)=ttns3w*(Dble(ttt)-Dble(trpl))/(Dble(ttt)-ttns4w)
      foewa(ttt) = fomult(dexp(foewaf(ttt)))
!   fonction calculant la derivee selon t de ln ew
      fodla(ttt)=ttns3w*(Dble(trpl)-ttns4w)/(Dble(ttt)-ttns4w)**2
!   fonction calculant l'humidite specifique saturante
      foqsa(ttt,prs)=Dble(eps1)/(dmax1(1.d0,Dble(prs)/foewa(ttt))-         &
       Dble(eps2))
!   fonction calculant la derivee de qsat selon t
      fodqa(qst,ttt)=Dble(qst)*(1.d0+Dble(delta)*Dble(qst))*fodla(ttt)
!   fonction calculant l'humidite relative
      fohra(qqq,ttt,prs)=dmin1(Dble(prs),foefq(qqq,prs))/foewa(ttt)
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!   definition of basic thermodynamic functions in mixed-phase mode
!     fff is the fraction of ice and ddff its derivative w/r to t
!     note: s.i. units are used
!           i.e. ttt in deg k, prs in pa
!          *** j. mailhot - jan. 2000 ***
!
!     saturation calculations in presence of liquid phase only
!     function for saturation vapor pressure (tetens)
      fesif(ttt)=ttns3i*(Dble(ttt)-Dble(trpl))/(Dble(ttt)-ttns4i)
      fesi(ttt)=fomult(dexp(fesif(ttt)))
!
      fdlesi(ttt)=ttns3i*(Dble(trpl)-ttns4i)/(Dble(ttt)-ttns4i)**2

      fesmx(ttt,fff) = (1.d0-Dble(fff))*foewa(ttt)+Dble(fff)*fesi(ttt)
      fesmxx(fff,fesi8,foewa8)=(1.d0-Dble(fff))*foewa8+Dble(fff)*fesi8

      fdlesmx(ttt,fff,ddff) = ( (1.d0-Dble(fff))*foewa(ttt)*fodla(ttt)     &
                            + Dble(fff)*fesi(ttt)*fdlesi(ttt)              &
                  + Dble(ddff)*(fesi(ttt)-foewa(ttt)) )/fesmx(ttt,fff)
      fdlesmxx(ttt,fff,ddff,foewa8,fesi8,fesmx8) =                         &
                              ( (1.d0-Dble(fff))*foewa8*fodla(ttt)         &
                            + Dble(fff)*fesi8*fdlesi(ttt)                  &
                  + Dble(ddff)*(fesi8-foewa8) )/fesmx8

      fqsmx(ttt,prs,fff) = Dble(eps1)/                                     &
              (dmax1(1.d0,Dble(prs)/fesmx(ttt,fff) ) - Dble(eps2)  )
      fqsmxx(fesmx8,prs) = Dble(eps1)/                                     &
              (dmax1(1.d0,Dble(prs)/fesmx8 ) - Dble(eps2)  )

      fdqsmx(qsm,dlemx) = Dble(qsm ) *(1.d0 + Dble(delta)* Dble(qsm ) )    &
                           * Dble(dlemx )
