***s/r statfld4 - calcule la moyenne, la variance, le rminimum et 
*                le maximum d'un champs et imprime le resultat.
*
      subroutine statfld4(nomvar,typvar,ip1,ip2,ip3,date,etiket,
     $     f,ni,nj,nk)
      implicit none
      character*4 nomvar
      character*2 typvar
      integer ip1,ip2,ip3,date
      character*12 etiket
* 
      integer ni,nj,nk
      real f(ni,nj,nk)
*
*OBJECT
*     calcule et imprime: la moyenne    (moy)
*                         la variance   (var)
*                         le minimum et le maximum
*     du champ f   
* 
*     arguments:
*         - f       - champ sur lequel on veut faire des statistiques
*         - n       - dimensions du champ f
*         - champ   - identification du champ
*         - no      - compteur 
*         - from    - identification du module d'ou on fait l'appel 
*
*METHOD
*
*EXTERNALS
*
*AUTHOR   Michel Desgagne                   Nov   1992
*
*Revision
* 001     M. Lepine, Mars  2003 -  appel a convip pour afficher les niveaux
*
*HISTORY
*
**
      integer i,j,k
      real sum,moy,var,rmin,rmax
      integer imin,jmin,kmin,imax,jmax,kmax,kind,dat2,dat3
      CHARACTER*15 Level
      REAL      rlevel
c--------------------------------------------------------------------
c
c ** On calcule la moyenne.
c
      sum = 0.0
      do 1 k=1,nk
         do 1 j=1,nj
            do 1 i=1,ni
         sum = sum + f(i,j,k)
 1    continue
      moy = sum / float(ni*nj*nk)
c
c ** On calcule la variance
c
      sum = 0.0
      do 2 k=1,nk
         do 2 j=1,nj
            do 2 i=1,ni
               sum = sum + ((f(i,j,k) - moy)*(f(i,j,k) - moy))
 2    continue
      var = sqrt (sum / float(ni*nj*nk))
c
c ** On identifie le minimum et le maximum.
c
      imin = 1
      jmin = 1
      kmin = 1
      imax = 1
      jmax = 1
      kmax = 1
      rmax = f(1,1,1)
      rmin = f(1,1,1)
c
      do 3 k=1,nk
         do 3 j=1,nj
            do 3 i=1,ni
               if (f(i,j,k) .gt. rmax) then
                  rmax  = f(i,j,k)
                  imax = i
                  jmax = j
                  kmax = k
               endif
               if (f(i,j,k) .lt. rmin) then
                  rmin  = f(i,j,k)
                  imin = i
                  jmin = j
                  kmin = k
               endif
 3    continue
*
      CALL convip(ip1,rlevel,kind,-1,level,.true.)
*      call newdate(date,dat2,dat3,-3);
*      print *,'Debug date=',date,dat2,dat3/100
c       
c ** On imprime
c 
C      write(6,10) nomvar,typvar,level,ip1,ip2,ip3,date,etiket,
      write(6,10) nomvar,typvar,level,ip2,ip3,date,etiket,
     $     moy,var,imin,jmin+(kmin-1)*nj,rmin,
     $     imax,jmax+(kmax-1)*nj,rmax
C 10   format (' ',a4,1x,a2,1x,a15,' (',i9,') ',i4,1x,i3,1x,i9,1x,a12,1x,
 10   format (' ',a4,1x,a2,1x,a15,1x,i4,1x,i3,1x,i9,1x,a12,1x,
     $     ' Mean:',e12.6,' StDev:',e12.6,
     $     '  Min:[(',i3,',',i3,'):',
     $     e10.4,']',' Max:[(',i3,',',i3,'):',
     $     e10.4,']')
c
c----------------------------------------------------------------
      return
      end
