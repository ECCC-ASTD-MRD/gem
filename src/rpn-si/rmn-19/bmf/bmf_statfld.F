#ifdef SELF_TEST
      program idiot
      real Z(30,20,10)
      do k=1,10
      do j=1,20
      do i=1,30
        z(i,j,k)=SQRT((i-16.)**2+(j-11.)**2+(k-6.)**2)
      enddo
      enddo
      enddo
      call bmf_statfld(Z,'test',0,'IDIOT',1,30,1,20,10,1,1,1,30,20,10)
      stop
      end
#endif
***s/r bmf_statfld - calcule la moyenne, la variance, le minimum et 
*                 le maximum d un champs et imprime le resultat.
*
      subroutine bmf_statfld (F_field, F_nv_S, F_no, F_from_S,
     $                     minx,maxx,miny,maxy,lnk,
     $                     F_i0,F_j0,F_k0,F_in,F_jn,F_kn)
*
* 
      implicit none
      character*(*) F_nv_S , F_from_S
      integer minx,maxx,miny,maxy,lnk,
     $        F_i0,F_j0,F_k0,F_in,F_jn,F_kn,F_no,unf
      real F_field(minx:maxx,miny:maxy,lnk)
*
*author
*     M. Valin - derived from original version of statfld (M. Desgagne)
*
*object
*     calcule et imprime: la moyenne    (moy)
*                         la variance   (var)
*                         le minimum et le maximum du champ f   
* 
*arguments
*  Name        I/O                 Description
*----------------------------------------------------------------
* F_field       I         Field to be operated on
* F_nv_S        I         User provided string to define F_field
* F_no          I         Usually the timestep #
* F_from_S      I         Usually the name of the calling subroutine
* F_i0,F_j0     I         Global lower-left indexes of the sub-domain
*                            on which to perform statistics
* F_in,F_jn     I         Global upper-right indexes of the sub-domain
*                            on which to perform statistics
* F_k0,F_kn     I         Range of levels on which to perform statistics
*----------------------------------------------------------------
*
*implicits
*
      integer i,j,k,imin,jmin,kmin,imax,jmax,kmax
      real*8 sum,sumd2,moy,var,mind,maxd,fijk,npt_8
      integer Lun_out
*
*--------------------------------------------------------------------
*
      Lun_out=0
      npt_8 = 1.0d0*((F_in-F_i0+1)*(F_jn-F_j0+1)*(F_kn-F_k0+1))
*
      sum   = 0.0
      sumd2 = 0.0
      imin  = F_i0
      jmin  = F_j0
      kmin  = F_k0
      imax  = F_in
      jmax  = F_jn
      kmax  = F_kn
      maxd  = F_field(F_in,F_jn,F_kn)
      mind  = F_field(F_i0,F_j0,F_k0)
*
      do k=F_k0,F_kn
      do j=F_j0,F_jn
      do i=F_i0,F_in
         fijk = F_field(i,j,k)
         sum = sum + fijk
         sumd2 = sumd2 + fijk*fijk
         if (fijk .gt. maxd) then
            maxd = fijk
            imax = i
            jmax = j
            kmax = k
         endif
         if (fijk .lt. mind) then
            mind = fijk
            imin = i
            jmin = j
            kmin = k
         endif
      end do
      end do
      end do
*      
      moy = sum / npt_8 
      var = max(0.d0,1.0d0*(sumd2 + moy*moy*npt_8 - 2*moy*sum) / npt_8)
      var = sqrt(var)
*       
* ** On imprime
*
         write(Lun_out,99) F_no,F_nv_S,moy,var,imin,jmin,kmin,mind,
     $                                imax,jmax,kmax,maxd,F_from_S 
*
 98   format (i4,a4,' Mean:',e14.7,' Var:',e14.7,
     $        ' Min:[(',i3,',',i3,',',i3,')',
     $        e14.7,']',' Max:[(',i3,',',i3,',',i3,')',
     $        e14.7,']',a6)
 99   format (i4,a4,' Mean:',e22.14,' Var:',e22.14,/
     $        ' Min:[(',i3,',',i3,',',i3,')',
     $        e22.14,']',' Max:[(',i3,',',i3,',',i3,')',
     $        e22.14,']',a6)
*
*----------------------------------------------------------------
*
      return
      end 
