!
      integer i,j,k,imin,jmin,kmin,imax,jmax,kmax,err,no, &
              nijk(7,Ptopo_numproc),tnijk(7,Ptopo_numproc),nw,i0,in,j0,jn
      integer*8 dim_L
      real*8 sum,sumd2,moy,var,mind,maxd,fijk,npt_8, &
             minmax(3,max(3,Ptopo_numproc)),tminmax(3,max(3,Ptopo_numproc))
!
!--------------------------------------------------------------------
!
      nijk  (:,:) = 0
      minmax(:,:) = 0.0d0

      i0 = max(F_i0 - Ptopo_gindx(1,Ptopo_myproc+1) + 1, 1   -west *G_halox)
      in = min(F_in - Ptopo_gindx(1,Ptopo_myproc+1) + 1, l_ni+east *G_halox)
      j0 = max(F_j0 - Ptopo_gindx(3,Ptopo_myproc+1) + 1, 1   -south*G_haloy)
      jn = min(F_jn - Ptopo_gindx(3,Ptopo_myproc+1) + 1, l_nj+north*G_haloy)

      if ((i0.le.l_ni).and.(in.ge.1).and. &
          (j0.le.l_nj).and.(jn.ge.1) ) nijk(7,Ptopo_myproc+1) = 1

      if (nijk(7,Ptopo_myproc+1).gt.0) then
         sum   = 0.0
         sumd2 = 0.0
         imin  = i0
         jmin  = j0
         kmin  = F_k0
         imax  = i0
         jmax  = j0
         kmax  = F_k0
         maxd  = F_field(i0,j0,F_k0)
         mind  = F_field(i0,j0,F_k0)
!
         do k=F_k0,F_kn
         do j=j0,jn
         do i=i0,in
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

         minmax(1,Ptopo_myproc+1) = maxd
         minmax(2,Ptopo_myproc+1) = mind
         minmax(3,1) = sum
         minmax(3,2) = sumd2
         minmax(3,3) = (in-i0+1)*(jn-j0+1)*(F_kn-F_k0+1)
         dim_L       = (in-i0+1)*(jn-j0+1)*(F_kn-F_k0+1)

         nijk  (1,Ptopo_myproc+1)= imax+ Ptopo_gindx(1,Ptopo_myproc+1)- 1
         nijk  (2,Ptopo_myproc+1)= jmax+ Ptopo_gindx(3,Ptopo_myproc+1)- 1
         nijk  (3,Ptopo_myproc+1)= kmax
         nijk  (4,Ptopo_myproc+1)= imin+ Ptopo_gindx(1,Ptopo_myproc+1)- 1
         nijk  (5,Ptopo_myproc+1)= jmin+ Ptopo_gindx(3,Ptopo_myproc+1)- 1
         nijk  (6,Ptopo_myproc+1)= kmin

      endif

      nw = 3*max(3,Ptopo_numproc)
      call rpn_comm_REDUCE ( nijk, tnijk, 7*Ptopo_numproc, &
                       "MPI_INTEGER","MPI_SUM",0,"grid",err )
      call rpn_comm_REDUCE ( minmax, tminmax, nw, &
                       "MPI_DOUBLE_PRECISION","MPI_SUM",0,"grid",err )

      if (Ptopo_myproc.eq.0) then

         imax  = tnijk (1,1)
         jmax  = tnijk (2,1)
         kmax  = tnijk (3,1)
         imin  = tnijk (4,1)
         jmin  = tnijk (5,1)
         kmin  = tnijk (6,1)
         maxd  = tminmax(1,1)
         mind  = tminmax(2,1)

         do i=1,Ptopo_numproc
            if ( tnijk (7,i) .gt. 0 ) then
               fijk = tminmax(1,i)
               if (fijk .gt. maxd) then
                  maxd = fijk
                  imax = tnijk (1,i)
                  jmax = tnijk (2,i)
                  kmax = tnijk (3,i)
               else if (fijk .eq. maxd) then
                  if (kmax.gt.tnijk(3,i)) then
                      imax = tnijk (1,i)
                      jmax = tnijk (2,i)
                      kmax = tnijk (3,i)
                  else if (kmax.eq.tnijk(3,i).and.jmax.gt.tnijk(2,i)) then
                      imax = tnijk (1,i)
                      jmax = tnijk (2,i)
                      kmax = tnijk (3,i)
                  endif
               endif
               fijk = tminmax(2,i)
               if (fijk .lt. mind) then
                  mind = fijk
                  imin = tnijk (4,i)
                  jmin = tnijk (5,i)
                  kmin = tnijk (6,i)
               else if (fijk .eq. mind) then
                  if (kmin.gt.tnijk(6,i)) then
                      imin = tnijk (4,i)
                      jmin = tnijk (5,i)
                      kmin = tnijk (6,i)
                  else if (kmin.eq.tnijk(6,i).and.jmin.gt.tnijk(5,i)) then
                      imin = tnijk (4,i)
                      jmin = tnijk (5,i)
                      kmin = tnijk (6,i)
                  endif
               endif
            endif
         end do

         sum   = tminmax(3,1)
         sumd2 = tminmax(3,2)
         npt_8 = tminmax(3,3)
         moy   = sum / npt_8 
         var   = max(0.d0,(sumd2 + moy*moy*npt_8 - 2*moy*sum) / npt_8)
         var   = sqrt(var)

         no=F_no
         if (F_rx.lt.8) then
         write(Lun_out,98) no,F_nv_S,moy,var,imin,jmin,kmin,mind, &
                        imax,jmax,kmax,maxd,F_from_S
         else
         write(Lun_out,99) no,F_nv_S,moy,var,imin,jmin,kmin,mind, &
                        imax,jmax,kmax,maxd,F_from_S
         endif

      endif

 98   format (i4,a4,' Mean:',1pe14.6,' Std:',1pe14.6, &
              ' Min:[(',i4,',',i4,',',i3,')', &
              1pe14.6,']',' Max:[(',i4,',',i4,',',i3,')', &
              1pe14.6,']',a6)
 99   format (i4,a4,' Mean:',1pe22.12,' Std:',1pe22.12,/ &
              ' Min:[(',i4,',',i4,',',i3,')', &
              1pe22.12,']',' Max:[(',i4,',',i4,',',i3,')', &
              1pe22.12,']',a6)
