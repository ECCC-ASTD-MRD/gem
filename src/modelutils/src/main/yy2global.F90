!COMP_ARCH=inteloneapi-2022.1.2 ; -suppress=-C
!COMP_ARCH=PrgEnv-intel-5.2.82 ; -suppress=-C
!COMP_ARCH=intel-19.0.3.199 ; -suppress=-C
!COMP_ARCH=PrgEnv-intel-6.0.5 ; -suppress=-C
!COMP_ARCH=PrgEnv-intel-6.0.5-19.0.5.281 ; -suppress=-C

subroutine yy2global()
   use, intrinsic :: iso_fortran_env, only: REAL64, INT64
   ! for forecast (P) or analysis (A) fields
   use vGrid_Descriptors, only: vgrid_descriptor,vgd_new,vgd_put,vgd_write
   use clib_itf_mod
   implicit none

   !author  Vivian Lee  Jun.2011

   !
   ! pgm -iyy [yyfile] -oglob [globfile] -inml [nmlfile] [-print]
   !
   ! option print- to print out table for interpolation/clipping choices
   ! Default is cubic and NO clipping
   !
   ! namelist format: MUST be 'nomvar:interp[:clip]'
   !
   ! Example:
   ! &interp
   !   varlist_S=
   !   'ES:LINEAR:CLIP',
   !   'HR:CUBIC:CLIP',
   ! /
   !      yin2global for nearest intcode=0
   !      yin2global for linear  intcode=1
   !      yin2global for cubic   intcode=2
   !

   !revision
   ! 103- Correction for linear  interpolation in file yin2gauss
   ! 104- Correction for nearest interpolation in file yin2gauss
   ! 105- Change etiket to take the one in Yin
   ! 106- Set global grid with points exact from Yin grid and fstecr rewrit=false
   !      because npas is not seen as a different field! and change routine
   !      yin2gauss to yin2global
   ! 107- Correction to definition of the global grid (for lat points to
   !      start at 0, ends at 360, for lon points to start and end -89.9,89.9)
   ! 108- To read a namelist to determine interpolation scheme and clipping
   !      and yinuv2global accepts interpolation code
   ! 109- To only copy Yin data into the Yin part and interpolate Yang part
   !      into the global Yin grid 
   !    - NOTE: will change poles to -90.0,90.0 when EZ is ready
   ! 110- To correct routine localise, int_lin_lag
   ! 111- To take either CORE or MODEL Yin/Yang grids where the copy of the Yin
   !      data is taken only from the Yin (45 to 315 deg lon, -45 to 45 deg lat)
   !      And to be able to read U grids (combined Yin and Yang grid in one)
   ! 112- Use Clib to accept a directory or a file in -iyin,-iyan,-iyy
   !      Also corrected to write out VV fields with rewrit=false like the others
   ! 113- Correction to last end points of the Yaxis and Xaxis
   !      Resolution will be uniform throughout except ends of Xaxis
   ! 114- End points on the Yaxis must be less that Abs(90.0): 0.95 *deltaY
   ! 115- Define myp3 (to 0), delete check for UU,VV, check to
   !      NOT allow different grid sizes and force error prints to STDERR
   ! 116- Add new feature to copy Z grids and their grid parameters to output
   !      and remove option to accept separate Yin, Yang files.

   include "rmnlib_basics.inc"

   integer, parameter :: maxnfiles=80
   integer, parameter :: maxvar=300
   integer, parameter :: nklemax=6
   integer, parameter :: nrectot=20000
   real(REAL64), parameter :: pole_8= 90.0d0
   real(REAL64), parameter :: near_pole_8= 89.99999999
   real,    parameter :: pole  = 90.0

   integer npos
   type(vgrid_descriptor) :: vgd

   character(len=256) :: defaut(nklemax),liste(0:nklemax-1),val(nklemax)
   character(len=256) :: yyfile_S,globfile_S,nml_S
   character(len=256) :: filelist_yy(maxnfiles)

   integer i,j,k,m,ier,iunout,iunyy,NYY,iunnml
   integer, allocatable, dimension(:) :: unit_yy
   integer err_iunyy,yyfiles
   integer ip1,ip2,ip3,ni,nj,nk,datev,key,len,mni,mnj,xni,xnj
   integer dateo,ni1,nj1,nk1,p1,p2,p3,yni,ynj,gni,gnj,intcode

   real(REAL64) :: deltat_8,h1_8,h2_8,deg2rad_8,dxla_8,dyla_8,rad2deg_8
   integer ig1,ig2,ig3,ig4,scrap,myig1,myig2,myig3,myig4
   integer myni,mynj,myp1,myp2,myp3,pni,pnj
   integer deet,npas,nbits,datyp,nsubgrids,sgid

   character(len=12) :: etiket,yyetiket
   character(len=4)  :: nomvar
   character(len=1)  :: typvar,grtyp,mygrtyp
   character(len=17) :: varlist_S(MAXVAR)
   character(len=4)  :: myvar_S(MAXVAR),myclip_S,print_S
   character(len=7)  :: myinterp_S
   logical myclip(MAXVAR),clipcode,print_L,ugrid_L,yy2enml_L
   integer myinterp(MAXVAR),myvarnum,interp_index,clip_index
   namelist /interp/ varlist_S

   integer, allocatable,dimension(:) :: list_yy,subgrid
   real, allocatable, dimension(:,:) :: work,gg,vgg,gwork,zwork,zzwork
   real, allocatable,dimension(:,:) :: workyin,vworkyin,workyinyan,vworkyinyan
   real(REAL64), allocatable, dimension(:,:) :: workyan_8,gg_8
   real(REAL64), allocatable, dimension(:,:) :: vworkyan_8,vgg_8
   real, allocatable, dimension(:) :: xpx,ypx,xglob,yglob
   real(REAL64),allocatable,dimension(:) :: xpx_8,ypx_8
   real(REAL64), allocatable,dimension(:) :: G_xglob_8,G_yglob_8
   integer, allocatable,dimension (:) :: yan_imx,yan_imy,yan_i,yan_j
   real(REAL64),allocatable,dimension (:) :: yan_t,yan_p
   real(REAL64),allocatable,dimension (:,:) :: yan_s
   real(REAL64) :: pi,s(2,2)

   !  s(1,1)=s(1,:)
   !  s(1,2)=s(2,:)
   !  s(2,1)=s(3,:)
   !  s(2,2)=s(4,:)

   pi = acos(-1.d0)
   deg2rad_8 = (acos(-1.0d0))/180.
   rad2deg_8 = 180./(acos(-1.0d0))
   npos=0
   iunout=0
   iunyy=0
   iunnml=0
   yy2enml_L=.false.

   ier = exdb('yy2global','116','NON')

   liste(:)  = (/'iyy.       ','inml.      ','oglob.     ','print.     ','liste      ','help       '/)
   defaut(:) = (/'yinyang.fst','yy2e.nml   ','glob.fst   ','1          ','           ','AIDE       '/)
   val(:)    = (/'yinyang.fst','yy2e.nml   ','glob.fst   ','0          ','           ','OK         '/)
   call ccard(liste,defaut,val,nklemax,npos)

   ! do i=1,nklemax
   !    write(6, '(a)') 'val='//val(i)
   ! enddo
   yyfile_S=val(1)
   nml_S=val(2)
   globfile_S=val(3)
   print_S=val(4)
   print_L=.false.
   if (trim(print_S).eq.'1') print_L=.true.

   print *,'yyfile_S=',trim(yyfile_S)
   print *,'namelist file=',trim(nml_S)
   print *,'globfile_S=',trim(globfile_S)

   ! Regardez les repertoires ou les fichiers
   if (clib_isdir(yyfile_S).gt.0) then
      !     print *,'unit',iunyy,' for ',trim(yyfile_S)
      print *,'directory found for U grid files'
      write(*, '(a)') 'directory found for yyfile_S'
      ier=clib_glob(filelist_yy,yyfiles,trim(yyfile_S)//'/*',maxnfiles)
      write(*, '(a,1x,i0)') 'yyfiles found =',yyfiles
      if (yyfiles.eq.0) then
         write(0, '(a)') 'Problem: No files found in -iyy : ABORT'
         call qqexit(-1)
      endif
      allocate(unit_yy(yyfiles))
      unit_yy(:)=0
      do i=1,yyfiles
         err_iunyy = fnom(unit_yy(i), filelist_yy(i),'RND+OLD+R/O',0)
         if (err_iunyy.ne.0) then
            write(0, '(a)') 'Problem with fnom on file '//trim(filelist_yy(i))//': ABORT'
            call qqexit(-1)
         endif
         !         print *,i,'yinfile=',trim(filelist_yy(i))
         ier = FSTOUV(unit_yy(i), 'RND')
      enddo
      ier=fstlnk(unit_yy,yyfiles)
      iunyy=unit_yy(1)
   else if (clib_isfile(yyfile_S).gt.0) then
      err_iunyy = fnom(iunyy, yyfile_S,'RND+OLD+R/O',0)
      if (err_iunyy.ne.0) then
         write(0, '(a)') 'Problem with file -iyy : ABORT'
         call qqexit(-1)
      endif
      ier = FSTOUV(iunyy, 'RND')
   else
      write(0, '(a)') 'ERROR: Opening these input file(s)'
      write(0, '(a)') 'OR ERROR in: -iyy  '//trim(yyfile_S)
      call qqexit(-1)
   endif

   ! Ouverture du fichier namelist
   myvarnum=0
   ier= fnom(iunnml, nml_S,'SEQ+OLD',0)
   if(ier.ne.0) then
      write(*, '(a)') 'WARNING: file opening failed: for '//trim(nml_S)
   else
      yy2enml_L = .true.
      do i=1,MAXVAR
         varlist_S(i)=''
      enddo
      rewind(iunnml)
      !    print *,'unit',iunnml,' for ',trim(nml_S)
      read (iunnml,nml=interp,end=100,err=100)

      if (print_L) then
         print *,'MY INTERPOLATION TABLE'
         write(*,900)
         write(*,901)
         write(*,900)
      endif
      do i=1,MAXVAR
         if (varlist_S(i).ne.'') then
            myinterp_S='CUBIC'
            myclip_S=''
            interp_index=index(varlist_S(i),":")
            clip_index=index(varlist_S(i)(interp_index+1:),":")
            call low2up(varlist_S(i)(1:interp_index-1),myvar_S(i))
            !           print *,'interp_index,clip_index= ',interp_index,clip_index
            if (clip_index.gt.0) then
               clip_index=clip_index+interp_index
               call low2up(varlist_S(i)(interp_index+1:clip_index-1),myinterp_S)
               call low2up(varlist_S(i)(clip_index+1:),myclip_S)
            else
               call low2up(varlist_S(i)(interp_index+1:),myinterp_S)
            endif
            if (trim(myvar_S(i)).eq.'VV') myvar_S(i)='UU'

            if (trim(myinterp_S).eq.'NEAREST')myinterp(i)=0
            if (trim(myinterp_S).eq.'LINEAR')myinterp(i)=1
            if (trim(myinterp_S).eq.'CUBIC')myinterp(i)=2
            if (trim(myclip_S).eq.'CLIP') myclip(i)=.true.
            if (print_L) &
                 write(*,902)trim(myvar_S(i)),trim(myinterp_S),myinterp(i),myclip(i)
         else
            myvarnum=i-1
            exit
         endif
      enddo

   endif
   if (print_L) write(*, '(a,1x,i0,1x,a,1x,i0)') 'Number of variables found in namelist=',myvarnum,' MAXVAR=',maxvar
   if (print_L) write(*,900)

100 continue

   if (ier.eq.0) then
      if (myvarnum.eq.0) write(*, '(a)') 'WARNING: Namelist interp not found in '//trim(nml_S)
   endif
   if (myvarnum.eq.0) write(*, '(a)') 'WARNING:Will use CUBIC interpolation, NO CLIPPING for ALL variables'
900 format('+',30('-'),'+')
901 format('|',1x,'NOMVAR',1x,'|',' INTERP  ','|',1x,'code| clip|')
902 format('|',1x,a4,3x,'| ',a7,1x,'|',1x,i3,1x,'|',1x,L3,1x,'|')



   !   Open the output file
   ier= fnom(iunout, globfile_S,'RND+NEW',0)
   if(ier.ne.0) then
      write(0, '(a)') 'ERROR: file opening failed: abort '//trim(globfile_S)
      call qqexit(ier)
   else
      ier = FSTOUV(iunout, 'RND')
      !     print *,'unit',iunout,' for ',trim(globfile_S)
   endif
 
   !   Get the complete list of records to do
   allocate(list_yy(nrectot))
   ier = fstinl(iunyy,NI,NJ,NK,-1,' ',-1,-1,-1,'P',' ',list_yy,NYY,nrectot)
   if (NYY.eq.0) then
      print *,'no forecast fields found'
      ier = fstinl(iunyy,NI,NJ,NK,-1,' ',-1,-1,-1,'A',' ',list_yy,NYY,nrectot)
      if (NYY.eq.0) then
         print *,'no analysis fields found'
         ier = fstinl(iunyy,NI,NJ,NK,-1,' ',-1,-1,-1,'C',' ',list_yy,NYY,nrectot)
         if (NYY.eq.0) then
            write(0, '(a)') 'ERROR: No fields found!!'
            call qqexit(-1)
         else
            print *,NYY, 'climatology fields found'
         endif
      else
         print *,NYY, 'analysis fields found'
      endif
   else
      print *,NYY, 'forecast fields found'
   endif
   sgid=-1
   ugrid_L=.false.
   myni = -1
   mynj = -1

   !   Obtain the grid description from YIN or YINYANG
   do m=1,NYY
      ier = FSTPRM(list_yy(m), dateo, deet, npas, NI, NJ, NK,nbits,&
           datyp, IP1,&
           IP2, IP3, TYPVAR, NOMVAR, ETIKET, grtyp, IG1, IG2, IG3,&
           IG4, scrap, scrap, scrap, scrap, scrap, scrap, scrap)
      if (grtyp.eq.'U') then
         sgid=ezqkdef(NI,NJ,grtyp,IG1, IG2, IG3,IG4,iunyy)
         ugrid_L=.true.
         exit
      endif
   enddo

   IF (UGRID_L) THEN !---setup for Ugrid to Global grid

      myni = NI
      mynj = NJ
      p1 = ig1
      p2 = ig2
      p3 = 0

      !   Create output tags
      myp1=ig1+1
      myp2=ig2+1
      myp3=0

      !   Get latlons for YIN PHI grid
      nsubgrids=ezget_nsubgrids(sgid)
      allocate(subgrid(nsubgrids))
      ier=ezget_subgridids(sgid,subgrid)

      !vl TEST SECTION for YIN-YANG axis creation
      !    dx=1.0
      !    yni = 271
      !    ynj = 91
      !    allocate(xpx(yni),ypx(ynj),xpx_8(yni),ypx_8(ynj))
      !    xpx(1)=45.0
      !    do i=1,yni
      !       xpx(i+1)=xpx(i)+dx
      !    enddo
      !    ypx(1)=-45.0
      !    do i=1,ynj
      !       ypx(i+1)=ypx(i)+dx
      !    enddo

      !vl REAL SECTION for obtaining YIN-YANG axis
      !    Get rotation for Global grid from Yin
      ier=ezgxprm(subgrid(1),yni,ynj,grtyp,ig1,ig2,ig3,ig4,&
           mygrtyp,myig1,myig2,myig3,myig4)
      allocate(xpx(yni),ypx(ynj),xpx_8(yni),ypx_8(ynj))
      ier=gdgaxes(subgrid(1),xpx,ypx)

      !vl CREATION for Global Axis starts here
      do i=1,yni
         xpx_8(i)=dble(xpx(i))*deg2rad_8
         !      print *,'xpx_8',i,xpx_8(i),xpx(i)
      enddo
      do i=1,ynj
         ypx_8(i)=dble(ypx(i))*deg2rad_8
         !      print *,'ypx_8',i,ypx_8(i),ypx(i)
      enddo

      ! Resolution of GLOBAL Grid is from Yin
      dxla_8=xpx(2)-xpx(1)
      dyla_8=ypx(2)-ypx(1)
      h1_8=xpx_8(2)-xpx_8(1)
      h2_8=ypx_8(2)-ypx_8(1)

      print *,'*'
      print *,'yinyang input grids: YNI=',YNI,' YNJ=',YNJ,xpx(2)-xpx(1),'deg'
      print *,'xpx(2)-xpx(1)=',dxla_8,'deg', h1_8, 'radians'
      print *,'ypx(2)-ypx(1)=',dyla_8,'deg', h2_8, 'radians'
      print *,'myp1=',myp1,'myp2=',myp2,'myp3=',myp3
      print *,'*'

      ! Extract core Yin X-axis from 45 to 315 deg
      do i=1,yni
         if (xpx(i).gt.45.0) exit
      enddo
      xni = i-1

      ! Extract core Yin Y-axis from 45 to 315 deg
      do i=1,ynj
         if (ypx(i).gt.-45.0) exit
      enddo
      xnj = i-1

      ! Mni = size of Yin Xaxis, Mnj = size Yin Yaxis
      mni = yni-(2*xni)
      mnj = ynj-(2*xnj)
 
      ! Pni,Pnj = No. of end points X,Y to complete the global grid.
      pni = (xpx(xni)/dxla_8) + 1
      pnj = ((ypx(xnj)+pole)/dyla_8) + 1

      !   print *,'pni=(xpx(xni)/dxla_8 +1)=',pni,'xpx(xni=',xpx(xni),'(pni-1)*dxla_8=',(pni-1)*dxla_8
      !   print *,'pnj=(ypx(xnj)/dyla_8 +1)=',pnj,'ypx(xnj=',ypx(xnj),'(pnj-1)*dyla_8=',(pnj-1)*dyla_8

      !   (Make sure the endpoints do not go past a certain value)
      if (xpx(xni) - (pni-1)*dxla_8 .gt. 0.0d0) pni = pni+1
      if (ypx(xnj) - (pnj-1)*dyla_8 .lt.-90.0d0) pnj = pnj+1

      !   print *,'final pnj=',pnj
      !   print *,'xpx(xni)',xpx(xni),'xpx(yni-xni+1)',xpx(yni-xni+1)
      !   print *,'ypx(xnj)',ypx(xnj),'ypx(ynj-xnj+1)',ypx(ynj-xnj+1)
      !   print *,'xni=',xni,'yni-xni+1=',yni-xni+1
      !   print *,'xnj=',xnj,'ynj-xnj+1=',ynj-xnj+1
      !   print *,'mni=',mni,'mnj=',mnj
      !   print *,'xpx(xni)=',xpx(xni),'xpx(xni+mni+1)=',xpx(xni+mni+1)
      !   print *,'ypx(xnj)=',ypx(xnj),'ypx(xnj+mnj+1)=',ypx(xnj+mnj+1)

      ! Gni,Gnj = No of points for the global grid
      gni = mni+pni*2
      gnj = mnj+pnj*2
      print *,'Global Yin grid: gni=',gni,' gnj=',gnj

      allocate (xglob(gni),G_xglob_8(gni+2))
      allocate (yglob(gnj),G_yglob_8(gnj+2))
      G_xglob_8(pni)=xpx(xni)
      G_yglob_8(pnj)=ypx(xnj)
      xglob(pni)=xpx(xni)
      yglob(pnj)=ypx(xnj)
 
      print *,'*  '
      print *,'Global takes rotation from YIN: GNI=',GNI,' GNJ=',GNJ
      print *,'*  '
      allocate(gg_8(gni-1,gnj),vgg_8(gni-1,gnj))
      allocate(gg(gni,gnj),vgg(gni,gnj),gwork(gni,gnj))

      ! Calculate the global grid to put YIN point exactly on the same grid
      ! print *,'The global X array'
      !   Left end pt of X axis
      do i=pni,2,-1
         G_xglob_8(i-1) = G_xglob_8(i)-dxla_8
         xglob(i-1) =G_xglob_8(i-1)
      enddo

      !   Copy Yin Xaxis points
      do i=pni+1,mni+pni+1
         G_xglob_8(i)=xpx(xni-pni+i)
         xglob(i)=G_xglob_8(i)
      enddo
      !   print *,'xglob',pni+1,'=',xglob(pni+1)
      !   print *,'xglob',mni+pni+1,'=',xglob(mni+pni+1)

      !   Right end pt of X axis
      do i=mni+pni+2,gni-1
         G_xglob_8(i) = G_xglob_8(i-1)+dxla_8
         xglob(i)=G_xglob_8(i)
      enddo
      !   print *,'xglob',mni+pni+2,'=',xglob(mni+pni+2)
      G_xglob_8(1) = 0.0d0
      G_xglob_8(gni) = 360.0d0
      xglob(1) = G_xglob_8(1)
      xglob(gni)=G_xglob_8(gni)

      do i=1,gni
         G_xglob_8(i)=G_xglob_8(i)*deg2rad_8
         !    print *,'G_xglob_8',i,G_xglob_8(i),' xglob=',xglob(i)
      enddo

      !   G_yglob_8(1) = -pole_8
      !   G_yglob_8(gnj) = pole_8
      !   yglob(1)= -pole
      !   yglob(gnj)= pole

      !   print *,'The global Y array'
      !   Bottom end pt of Y axis
      yglob(pnj) =G_yglob_8(pnj)
      do j=pnj,2,-1
         G_yglob_8(j-1) = G_yglob_8(j)-dyla_8
         yglob(j-1) =G_yglob_8(j-1)
      enddo
      !   Copy Yin Yaxis points
      do j=pnj+1,mnj+pnj+1
         G_yglob_8(j)=ypx(xnj-pnj+j)
         yglob(j)=G_yglob_8(j)
      enddo
      !   print *,'yglob',pnj+1,'=',yglob(pnj+1)
      !   print *,'yglob',mnj+pnj+1,'=',yglob(mnj+pnj+1)

      !   Top end pt of Y axis
      do j=mnj+pnj+2,gnj
         G_yglob_8(j) = G_yglob_8(j-1)+dyla_8
         yglob(j) =G_yglob_8(j)
      enddo
      if (G_yglob_8(1)  .lt. -near_pole_8 .or. yglob(1)  .lt. -near_pole_8 .or. &
           G_yglob_8(gnj).gt.  near_pole_8 .or. yglob(gnj).gt.  near_pole_8 ) then
         G_yglob_8(1)   = G_yglob_8(2)    - 0.95*dyla_8
         G_yglob_8(gnj) = G_yglob_8(gnj-1)+ 0.95*dyla_8
         yglob(  1)= G_yglob_8(1  )
         yglob(gnj)= G_yglob_8(gnj)
      endif

      do j=1,gnj
         G_yglob_8(j)=G_yglob_8(j)*deg2rad_8
         !    print *,'G_yglob_8',j,G_yglob_8(j),'yglob=',yglob(j)
      enddo
      !   print *,'Endpoints in degrees'
      !   print *,'xglob(1)=',xglob(1)
      !   print *,'xglob(gni)=',xglob(gni)
      !   print *,'yglob(1)=',yglob(1)
      !   print *,'yglob(gnj)=',yglob(gnj)

      !   Allocate and precalculate values for interpolation to Yang grid

      len=0
      pi = acos(-1.d0)
      do j=1,gnj
         do i=1,gni-1
            if ( ( abs(G_xglob_8(i)).lt.xpx_8(xni).or.  &
                 abs(G_xglob_8(i)).gt.xpx_8(xni+mni+1) ).or. &
                 (    G_yglob_8(j) .lt.ypx_8(xnj)  .or.   &
                 G_yglob_8(j) .gt.ypx_8(xnj+mnj+1) ) )  then
               len=len+1
            endif
         enddo
      enddo
      allocate(yan_s(4,len),yan_t(len),yan_p(len))
      allocate(yan_imx(len),yan_imy(len),yan_i(len),yan_j(len))

      k=0
      do j=1,gnj
         do i=1,gni-1
            if ( ( abs(G_xglob_8(i)).lt.xpx_8(xni).or.  &
                 abs(G_xglob_8(i)).gt.xpx_8(xni+mni+1) ).or. &
                 (    G_yglob_8(j) .lt.ypx_8(xnj)  .or.   &
                 G_yglob_8(j) .gt.ypx_8(xnj+mnj+1) ) )  then
               k=k+1
               yan_i(k)=i
               yan_j(k)=j
               call smat(s,yan_t(k),yan_p(k),G_xglob_8(i)-pi,G_yglob_8(j))
               yan_t(k)=yan_t(k)+pi
               yan_s(1,k)=s(1,1)
               yan_s(2,k)=s(1,2)
               yan_s(3,k)=s(2,1)
               yan_s(4,k)=s(2,2)
               call localise (yan_imx(k),yan_imy(k),yan_t(k),yan_p(k),&
                    xpx_8,ypx_8,h1_8,h2_8,1,1)
            endif
         enddo
      enddo

      !   These work fields are for UU,VV calculations
      allocate(workyinyan(yni,ynj*2),vworkyinyan(yni,ynj*2))
      allocate (vworkyin(yni,ynj),vworkyan_8(yni,ynj))
      allocate (workyin(yni,ynj),workyan_8(yni,ynj))
      allocate (work(yni,ynj))

   ENDIF !---setup for Ugrid to Global grid

   print *,'Start re-assemble : Yin+Yang -> Global and add Z grids'
   print *,'*'

   !  Now interpolate fields
   DO_FIELDS: DO m=1,NYY
      ier = FSTPRM(list_yy(m), dateo, deet, npas, NI, NJ, NK,nbits,&
           datyp, IP1,IP2, IP3, TYPVAR, NOMVAR, ETIKET, &
           grtyp, IG1, IG2, IG3, IG4, &
           scrap, scrap, scrap, scrap, scrap, scrap, scrap)
      if (grtyp.eq.'Z') then
         allocate (zwork(ni,nj),zzwork(ni,nj))
         ier=fstluk(zwork,list_yy(m),NI,NJ,NK)
         ier = FSTECR(zwork, ZZWORK, -nbits, iunout, dateo,&
              deet, npas, ni, nj,1,&
              ip1, ip2, ip3, typvar,nomvar, etiket, 'Z',&
              ig1, ig2, ig3, ig4, datyp, .false.)
         deallocate(zwork,zzwork)
         cycle
      else if (myni.ne.ni .or. mynj.ne.nj) then
         write(0, '(a)') 'ERROR:' //trim(NOMVAR),' is not same grid size!!'
         call qqexit(-1)
      endif
      deltat_8=(dble(deet)*dble(npas))/3600.0
      call incdatr(datev,dateo,deltat_8)
      yyetiket=etiket(1:len_trim(etiket))
      if (typvar.eq.'C') datev=-1

      !  The size of these work fields depend if they are U,V, or PHI grids
      ier=fstluk(workyinyan,list_yy(m),NI,NJ,NK)
      workyin(:,:)=workyinyan(1:yni,1:ynj)
      work(:,:)=workyinyan(1:yni,ynj+1:ynj*2)
      do j=1,ynj
         do i=1,yni
            workyan_8(i,j)=dble(work(i,j))
         enddo
      enddo

      ! Special treatment for UU and VV variables
      if (nomvar.eq."UU") then
         intcode=2
         clipcode=.false.
         do i=1,myvarnum
            if (trim(nomvar).eq.trim(myvar_S(i))) then
               intcode= myinterp(i)
               clipcode=myclip(i)
            endif
         enddo
         key = fstinf(iunyy,NI1,NJ1,NK1,datev,' ',ip1,ip2,ip3,typvar,"VV")
         if(key.lt.0) then
            write(0, '(a)') 'ERROR: cannot find VV in yin file: abort'
            call qqexit(key)
         endif

         ier=fstluk(vworkyinyan,key,NI,NJ,NK)
         vworkyin(:,:)=vworkyinyan(1:yni,1:ynj)
         work(:,:)=vworkyinyan(1:yni,ynj+1:ynj*2)
         do j=1,ynj
            do i=1,yni
               vworkyan_8(i,j)=dble(work(i,j))
            enddo
         enddo
         call yanuv2global(gg_8,vgg_8,gni-1,gnj, workyan_8,vworkyan_8,&
              yni,ynj,xpx_8,ypx_8,yan_i,yan_j,yan_imx,yan_imy,yan_t,yan_p,yan_s,&
              len,h1_8,h2_8,intcode,clipcode)
         do j=1,gnj
            do i=1,gni-1
               gg(i,j)=real(gg_8(i,j))
               vgg(i,j)=real(vgg_8(i,j))
            enddo
            gg(gni,j) = gg(1,j)
            vgg(gni,j) = vgg(1,j)
         enddo

         ! Copy Yin data into global Yin array
         do j=1,mnj+2
            do i=1,mni+2
               gg(pni+i-1,pnj+j-1)=workyin(xni+i-1,xnj+j-1)
               vgg(pni+i-1,pnj+j-1)=vworkyin(xni+i-1,xnj+j-1)
            enddo
         enddo


         ier = FSTECR(gg, gWORK, -nbits, iunout, dateo,&
              deet, npas, gni, gnj,1,&
              ip1, ip2, ip3, typvar,'UU', yyetiket, 'Z',&
              myp1, myp2, myp3, 0, datyp, .false.)

         ier = FSTECR(vgg, gWORK, -nbits, iunout, dateo,&
              deet, npas, gni, gnj,1,&
              ip1, ip2, ip3, typvar,'VV', yyetiket, 'Z',&
              myp1, myp2, myp3, 0, datyp, .false.)

      else if (nomvar.ne."VV".and.nomvar.ne."UT1".and.nomvar.ne."VT1") then
         intcode=2
         clipcode=.false.
         do i=1,myvarnum
            if (trim(nomvar).eq.trim(myvar_S(i))) then
               intcode= myinterp(i)
               clipcode=myclip(i)
            endif
         enddo

         call yan2global(gg_8,gni-1,gnj,workyan_8,yni,ynj,&
              xpx_8,ypx_8,yan_i,yan_j,yan_imx,yan_imy,yan_t,yan_p,&
              len,h1_8,h2_8,intcode,clipcode)
         do j=1,gnj
            do i=1,gni-1
               gg(i,j)=real(gg_8(i,j))
            enddo
            gg(gni,j) = gg(1,j)
         enddo
         ! Copy Yin data into global Yin array
         do j=1,mnj+2
            do i=1,mni+2
               gg(pni+i-1,pnj+j-1)=workyin(xni+i-1,xnj+j-1)
            enddo
         enddo

         ier = FSTECR(gg, gWORK, -nbits, iunout, dateo,&
              deet, npas, gni, gnj,1,&
              ip1, ip2, ip3, typvar,nomvar, yyetiket, 'Z',&
              myp1, myp2, myp3, 0, datyp, .false.)
         ! IEEE output
         !          ier = FSTECR(gg, gWORK, -32, iunout, dateo,&
         !                       deet, npas, gni, gnj,1,&
         !                       ip1, ip2, ip3, typvar,nomvar, yyetiket, 'Z',&
         !                       myp1, myp2, myp3, 0, 5, .false.)
      endif
   enddo DO_FIELDS
   if (ugrid_L) deallocate (workyin,vworkyin,workyan_8,vworkyan_8,work)

   ! Write out the tictic,tactac and toctoc for new global Z grid
   ! Only one U grid is assumed in input file which is interpolated to new Z grid
   if (ugrid_L) then
      ier = FSTECR(xglob, gWORK, -32, iunout, dateo,&
           0, 0, gni, 1,1,&
           myp1, myp2, myp3, 'X','>>', yyetiket, mygrtyp,&
           myig1, myig2, myig3, myig4, 5, .true.)
      ier = FSTECR(yglob, gWORK, -32, iunout, dateo,&
           0, 0, 1, gnj,1,&
           myp1, myp2, myp3, 'X','^^', yyetiket, mygrtyp,&
           myig1, myig2, myig3, myig4, 5, .true.)
      key = fstinf(iunyy,NI1,NJ1,NK1,-1,' ',p1,p2,p3,' ','!!')
      if (key.ge.0) then
         ier = vgd_new(vgd,unit=iunyy,format='fst',ip1=p1,ip2=p2)
         ier = vgd_put(vgd,key='IP_1 - record ip1',value=myp1)
         ier = vgd_put(vgd,key='IP_2 - record ip2',value=myp2)
         ier = vgd_write(vgd,unit=iunout,format='fst')
      endif
      deallocate(xglob,yglob)
   endif

   ! Write out the tictac pairs for the other Z grid(s)
   ! Assume multiple Z grids from input file to be copied (no interp) over

   ! Re-using variable list_yy, NYY to obtain number of tictac pair records
   ier = fstinl(iunyy,NI,NJ,NK,-1,' ',-1,-1,-1,'X','>>',list_yy,NYY,nrectot)
   if (NYY.gt.0) then
      do m=1,NYY
         ier = FSTPRM(list_yy(m), dateo, deet, npas, NI, NJ, NK,nbits,&
              datyp, IP1,IP2, IP3, TYPVAR, NOMVAR, ETIKET, &
              grtyp, IG1, IG2, IG3, IG4, &
              scrap, scrap, scrap, scrap, scrap, scrap, scrap)
         allocate(xglob(ni),zwork(ni,nj))
         ier=fstluk(xglob,list_yy(m),NI,NJ,NK)
         ier = FSTECR(xglob, ZWORK, -nbits, iunout, dateo,&
              deet, npas, ni, nj,1,&
              ip1, ip2, ip3, typvar,nomvar, etiket, grtyp,&
              ig1, ig2, ig3, ig4, datyp, .false.)
         deallocate(xglob,zwork)
         key = fstinf(iunyy,NI1,NJ1,NK1,-1,' ',ip1,ip2,ip3,typvar,"^^")
         ier = FSTPRM(key, dateo, deet, npas, NI, NJ, NK,nbits,&
              datyp, IP1,IP2, IP3, TYPVAR, NOMVAR, ETIKET, &
              grtyp, IG1, IG2, IG3, IG4, &
              scrap, scrap, scrap, scrap, scrap, scrap, scrap)
         allocate(yglob(nj),zwork(ni,nj))
         ier=fstluk(yglob,key,NI,NJ,NK)
         ier = FSTECR(yglob, ZWORK, -nbits, iunout, dateo,&
              deet, npas, ni, nj,1,&
              ip1, ip2, ip3, typvar,nomvar, etiket, grtyp,&
              ig1, ig2, ig3, ig4, datyp, .false.)
         deallocate(zwork,yglob)
         key = fstinf(iunyy,NI1,NJ1,NK1,-1,' ',ip1,ip2,0,' ','!!')
         if (key.ge.0) then
            ier = vgd_new(vgd,unit=iunyy,format='fst',ip1=ip1,ip2=ip2)
            ier = vgd_put(vgd,key='IP_1 - record ip1',value=ip1)
            ier = vgd_put(vgd,key='IP_2 - record ip2',value=ip2)
            ier = vgd_write(vgd,unit=iunout,format='fst')
         endif
      enddo
   endif

   ! Deallocate the rest of the variables and close the files
   if (ugrid_L) deallocate(xpx_8,ypx_8,xpx,ypx,gg_8,vgg_8,gg,vgg)


   ier=fstfrm(iunyy)
   ier=fstfrm(iunout)
   ier = fclos(iunyy)
   ier = fclos(iunout)
   if (yy2enml_L) ier = fclos(iunnml)
   ier = exfin('yy2global','116','NON')
   stop

end subroutine yy2global


subroutine yan2global(gg,gni,gnj,G2,ni,nj,x,y,ix,jx,&
     imx,imy,t,p,len,udx,udy,interp,clip)
   use, intrinsic :: iso_fortran_env, only: REAL64
   !author Abdessamad Qaddouri Nov.2009
   !
   !revision 103- V.Lee Jun.2010 - corr for linear int.
   !revision 104- V.Lee Sep.2011 - corr for nearest int.
   !revision 106- V.Lee Jan.2012 - take full Yin grid onto global
   !revision 107- V.Lee May.2012 - allow choice of interpolation for U,V
   !revision 109- V.Lee May.2012 - only interpolate Yang values
   !revision 110- V.Lee Aug.2012 - correct int_lin_lag routine
   !revision 111- V.Lee Feb.2013 - rename subroutines with prefix yy
   !
   ! gg - sortie sur grille defini gx,gy
   ! gni,gnj - dim  pour gg
   ! g2   - les champs de yan
   ! ni,nj   - dim pour yin,yan
   ! x,y     - coord de yinyan in radians
   ! gx,gy   - coord de gg in radians
   ! udx,udy  - h1,h2 dans yin,yan
   ! interp   - 0-nearest, 1-linear, 2-cubic
   !
   implicit none
   include "rmnlib_basics.inc"

   integer len,k
   integer gni,gnj,ni,nj,imx(len),imy(len),ix(len),jx(len),i,j,interp
   logical clip
   real(REAL64) ::  gg(gni,gnj)
   real(REAL64) ::  g2(ni,nj)
   real(REAL64) ::  x(ni),y(nj)

   real(REAL64) :: t(len),p(len)
   real(REAL64) ::  udx,udy

   ! interpolate geop to output grid
   !     print *,'x(1)=',x(1),' y(1)=',y(1)
   !     print *,'x(ni)=',x(ni),' y(nj)=',y(nj)
   !     print *,'udx=',udx,'udy=',udy

   !in Yang grid
   do k=1,len
      i=ix(k)
      j=jx(k)
      if (interp.eq.0) then
         call yy_int_near_lag(gg(i,j),g2,imx(k),imy(k),1,1, &
              ni,nj,t(k),p(k),x,y,udx,udy)
      elseif (interp.eq.1) then
         call yy_int_lin_lag(gg(i,j),g2,imx(k),imy(k),1,1,  &
              ni,nj,t(k),p(k),x,y)
      elseif (interp.eq.2) then
         call yy_int_cub_lag(gg(i,j),g2,imx(k),imy(k),1,1,  &
              ni,nj,t(k),p(k),x,y)
      else
         write(0, '(a)') 'Error interpolation code!!'
         call qqexit(-1)
      endif
      if (clip) then
         gg(i,j)=dmax1(gg(i,j),0.0d0)
      endif
   enddo

   return
end subroutine yan2global


subroutine yanuv2global(gu,gv,gni,gnj,U2,V2,ni,nj, x,y, &
     ix,jx,imx,imy,t,p,s,len,udx,udy,interp,clip)
   use, intrinsic :: iso_fortran_env, only: REAL64
   !revision V.Lee Nov.2009
   !
   ! gu,gv - sortie sur grille gx,gy
   ! gni,gnj - dim  pour gg
   ! u2,v2   - les champs de yan
   ! ni,nj   - dim pour yin,yan
   ! x,y     - coord de yinyan in radians
   ! gx,gy   - coord de gg in radians
   ! udx,udy  - h1,h2 dans yin,yan
   !
   implicit none
   include "rmnlib_basics.inc"

   integer len,k
   integer gni,gnj,ni,nj,imx(len),imy(len),i,j,interp,ix(len),jx(len)
   logical clip
   real(REAL64) ::  gu(gni,gnj),gv(gni,gnj)
   real(REAL64) ::  u2(ni,nj),v2(ni,nj)
   real(REAL64) ::  gu2(gni,gnj),gv2(gni,gnj)
   real(REAL64) ::  x(ni),y(nj)
   real(REAL64) ::  s(4,len),t(len),p(len)
   real(REAL64) ::  udx,udy

   ! interpol geop to defined output grid

   do k=1,len
      i=ix(k)
      j=jx(k)

      if (interp.eq.0) then
         call yy_int_near_lag(gu2(i,j),u2,imx(k),imy(k),1,1, &
              ni,nj,t(k),p(k),x,y,udx,udy)
         call yy_int_near_lag(gv2(i,j),v2,imx(k),imy(k),1,1, &
              ni,nj,t(k),p(k),x,y,udx,udy)
      elseif (interp.eq.1) then
         call yy_int_lin_lag(gu2(i,j),u2,imx(k),imy(k),1,1,  &
              ni,nj,t(k),p(k),x,y)
         call yy_int_lin_lag(gv2(i,j),v2,imx(k),imy(k),1,1,  &
              ni,nj,t(k),p(k),x,y)
      elseif (interp.eq.2) then
         call yy_int_cub_lag(gu2(i,j),u2,imx(k),imy(k),1,1, &
              ni,nj,t(k),p(k),x,y)
         call yy_int_cub_lag(gv2(i,j),v2,imx(k),imy(k),1,1, &
              ni,nj,t(k),p(k),x,y)
      endif
      gu(i,j)=s(1,k)*gu2(i,j) + s(2,k)*gv2(i,j)
      gv(i,j)=s(4,k)*gv2(i,j) + s(3,k)*gu2(i,j)
      if (clip) then
         gu(i,j)=dmax1(gu(i,j),0.0d0)
         gv(i,j)=dmax1(gv(i,j),0.0d0)
      endif
   enddo

   return
end subroutine yanuv2global


subroutine yy_int_cub_lag(FF,F,Imx,Imy,Ni,Nj,Nx,Ny,Xi,Yi,x,y)
   use, intrinsic :: iso_fortran_env, only: REAL64
   !author Abdessamad Qaddouri Nov.2009
   implicit none
   include "rmnlib_basics.inc"

   integer Ni,Nj,Imx(Ni,Nj),Imy(Ni,Nj),Nx,Ny
   integer i,j,Mx(Ni,Nj),My(Ni,Nj)
   real(REAL64) ::  W1,W2,W3,W4,X1,XX,X2,X3,X4
   integer Im, Jm
   real(REAL64) :: YY,y1,y2,y3,y4,FF(Ni,Nj)
   real(REAL64) :: F(Nx,Ny),fx1,fx2,fx3,fx4,x(*),y(*)
   real(REAL64) :: Xi(Ni,Nj),Yi(Ni,Nj)

   !
   do j =1,Nj
      do i= 1,Ni
         Mx (i,j)= min( max(Imx(i,j)-1,1), Nx - 3 )
         My (i,j)= min( max(Imy(i,j)-1,1), Ny - 3 )
      enddo
   enddo

   Do j=1, Nj
      Do i = 1, Ni
         Im = Mx(i,j)
         Jm = My(i,j)
         X1  = X(im)
         X2  = X(im+1) - X1
         X3  = X(im+2) - X1
         X4  = X(im+3) - X1
         XX  = Xi(i,j)   - X1
         !
         W1  = (XX-X2)*(XX-X3)*(XX-X4)/(-X2*X3*X4)
         W2  =  XX    *(XX-X3)*(XX-X4)/(X2*(X2-X3)*(X2-X4))
         W3  =  XX    *(XX-X2)*(XX-X4)/(X3*(X3-X2)*(X3-X4))
         W4  =  XX    *(XX-X2)*(XX-X3)/(X4*(X4-X2)*(X4-X3))

         Fx1 = W1*F(Im ,Jm)+W2*F(Im + 1,Jm )+ W3* F( Im + 2, Jm )+ &
              W4 * F( Im + 3, jm )
         Fx2 =W1*F(Im ,Jm +1)+W2*F(Im +1,jm +1)+W3*F(Im +2,Jm +1 )+ &
              W4 * F( Im + 3, Jm + 1 )
         Fx3 =W1*F(Im,Jm+2)+W2*F(Im+1,jm+2)+W3*F(Im+2,Jm+2)+ &
              W4 * F( Im + 3, Jm + 2 )
         Fx4 =W1*F(Im,jm+3)+W2*F(Im+1,jm+3)+W3*F(Im +2,Jm+3)+ &
              W4 * F( Im + 3, Jm + 3 )

         Y1  = Y(Jm)
         Y2  = Y(Jm+1) - Y1
         Y3  = Y(Jm+2) - Y1
         Y4  = Y(Jm+3) - Y1
         YY  = Yi(i,j)   - Y1
         !
         W1  = (yy-y2)*(yy-y3)*(yy-y4)/(-y2*y3*y4)
         W2  =  yy    *(yy-y3)*(yy-y4)/(y2*(y2-y3)*(y2-y4))
         W3  =  yy    *(yy-y2)*(yy-y4)/(y3*(y3-y2)*(y3-y4))
         W4  =  yy    *(yy-y2)*(yy-y3)/(y4*(y4-y2)*(y4-y3))

         FF(i,j)= W1*fx1+W2*fx2+W3*fx3+W4*fx4

      enddo
   enddo
   return
end subroutine yy_int_cub_lag


subroutine yy_int_lin_lag(FF,F,Imx,Imy,Ni,Nj,Nx,Ny,Xi,Yi,x,y)
   use, intrinsic :: iso_fortran_env, only: REAL64
   !author Abdessamad Qaddouri Nov.2009
   implicit none
   include "rmnlib_basics.inc"

   integer Ni,Nj,Imx(Ni,Nj),Imy(Ni,Nj),Nx,Ny
   integer i,j,Mx(Ni,Nj),My(Ni,Nj)

   integer Im, Jm
   real(REAL64) :: FF(Ni,Nj)
   real(REAL64) :: F(Nx,Ny),x(*),y(*)
   real(REAL64) :: Xi(Ni,Nj),Yi(Ni,Nj)
   real(REAL64) :: betax,betax1,betay,betay1

   do j =1,Nj
      do i= 1,Ni
         Mx (i,j)= min( max(Imx(i,j),1), Nx - 1 )
         My (i,j)= min( max(Imy(i,j),1), Ny - 1 )
      enddo
   enddo

   Do j=1, Nj
      Do i = 1, Ni
         Im = Mx(i,j)
         Jm = My(i,j)
         ! Min and max over fraction values to be non-negative and not gt 1.0
         betax= (Xi(i,j)-x(Im))/(X(im+1) - X(im))
         betax1= (1.0d0-betax)
         betay= (Yi(i,j)-y(Jm))/(Y(Jm+1) - Y(jm))
         betay1= (1.0d0-betay)

         FF(i,j)= betay1*(betax1*F(Im,Jm)+betax*F(Im+1,Jm))+ &
              betay*(betax1*F(Im,Jm+1)+betax*F(Im+1,Jm+1))
      Enddo
   Enddo
   return
end subroutine yy_int_lin_lag


subroutine yy_int_near_lag(FF,F,Imx,Imy,Ni,Nj,Nx,Ny,Xi,Yi,x,y,h1,h2)
   use, intrinsic :: iso_fortran_env, only: REAL64
   !author Abdessamad Qaddouri Nov.2009
   implicit none
   include "rmnlib_basics.inc"

   integer Ni,Nj,Imx(Ni,Nj),Imy(Ni,Nj),Nx,Ny
   integer i,j,Mx(Ni,Nj),My(Ni,Nj)

   integer Im, Jm
   real(REAL64) :: FF(Ni,Nj)
   real(REAL64) :: F(Nx,Ny),x(*),y(*)
   real(REAL64) :: Xi(Ni,Nj),Yi(Ni,Nj),h1,h2
   real(REAL64) :: betax,betax1,betay,betay1
   !
   do j =1,Nj
      do i= 1,Ni
         Mx (i,j)= min( max(Imx(i,j),1), Nx-1 )
         My (i,j)= min( max(Imy(i,j),1), Ny-1 )
      enddo
   enddo
   do j=1, Nj
      do i = 1, Ni
         Im = Mx(i,j)
         Jm = My(i,j)
         betax= ( Xi(i,j)-x(Im))/h1
         betax1= (1.0-betax)
         if (betax.le.betax1) then
            betax=0.0d0
            betax1=1.0d0
         else
            betax=1.0d0
            betax1=0.0d0
         endif
         betay=(Yi(i,j)-y(Jm))/h2
         betay1=1.0-betay
         if (betay.le.betay1) then
            betay=0.0d0
            betay1=1.0d0
         else
            betay=1.0d0
            betay1=0.0d0
         endif
         FF(i,j)= betay1*(betax1*F(Im,Jm)+betax*F(Im+1,Jm))+ &
              betay*(betax1*F(Im,Jm+1)+betax*F(Im+1,Jm+1))

      enddo
   enddo
   return
end subroutine yy_int_near_lag

