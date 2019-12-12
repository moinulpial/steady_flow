

	complex alfa(15500,550),alfadx(15500,550)
      REAL*8 GRAV,DOMEGA,SPEED,DX,DY,DZ,GRENR,GRENI,DUM

     *      ,SGR,SGI,SXGR,SXGI,angle,SYR,SYI,SZGR,SZGI,szr,szi

      COMMON/PARAMT/ GRAV,DOMEGA,SPEED,DX,DY,DZ,DUM

      COMMON/RESULT/ GRENR(4),GRENI(4)


      COMMON/WAVE  / NTW,TWAVE(5),IWW,FN,FORWRD

      DIMENSION ENXYZ(3),X(15000),Y(15000),Z(15000),
     &         XG(500),YG(500),ZG(500)
C---Parameter set


	open(4,file='output/centroid',status='unknown')
	open(5,file='data/fine_surf2',status='unknown')


	open(7,file='output/grn_surf',status='unknown')
	open(44,file='output/dxgrn1',status='unknown')


c
	Write(*,*)'how many points in surface'
	Write(*,*)'FOR PRESENT CASE, PLEASE WRITE  2772'
	read(*,*)nsurp
	Write(*,*)'how many elements in body'
	Write(*,*)'FOR PRESENT CASE, PLEASE WRITE  250'
	read(*,*)nelem




cc	fn= froude number
cc	alg= length of the ship
CC	DOMEGA= FREQUENCY OF OSCILLATION


	PI  = 3.141592653589793238D0
c	alpha=26.565
	grav=9.8d0
	domega=.001
c	domega=2.155
	fn=0.2670
c	fn=0.50
	alg=2.0

	speed=fn*sqrt(grav*alg)


	do 1 i=1,nsurp
	do 1 j=1,nelem
	alfa(i,j)=(0.0,0.0)
	alfadx(i,j)=(0.0,0.0)
1	continue


	do 11 i = 1,nsurp
c	dr=dr-0.1
	read(5,*)x(i),y(i),z(i)
c	z(i)=0.0
	Write(*,*)'starting calculations for points' ,i
C	Write(*,*)x(I),y(I),z(I)
	do 111 j=1,nelem

	read(4,*)xg(J),yg(J),zg(J)

	XGZAI=X(I)-XG(J)
	YETA=Y(I)-YG(J)
	ZZ=Z(I)+ZG(J)
	ZZ1=Z(I)-ZG(J)	

c	WRITE(1,*)XGZAI,YETA,ZZ,ZZ1

	CALL GRTRM2(XGZAI,YETA,ZZ,ZZ1,grr,gri,sxr,sxi)



c	x=dr*cos(alpha*pi/180.0d0)
c	y=dr*sin(alpha*pi/180.0d0)
c	z=0.0
c	call grtrm2(x,y,z,grr,gri)


	alfa(i,j)=cmplx(grr,gri)
	alfadx(i,j)=cmplx(sxr,sxi)

	
C	Write(6,10)x,y,z,xg,yg,zg,GRR,GRI
	write(7,*) alfa(i,j)
	write(44,*) alfadx(i,j)


   10 format(6(2x,f8.4),2f18.11)
		
  111 continue
	Write(*,*)'finished calculations for points' ,i
	Write(*,*)' ======================================== '
	Write(*,*)
	rewind 4
c	dr=dr-0.02
   11 continue
	stop
	end 



cz	SUBROUTINE GRTRM2(dxx,dyy,dzz,x0,y0,z0,grr,gri)
cz	SUBROUTINE GRTRM2(dxx,dyy,dzz,x0,y0,z0,grr,gri,sxr,sxi)

	SUBROUTINE GRTRM2(XGZAI,YETA,ZZ,ZZ1,grr,gri,sxr,sxi)

      REAL*8 GRAV,DOMEGA,SPEED,DX,DY,DZ,GRENR,GRENI,DUM

     *      ,SGR,SGI,SXGR,SXGI,angle,SYR,SYI,SZGR,SZGI,szr,szi

c      COMMON/FILE  / KEY(72),ITI,ITO
      COMMON/PARAMT/ GRAV,DOMEGA,SPEED,DX,DY,DZ,DUM

      COMMON/RESULT/ GRENR(4),GRENI(4),agreen(6000,100),agr(6000,100)

c      COMMON/IWAVE / DUMMY(3),WD(2)

      COMMON/WAVE  / NTW,TWAVE(5),IWW,FN,FORWRD,x,y,z

      DIMENSION  ENXYZ(3)

C---Parameter set

c	open(3,file='anext.out',status='unknown')	
c	open(3,file='agr.out',status='unknown')	
c	open(4,file='agrx.out',status='unknown')
c	open(5,file='agry.out',status='unknown')
c	open(6,file='agrz.out',status='unknown')

	PI  = 3.141592653589793238D0
	

      P4=-0.07957747d0

C----------------------------------------------------------------

C     DX    - (x - guzai)

C     DY    - (y - eta  )

C     DZ    - (z + zeta )

C     DOMEGA - circuler frequency of wave motion

C     ENXYZ - direction cosine of Normal with respect to (x,y,z) coord.

C     SS    - area of element

C------------------- output parameter ---------------------------

C     GRR   - real part of potential

C     GRI   - imaginary part of potential

C     GRDNR - real part of Normal-derivative

C     GRDNI - imagnary part of Normal-derivative

C----------------------------------------------------------------

c-----the next potion is the main prog... by zaman 12/12/02


	if(domega .eq.0) then 
	go to 1

	end if


      IF (SPEED.LE.0.005D0) SPEED=0.005D0


C---Cal. Potential & Derivative of unit Source
c	write(5,*)'point , real and imaginary part of green function,// ' 
    
c	do 2003 jj=0,180,5	
c	jj=31.0
	
ct 21/1 for gnu plot
c	dr=55.0d0
c	dr = 5.0d0

c	do 110 i=1,550
	
c	dr=dr-0.2d0
c	alpha=jj-1
c	alpha=jj

c	dxx=dx*dcos(alpha*pi/180.0d0)
c	dyy=dy*dsin(alpha*pi/180.0d0)
	

czapril	dzz=0.0d0+z
c	dzz=dzz+z
c	write(*,*)dxx,dyy,dzz

CZ	dx=dxx-x0
	
CZ	dy=dyy-y0

CZ	dz=dzz+z0

CZ	dz1=dzz-z0

	DX =DBLE(XGZAI)

      DY =DBLE(YETA)

	DZ =DBLE(ZZ)

      DZ1 =DBLE(ZZ1)



      CALL GREEN
czam the next portion is blocked by zaman------

      SGR = GRENR(1)

      SGI = GRENI(1)

	SXGR= GRENR(2)

	SXGI= GRENI(2)
	
	SYGR=GRENR(3)
		
	SYGI=GRENI(3)

	SZGR=GRENR(4)
      	
	SZGI= GRENI(4)

	

	r1=dsqrt(dx*dx+dy*dy+dz1*dz1)
	r2=dsqrt(dx*dx+dy*dy+dz*dz)

	grr=(1/r1)-(1/r2)+sgr
	gri=sgi


c	Calculation start for derivatives....
	
	rx1=r1**3
	rx2=r2**3
C	

	sxr=(-dx/rx1)-(-dx/rx2)+sxgr
	sxi=sxgi

	ry1=r1**3
	ry2=r2**3
	

	syr=(-dy/ry1)-(-dy/ry2)+sygr
	syi=sygi

	rz1=r1**3
	rz2=r2**3


	szr=(-dz1/rz1)-(-dz/rz2)+szgr
	szi=szgi


c	write(3,10)domega,speed,dx,dy,dz,grr,gri
c	write(3,10)domega,speed,dx,dy,dr,grr,gri
c	write(3,10)domega,speed,alpha,dx,dy,dr,SXR,SXI
c	write(5,10)domega,speed,alpha,dx,dy,dr,SYR,SYI
c	write(6,10)domega,speed,alpha,dx,dy,dr,SZR,SZI
c	write(3,10)domega,speed,alpha,dr,grr,gri
c	write(3,10)speed,dxx,dyy,x0,y0,grr,gri
c	write(*,*)dxx,dyy,grr,gri,sgr,sgi


c	write(jj,12)i,alpha,domega,speed,dr,grr,gri


   10 format(2f6.3,2x,3(2x,f8.4),2f25.15,x)
   11	format(6(f8.3,2x))
c   10	format(i6,2(x,f8.4),2f25.15)
   12 format(i6,3(2x,f7.3),f10.6,x,2f25.15,x)

	

	
110   continue
c	write(3,*)
2003	continue
		
	go to 2

1	write(*,*)' bad data, please check the value of omega'
	continue


2	return

      END

C

C
C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

C   Name    -- Subroutine GREEN

C   Purpose -- Calculation of the velocity Potential

C              of a Translating,Pulsating Source ( Wave Term only )

C   Notes   -- Double Precision

C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      SUBROUTINE GREEN
	
      IMPLICIT REAL*8 (A-H,O-Z)

      COMPLEX*16 ICONST

      COMMON/PARAMT/ GRAV,DOMEGA,SPEED,DX,DY,DZ,GRBATA

      COMMON/CONST / ICONST,PI,angle

      COMMON/INTGRL/ GAMMA,POLE1,POLE2

      COMMON/AUTPMT/ SRBASE(4,32),SIBASE(4,32),BOTTOM,PERIOD,EPS,IPHASE

     *              ,NSECT

      COMMON/COUNT/JDIV

      COMMON/SELECT/ISEL1

      COMMON/RESULT/ GRENR(4),GRENI(4)
czam-----
c	dimension SRBASE(4,32),SIBASE(4,32)

	
	 PI    = 3.141592653589793238D0
c	theta=angle*pi/180.0d0

      ICONST= (0.D0,1.D0)

      DATA PI2 /1.570796326794896619D0/

C---Initial Value Set

C-*-Definition of Values ------------------------------

C     DOMEGA    - Frequency of Oscillation (rad/s)

C     SPEED    - Steady Forward Speed of Source (m/s)

C------------------------------------------------------

     
C---Information of Cal.Condition
c	write(*,*)' start of green'

      GRBATA =DOMEGA*SPEED/GRAV

      IF(GRBATA.LT.0.25D0) GO TO 1

      GAMMA=DACOS(1.D0/(GRBATA*4.D0))

      GO TO 2

    1 GAMMA=0.D0

    2 CONTINUE

C---Quadrature of Integrand

        JDIV=0

        DO 20 IDERIV=1,4

          GRENR(IDERIV)=0.D0

          GRENI(IDERIV)=0.D0

   20   CONTINUE

        IF (GRBATA.GE.0.25D0) GO TO 30

C-1)-for Bata < 0.25 -----------------------------------------

         CALL QUADR(0.D0,PI,1,0.04D0,8)

         GO TO 40

C-2)-for Bata > 0.25 -----------------------------------------

   30    CONTINUE

         POT1=GAMMA

         POT2=GAMMA*0.95D0

         POT3=GAMMA+(PI2-GAMMA)*0.05D0

         A=POT2

         B=GAMMA-POT2

         C=POT3-GAMMA

         D=PI2-POT3

         CALL QUADR(0.D0,A,2,0.04D0,4)

         CALL QUADR(POT2,B,2,0.04D0,8)

         CALL QUADR(POT1,C,1,0.04D0,8)

         CALL QUADR(POT3,D,1,0.04D0,4)

         CALL QUADR(PI2,PI2,1,0.04D0,4)

   40 CONTINUE

C---Normal Return

      RETURN
	

      END

C






C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

C   Name    -- Subroutine QUADR

C   Purpose -- Quadrature of Theta range by Gauss-Legendre method.

C   Notes   -- Double Precision

C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      SUBROUTINE QUADR(BTM,PRD,IS1,EPSB,NS)

      IMPLICIT REAL*8 (A-H,O-Z)

      COMMON/PARAMT/ GRAV,DOMEGA,SPEED,DX,DY,DZ,GRBATA

      COMMON/AUTPMT/ SRBASE(4,32),SIBASE(4,32),BOTTOM,PERIOD,EPS,IPHASE

     *              ,NSECT

      COMMON/SELECT/ ISEL1

      COMMON/RESULT/ GRENR(4),GRENI(4)

      DIMENSION RESR(4),RESI(4),SECSR(4),SECSI(4)


czam--------------------


C---

      DO 10 IDERIV=1,4

        RESR(IDERIV)=0.D0

        RESI(IDERIV)=0.D0

   10 CONTINUE

C---Phase1 - The Period is divided into subinterval.

C            Max number of the partitions is 8.

      ISEL1 =IS1

      NSECT =NS

      IPHASE=1

      ISECT =1

      BOTTOM=BTM

      PERIOD=PRD

      EPS   =EPSB

      CALL AUTDIV( ISECT,SECSR,SECSI )

      IF ( IPHASE .EQ. 2 ) GO TO 40

        DO 30 IDERIV=1,4

          RESR(IDERIV)=SECSR(IDERIV)

          RESI(IDERIV)=SECSI(IDERIV)

   30   CONTINUE

        GO TO 60

C---Phase2 - The section is divided into subintervals.

C            Max number of the partitions is 32.

   40   PERIOD=PERIOD/NSECT

C       EPS   =EPS/NSECT

      DO 50 ISECT=1,NSECT

        CALL AUTDIV( ISECT,SECSR,SECSI )

        DO 50 IDERIV=1,4

          RESR(IDERIV)=RESR(IDERIV)+SECSR(IDERIV)

          RESI(IDERIV)=RESI(IDERIV)+SECSI(IDERIV)

   50   CONTINUE

C---Normal Return
czam--
c	write(4,*)'this is real and imag part of G'

   60   DO 70 IDERIV=1,4

          GRENR(IDERIV)=GRENR(IDERIV)+RESR(IDERIV)

          GRENI(IDERIV)=GRENI(IDERIV)+RESI(IDERIV)
czam-------check it 02/7/5 at 8.05pm

	
c	write(*,*)GRENR(IDERIV),GRENI(IDERIV)
c	write(*,*)dx,GRENR(3),GRENI(3)
czam------

   70   CONTINUE

      RETURN


      END


C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

C   Name    -- Subroutine AUTDIV

C   Purpose -- Divide into subintervals by automatic processing.

C              Quadrate Integrand by Gauss-Legendre Method.

C   Notes   -- Double Precision

C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      SUBROUTINE AUTDIV( ISECT,SECSR,SECSI )

      IMPLICIT REAL*8(A-H,O-Z)

      COMMON/AUTPMT/ SRBASE(4,32),SIBASE(4,32),BOTTOM,PERIOD,EPS,IPHASE

     *              ,NSECT

      COMMON/SELECT/ ISEL1

      DIMENSION  SECSR(4),SECSI(4)

     *          ,DIVSR(4,32),DIVSI(4,32),SRERR(4),SIERR(4)

c	write(6,*)'2,start of auto div'		

      GO TO ( 1,2 ) IPHASE

C---Preparing Process for Phase1

    1 IDIV=1

      TOP=BOTTOM+PERIOD

      CALL GAUSLD(BOTTOM,TOP,DIVSR,DIVSI,IDIV)

      DO 5  IDERIV=1,4

        SRBASE(IDERIV,ISECT)=DIVSR(IDERIV,IDIV)

        SIBASE(IDERIV,ISECT)=DIVSI(IDERIV,IDIV)

    5 CONTINUE

C---Main Process.

    2 NDIV=2

      H=PERIOD/2.D0

   10 DO 20 IDERIV=1,4

        SECSR(IDERIV)=0.D0

        SECSI(IDERIV)=0.D0

   20 CONTINUE

      UNDER=BOTTOM+PERIOD*(ISECT-1)

      OVER =UNDER

      DO 30 IDIV=1,NDIV

        UNDER=OVER

        OVER =UNDER+H

        CALL GAUSLD(UNDER,OVER,DIVSR,DIVSI,IDIV)

        DO 40 IDERIV=1,4

          SECSR(IDERIV)=SECSR(IDERIV)+DIVSR(IDERIV,IDIV)

          SECSI(IDERIV)=SECSI(IDERIV)+DIVSI(IDERIV,IDIV)

   40   CONTINUE

   30 CONTINUE

      DO 50 IDERIV=1,4

        IF ((SRBASE(IDERIV,ISECT).EQ.0.D0).AND.(SECSR(IDERIV).EQ.0.D0))

     *      GO TO 50

        IF ((SRBASE(IDERIV,ISECT).EQ.0.D0).AND.(SECSR(IDERIV).NE.0.D0))

     *      GO TO 60

        IF ((SIBASE(IDERIV,ISECT).EQ.0.D0).AND.(SECSI(IDERIV).EQ.0.D0))

     *      GO TO 50

        IF ((SIBASE(IDERIV,ISECT).EQ.0.D0).AND.(SECSI(IDERIV).NE.0.D0))

     *      GO TO 60

        SRERR(IDERIV)=DABS( SECSR(IDERIV)/SRBASE(IDERIV,ISECT) )

        SIERR(IDERIV)=DABS( SECSI(IDERIV)/SIBASE(IDERIV,ISECT) )

        SRERR(IDERIV)=DABS( SRERR(IDERIV)-1.D0 )

        SIERR(IDERIV)=DABS( SIERR(IDERIV)-1.D0 )

        IF ( SRERR(IDERIV) .GT. EPS )  GO TO 60

        IF ( SIERR(IDERIV) .GT. EPS )  GO TO 60

   50 CONTINUE

      GO TO 130

   60 CONTINUE

      GO TO ( 70,80 ) IPHASE

   70 IF ( NDIV .EQ. NSECT) GO TO 100

   80 IF ( NDIV .EQ. 32 )   GO TO 120

      NDIV=NDIV*2

      H=H/2.D0

      DO 90 IDERIV=1,4

        SRBASE(IDERIV,ISECT)=SECSR(IDERIV)

        SIBASE(IDERIV,ISECT)=SECSI(IDERIV)
czam----------	
c	write(4,*)SRBASE(IDERIV,ISECT),SIBASE(IDERIV,ISECT)
czam----------

   90 CONTINUE

      GO TO 10

C---Result of Integral don't satisfy the prescribed precision.

C    Preparing Process for Phase2

  100 IPHASE=2
c	write(4,*)'SRBASE(IDERIV,ISECT),SIBASE(IDERIV,ISECT)'


      DO 110 IDIV=1,NSECT

        ISECT=IDIV

        DO 110 IDERIV=1,4

          SRBASE(IDERIV,ISECT)=DIVSR(IDERIV,IDIV)

          SIBASE(IDERIV,ISECT)=DIVSI(IDERIV,IDIV)
c		write(4,*)SRBASE(IDERIV,ISECT),SIBASE(IDERIV,ISECT)


  110 CONTINUE

      GO TO 130

C---Error Message

  120 CONTINUE

C     WRITE(*,*) '+++ Warning : Result of Integral ',

C    *           'didn''t satisfy the prescribed precision. +++'

C    *           ,ISEL1

      DO 125 IDERIV=1,4

C       WRITE(*,*) SRERR(IDERIV),SIERR(IDERIV)

  125 CONTINUE

C---Normal Return

  130 RETURN

      END


	

c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

C   Name    -- Subroutine GAUSLD

C   Purpose -- Quadrature of Integrand on Theta-range

C              by Gauss-Legendre Method

C   Notes   -- Double Precision

C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      SUBROUTINE GAUSLD(A,B,SUMR,SUMI,IDIV)

      IMPLICIT REAL*8(A-H,O-Z)

      COMPLEX*16 GRIDC,GRID

      COMMON /COUNT / JDIV

      DIMENSION SUMR(4,32),SUMI(4,32),GRID(4),X(4),Y(4),U(8),AA(8)

C---Set Values of Weight

      DATA X / 0.96028 98564 97536 20D0,0.79666 64774 13626 70D0,

     *         0.52553 24099 16329 00D0,0.18343 46424 95649 80D0 /

      DATA Y / 0.10122 85362 90376 30D0,0.22238 10344 53374 50D0,

     *         0.31370 66458 77887 30D0,0.36268 37833 78362 00D0 /




C---Cal. Process

      C1=(B+A)/2.D0

      C2=(B-A)/2.D0

      DO 1 J=1,4

        U(J)  =C1-C2*X(J)

        U(9-J)=C1+C2*X(J)

        AA(J)  =Y(J)

        AA(9-J)=Y(J)

    1 CONTINUE

      DO 2 IDERIV=1,4

        SUMR(IDERIV,IDIV)=0.D0

        SUMI(IDERIV,IDIV)=0.D0

    2 CONTINUE

      DO 3 J=1,8
CZAM	DO 3 J=1,1
	
        VAL=U(J)

        CALL GREENF (VAL,GRID)

       DO 3 IDERIV=1,4
        
         GRIDC=GRID(IDERIV)

          GRIDR=DREAL(GRIDC)

          GRIDI=DIMAG(GRIDC)

          SUMR(IDERIV,IDIV)=SUMR(IDERIV,IDIV)+AA(J)*GRIDR

          SUMI(IDERIV,IDIV)=SUMI(IDERIV,IDIV)+AA(J)*GRIDI

    3   CONTINUE
c      write(4,*)'SUMR(IDERIV,IDIV),SUMI(IDERIV,IDIV)'

      DO 4 IDERIV=1,4

        SUMR(IDERIV,IDIV)=SUMR(IDERIV,IDIV)*C2

        SUMI(IDERIV,IDIV)=SUMI(IDERIV,IDIV)*C2

c	write(4,*)SUMR(IDERIV,IDIV),SUMI(IDERIV,IDIV)

    4 CONTINUE

C---Normal Return

      JDIV=JDIV+1

      RETURN

      END

C

C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

C   Name    -- Subroutine GREENF

C   Purpose -- Cal. Integrand of Theta range Quadrature

C   Notes   -- Double Precision

C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      SUBROUTINE GREENF (THETA,GRID)

      IMPLICIT REAL*8 (A-H,O-Z)

      COMPLEX*16 GRID,GRIDPV,GRIDDT,QWBAS1,QWBAS2,QWNUM1,QWNUM2,QWNUM3,

     *           QWDEN,HBAS1,HBAS2,HBAS3,K1,K2,FKJ,ICONST,part3,

     *           H1,H2,H3,QW1,QW2

      COMMON/PARAMT/ GRAV,DOMEGA,SPEED,DX,DY,DZ,GRBATA

      COMMON/CONST / ICONST,PI

      COMMON/SELECT/ISEL1

      COMMON/INTGRL/ GAMMA

      DIMENSION GRID(4),GRIDPV(4),GRIDDT(4),

     *          QWBAS1(4),QWBAS2(4),HBAS1(2),HBAS2(2),HBAS3(2)

C---Cal.Process

      PART1=DCOS(THETA)

      PART2=DSIN(THETA)

      PART3=1.D0-4.D0*GRBATA*PART1

      PART3=CDSQRT(PART3)



	K1=FKJ(1,THETA)

      K2=FKJ(2,THETA)

C-*-Cal.selecter set

      JSEL=0

        IF (ISEL1.EQ.2) THEN

         IF (THETA.LE.GAMMA) JSEL=1

        END IF

C-<1>-Principal Value Integral Part of Green4

        WPART1=DX*PART1

        WPART2=DY*PART2

        W1    =WPART1+WPART2

        W2    =WPART1-WPART2

        CALL QWBASE(K1,W1,QW1,QW2,JSEL)

             QWBAS1(1)=QW1

             QWBAS2(1)=QW2

        CALL QWBASE(K1,W2,QW1,QW2,JSEL)

             QWBAS1(2)=QW1

             QWBAS2(2)=QW2

        CALL QWBASE(K2,W1,QW1,QW2,JSEL)

             QWBAS1(3)=QW1

             QWBAS2(3)=QW2

        CALL QWBASE(K2,W2,QW1,QW2,JSEL)

             QWBAS1(4)=QW1

             QWBAS2(4)=QW2

        QWNUM1=QWBAS1(1)+QWBAS1(2)-QWBAS1(3)-QWBAS1(4)

        QWNUM2=QWBAS2(1)+QWBAS2(2)-QWBAS2(3)-QWBAS2(4)

        QWNUM3=QWBAS2(1)-QWBAS2(2)-QWBAS2(3)+QWBAS2(4)
c	QWNUM3=QWBAS2(1)+QWBAS2(2)-QWBAS2(3)-QWBAS2(4)


        QWDEN =PART3*PI

        GRIDPV(1)=             QWNUM1/QWDEN

        GRIDPV(2)=ICONST*PART1*QWNUM2/QWDEN

        GRIDPV(3)=ICONST*PART2*QWNUM3/QWDEN

        GRIDPV(4)=             QWNUM2/QWDEN

        IF  (JSEL.EQ.1) THEN

          DO 10 IDERIV=1,4

           GRIDDT(IDERIV)=(0.D0,0.D0)

   10     CONTINUE

          GO TO 100

        END IF

C-<2>-Detour Integral Part of Green4

        CALL HBASE (K1,THETA,H1,H2,H3)

             HBAS1(1)=H1

             HBAS2(1)=H2

             HBAS3(1)=H3

        CALL HBASE (K2,THETA,H1,H2,H3)

             HBAS1(2)=H1

             HBAS2(2)=H2

             HBAS3(2)=H3

        IF ( THETA .GT. 0.5D0*PI ) GO TO 1

          GRIDDT(1)= 2.D0*ICONST      *(HBAS1(1)+HBAS1(2))/PART3

          GRIDDT(2)=-2.D0       *PART1*(HBAS2(1)+HBAS2(2))/PART3

          GRIDDT(3)=-2.D0*ICONST*PART2*(HBAS3(1)+HBAS3(2))/PART3
czam	GRIDDT(3)=2.D0*ICONST*PART2*(HBAS3(1)+HBAS3(2))/PART3


          GRIDDT(4)= 2.D0*ICONST      *(HBAS2(1)+HBAS2(2))/PART3

            GO TO 2

    1     GRIDDT(1)= 2.D0*ICONST      *(HBAS1(1)-HBAS1(2))/PART3

          GRIDDT(2)=-2.D0       *PART1*(HBAS2(1)-HBAS2(2))/PART3

          GRIDDT(3)=-2.D0*ICONST*PART2*(HBAS3(1)-HBAS3(2))/PART3
czam	GRIDDT(3)=2.D0*ICONST*PART2*(HBAS3(1)-HBAS3(2))/PART3


          GRIDDT(4)= 2.D0*ICONST      *(HBAS2(1)-HBAS2(2))/PART3

    2   CONTINUE

C-<3>-Integrand Green Function

  100   DO 3 IDERIV=1,4

          GRID(IDERIV)=GRIDPV(IDERIV)+GRIDDT(IDERIV)
czam	write(*,*) GRID(IDERIV)
    3   CONTINUE

C---Normal Return

      RETURN
      END

C

C

C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

C   Name    -- Subroutine QWBASE

C   Purpose -- Cal. Basic Part of Q(w) in Integrand ( Theta range )

C   Notes   -- Double Precision

C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      SUBROUTINE QWBASE (KJ,W,QWBAS1,QWBAS2,JSEL)

      IMPLICIT REAL*8(A-H,O-Z)

      COMPLEX*16 QWBAS1,QWBAS2,KJ,IJ,QJ,ICONST,EXPITG,

     *           PART1,PART2,PART3,PART4

      COMMON/PARAMT/ GRAV,DOMEGA,SPEED,DX,DY,DZ,GRBATA

      COMMON/CONST / ICONST,PI
c      write(6,*)'4,start of qwbase'


C---Cal.Process

      IF (JSEL.EQ.1) GO TO 200

C-*-Cal. for Type 1

      QJ=-KJ*(DABS(DZ)-ICONST*W)

      IF (W) 10,20,30

   10  SGNW=-1.D0

         GO TO 40

   20  SGNW= 0.D0

         GO TO 40

   30  SGNW= 1.D0

   40 CONTINUE

      PART1=CDEXP(QJ)

      PART2=ICONST*PI*PART1*SGNW

      PART3=1.D0/QJ

      PART4=EXPITG(QJ)

      IJ   =PART2+PART4-PART3

      GO TO 300

C---Cal. for Type 2

  200 RKJR=DREAL(KJ)

      RKJI=DIMAG(KJ)

      DZABS=DABS(DZ)

      QJR =-RKJI*W-RKJR*DZABS

      QJI = RKJR*W-RKJI*DZABS

      ZETAR= DZABS

      ZETAI= W

      ARG1 =DATAN2(ZETAI,ZETAR)

      ARG2 =DATAN2(RKJI,RKJR)

      QJ=-KJ*(DABS(DZ)-ICONST*W)

      IJ=EXPITG(QJ)-1.D0/QJ

C----

      IF (W) 50,60,70

   50 SGNW=-1.D0

        IF ((0.D0.GT.ARG2).AND.(ARG2.GT.ARG1)) GO TO 80

           GO TO 300

   60 SGNW= 0.D0

           GO TO 300

   70 SGNW=+1.D0

        IF ((ARG1.GT.ARG2).AND.(ARG2.GT.0.D0)) GO TO 80

           GO TO 300

   80 PART1=CDEXP(QJ)

      PART2=2.D0*ICONST*PI*PART1*SGNW

      IJ=IJ+PART2

C---Conclusion

  300 QWBAS1=KJ*IJ

      QWBAS2=KJ*KJ*IJ

C---Normal Return

      RETURN

      END

C

C

C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

C   Name    -- Subroutine HBASE

C   Purpose -- Cal. Basic Part of h(theta,k) in Integrand.

C   Notes   -- Double Precision

C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      SUBROUTINE HBASE (KJ,THETA,HBAS1,HBAS2,HBAS3)

      IMPLICIT REAL*8(A-H,O-Z)

      COMPLEX*16 KJ,PART3,PART4,PART5,PART6,HBAS1,HBAS2,HBAS3

      COMMON/PARAMT/ GRAV,DOMEGA,SPEED,DX,DY,DZ,GRBATA

C---Cal.Process of Green Function H
c	write(6,*)'5,start of sub hbase'

      PART1 =DCOS(THETA)

      PART2 =DSIN(THETA)

      PART3R=KJ*DZ

      PART3I=KJ*DX*PART1

      PART3 =DCMPLX(PART3R,PART3I)

      PART3 =CDEXP(PART3)

      PART4 =KJ*DY*PART2

      PART5 =CDCOS(PART4)

      PART6 =CDSIN(PART4)

      HBAS1 =   KJ*PART3*PART5

      HBAS2 =KJ*KJ*PART3*PART5

      HBAS3 =KJ*KJ*PART3*PART6


C---Normal Return

      RETURN

      END

C

C

C

C   Name    -- Function FKJ
c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

C   Name    -- Function EXPITG

C   Purpose -- Cal. exp(Qj) * E(Qj)

C                 E(x): Complex Exponential Integral

C   Notes   -- Double Precision

C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	
      COMPLEX FUNCTION EXPITG*16 (QJ)

      IMPLICIT REAL*8 (A-H,O-Z)

      COMPLEX*16 QJ,TERM1,TERM2,TBASE,EXITG1,EXITG2,ICONST

      COMMON/CONST/ICONST

C---Euler Constant +-----+-----+-----+
c	write(6,*)'6,start of k1,k2'


      DATA EULER /0.57721 56649 01532 9D0/

C---Cal. Process

      JCVG= 0

      QJR= DREAL(QJ)

      QJI= DIMAG(QJ)

      QJR2= QJR**2

      QJI2= QJI**2

      QJABS = DSQRT(QJR2+QJI2)

      NLIMIT=IDINT(QJABS)

      IF ( QJABS .LE.  12.D0 ) GO TO 10

      IF ( QJABS .GT.  12.D0 ) GO TO 20

C-<1>-Complete Expansion.

   10  KCAL=1

       EXITG1=-QJ

       TERM2 =-QJ

       N=2

   11  TERM2 = TERM2*(-QJ)/N

       TERM1 = TERM2/N

       EXITG2=EXITG1+TERM1

       GO TO 30

   12  EXPITG= EXITG2+EULER+CDLOG(QJ)

       EXPITG=-EXPITG

       EXPITG= EXPITG*CDEXP(QJ)

       GO TO 40

C-<2>-Asymptonic Expansion

   20  KCAL=2

       EXITG1=1.D0

       TBASE =1.D0/(-QJ)

       TERM1 =1.D0

       N=1

   21  TERM1=TERM1*TBASE*N

       EXITG2=EXITG1+TERM1

       IF ( N .EQ. NLIMIT ) GO TO 22

       GO TO 30

   22  EXPITG=EXITG2/QJ

       GO TO 40

C-<*>-Check of Convergence

   30  IF ( EXITG2 .EQ. EXITG1 ) JCVG=JCVG+1

       IF ( JCVG .EQ. 2 )    GO TO (12,22),KCAL

         EXITG1=EXITG2

         N=N+1

         GO TO (11,21),KCAL

   40 IF ( QJI .EQ. 0.D0 ) EXPITG=DREAL(EXPITG)

C---Normal Return

      RETURN

      END

C

C

C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C   Purpose -- Cal. Singular Point-k1,k2

C   Notes   -- Double Precision

C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
      COMPLEX FUNCTION FKJ*16 (J,THETA)

      IMPLICIT REAL*8(A-H,O-Z)

      COMPLEX*16 PART3,PART4,PART41,PART42

      COMMON/PARAMT/ GRAV,DOMEGA,SPEED,DX,DY,DZ,GRBATA

C---Cal.Process
c	write(6,*)'7,start of exp ie "fkj"'


      PART1=DCOS(THETA)

      PART2=(DOMEGA**2)/GRAV

      PART3=1.D0-4.D0*GRBATA*PART1
      PART3=CDSQRT(PART3)

      PART41=1.D0+PART3
	PART42=1.D0-PART3
	
      PART5=2.D0*GRBATA*PART1

      GO TO (1,2),J

c    1 FKJ=4.D0*PART2/(PART41**2)
    1 FKJ=PART2*((PART42/PART5)**2)
      GO TO 3

    2 FKJ=PART2*((PART41/PART5)**2)

C---Normal Return

    3 CONTINUE

      RETURN

      END

C