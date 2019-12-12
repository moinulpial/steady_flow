C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

C   Purpose --There are 7 programs under this folder
c              
c		      1. ggfile.f (reads input data and generates output 
c							for the following programs...)

c	    	  2. mggren.f (calculates the Green function (G,Gn) for 
c							the body surface)

c			  3. svgren_surf.f (calculates the source densities(sigma) 
c								for the body surface by using the G,Gn
c								 (output of mggren.f))

c			  4. surface_grn.f(	calculates the Green function for 
c							      the free surface)

c			  5. phi.f	(calculates the velocity potential for the free 
c						surface and free surface elevation by using source 
c						dencities(output of svgren_surf.f) and free surface 
c						Green function(output of surface_grn.f))

c			  6. gnu.f  ( generates output for  ploting graph)

c			  7. graphplot.f( generates output  according 
c													to "gnuplot" software)   	

C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
c*********************************************************************************************    
C																						   *	
c	THIS IS THE SECOND PROGRAM FOR THE CALCULATION  OF STEADY FREE SURFACE FLOW PATTERN    *
C																						   *	
C	       THIS PROGRAM TAKES INPUT FROM GGFILE.F AND CALCULATES GREEN FUNCTION            *
C																						   *		
C		            	IT GENERATES INPUT FOR  SVGRN_SURF.F							   *	
c*********************************************************************************************
c   
C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

C   Purpose -- Cal. Coefficient matrix (G,Gn)

C   Notes   --1.the value of Singular Part ( Term of Rankin's

C               Source Potential ) is cal.ed by the Hess & Smith 's

C               method.

C             2.the value of Non Singular Part ( Wave Term Potential )

C               is calculated about each Panel's Centroid.

C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

c----------------------------------------------------------------

C     ALFA - Normal-derivatine of Green coefficient

C     BETA - Green coefficient

C	BETA2- X DERIVATIVE OF GREEN COEFFICIENT

C----------------------------------------------------------------

	PARAMETER (PAR=1000)

      COMPLEX ALFA,BETA,BETA2,AB,XYZN,aalfa

      COMMON /FILE / KEY(72),ITI,ITO,IFILE1,IFILE4,IFILE2,IFILE3

     *              ,IFILE5,IFILE6,IFILE7,IFILE8,IFILE9

      COMMON /IWAVE/ RK,V,H,W(9),ISYM(2)

      COMMON /WAVE / NTW,TWAVE(5),IWW

      COMMON /PANEL/XVALUE(4,PAR),YVALUE(4,PAR),ZVALUE(4,PAR)

      COMMON /MUGEN/ INF


	DIMENSION XG(PAR),YG(PAR),ZG(PAR),XYZN(PAR,3),X(par),Y(par),
	
     *	Z(par),XLENG(PAR),YLENG(PAR),zleng(PAR)

     *    ,PGRR(4),PGRI(4),PGRDNR(4),PGRDNI(4),PGRDXR(4),PGRDXI(4)

     *     ,CXI(5),CETA(5),TR(9),UMU(50),CA(50),ENXYZ(3),DMY(5),

     *	  BETA(par,par),ALFA(par,par),BETA2(par,par)

            


	OPEN(1,FILE='output/centroid',STATUS='UNKNOWN')
	OPEN(2,FILE='output/ggfile',STATUS='UNKNOWN')

	OPEN(4,FILE='output/alfa_sig',STATUS='UNKNOWN')
	OPEN(5,FILE='output/beta_sig',STATUS='UNKNOWN')
	OPEN(6,FILE='output/beta2',STATUS='UNKNOWN')

	open(7,file='output/enxyz',status='unknown')
	open(8,file='output/length',status='unknown')

c	write(*,*)'nelem = ?'
c	read(*,*) nelem

	write(*,*)'please check subroutine grtrm2' 
	write(*,*)'before running this program '
	write(*,*)
	write(*,*)'nelem = ?'
	read(*,*) nelem
	ndiv=nelem
	idiv=1


C---Parameter Set.

      DATA PAI /3.14159265358979/

C---Parameter Initialize

      DO 1 I=1,NDIV

        DO 1 J=1,NDIV

          ALFA(I,J) =(0.0,0.0)

          BETA(I,J) =(0.0,0.0)

          BETA2(I,J)=(0.0,0.0)

    1 CONTINUE

      FAI   =0.0

      FAIDN =0.0

      FAIDX =0.0

      GRR   =0.0

      GRI   =0.0

      GRDNR =0.0

      GRDNI =0.0

      GRDXR =0.0

      GRDXI =0.0

      CGRR   =0.0

      CGRI   =0.0

      CGRDNR =0.0

      CGRDNI =0.0

      CGRDXR =0.0

      CGRDXI =0.0

      DO 50 J1=1,4

        PGRR(J1)  =0.0

        PGRI(J1)  =0.0

        PGRDNR(J1)=0.0

        PGRDNI(J1)=0.0

        PGRDXR(J1)=0.0

        PGRDXI(J1)=0.0

   50 CONTINUE




	DO 300 I=1,nelem
	read(1,*)xg(i),yg(i),zg(i)
	READ(7,11) XYZN(I,1),XYZN(I,2),XYZN(I,3)
	READ(8,*)YLENG(I)
 11    format(3(2f12.8))
 12    format(3f12.8)
300	continue	



        DO 1250 M=1,NDIV


	write(*,*)'calculating for element ',m
c	write(*,*)

	READ(2,*) II,(CXI(J),J=1,4),(CETA(J),J=1,4),(TR(J),J=1,9)
     *                ,SS,XXI,XYI,YYI,EL,HMAX
	

c          K = (IJ-1) * NDIV + M
				K = M

          DO 1200 L=1,NDIV

            DO 450 KKK=1,3
	
              ENXYZ(KKK) = XYZN(L,KKK)

  450       CONTINUE

c	 write(*,*)' Element',m, 'and',l

C

C---Calculation of Green function

C-*-Cal. Term1

      CALL GRTRM1(XG(L),YG(L),ZG(L),XG(K),YG(K),ZG(K),TR,EL,SS


     *            ,ENXYZ,XXI,XYI,YYI,CXI,CETA,FAI,FAIDN,FAIDX)



C-*-Cal. Term2

C--1)-Cal. Potential of the Panel's Centroid.

      XGZAI=XG(L)-XG(K)

      YETA =YG(L)-YG(K)

      ZZ1  =ZG(L)+ZG(K)




      CALL GRTRM2(m,yleng(m),XGZAI,YETA,ZZ1,ENXYZ,SS,CGRR,CGRI

     *            ,CGRDNR,CGRDNI,CGRDXR,CGRDXI)

C--2)-Cal. Potential of the Panel's 4 nodes.

      GO TO 510


C--3)-Mean Value of Green Function on the Panel.

  510 GRR  =CGRR

      GRI  =CGRI

      GRDNR=CGRDNR

      GRDNI=CGRDNI

      GRDXR=CGRDXR

      GRDXI=CGRDXI

c      IF ( KEY(70).EQ.0 ) GO TO 560

	GO TO 560

      DO 550 J3=1,4

        GRR  = GRR  +PGRR(J3)

        GRI  = GRI  +PGRI(J3)

        GRDNR= GRDNR+PGRDNR(J3)

        GRDNI= GRDNI+PGRDNI(J3)

        GRDXR= GRDXR+PGRDXR(J3)

        GRDXI= GRDXI+PGRDXI(J3)

  550 CONTINUE

        GRR  = GRR  / 5.0

        GRI  = GRI  / 5.0

        GRDNR= GRDNR/ 5.0

        GRDNI= GRDNI/ 5.0

        GRDXR= GRDXR/ 5.0

        GRDXI= GRDXI/ 5.0

C---Green Function Matrix

  560 IF (L.EQ.K) FAIDX=FAIDX+0.5*ENXYZ(1)

      AL =FAIDN+GRDNR

      BE =FAI  +GRR

      BE2=FAIDX+GRDXR

      BETA (L,M)=CMPLX(BE ,GRI  )

      BETA2(L,M)=CMPLX(BE2,GRDXI)

      ALFA (L,M)=CMPLX(AL ,GRDNI)

  850 IF(L.NE.K) GO TO 1450

      AB=(0.5,0.0)

 1100 ALFA(L,M)=AB+ALFA(L,M)

c	WRITE(4,*) ALFA (L,M)

c      WRITE(5,*) BETA (L,M)



c
 1450 CONTINUE



 1200 CONTINUE
c	write(*,*)'ok'
	

 1250 CONTINUE

C---Results are written into Work Files

	WRITE(4,*) ((ALFA (I,J),i=1,NDIV),j=1,NDIV)

      WRITE(5,*) ((BETA (I,J),J=1,NDIV),I=1,NDIV)

      WRITE(6,*) ((BETA2(I,J),J=1,NDIV),I=1,NDIV)
		
		
113	format(4(2f16.10))

 5000 FORMAT(1H ,5X,5HERROR,10X,8HISYM(1)=,I5,3X,8HISYM(2)=,I5)

 5109 FORMAT(1H0,3X,8H* ALFA *,3X,1H(,I2,1H,,I2,8H) MATRIX / (4X

     1      ,4(2E12.4,2X))  )

 5209 FORMAT(1H0,3X,8H* BETA *,3X,1H(,I2,1H,,I2,8H) MATRIX  / (4X

     1      ,4(2E12.4,2X))   )

C---Normal Return

c      RETURN
	stop

      END

C

C

C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

C   Name    -- Subroutine GRTRM1

C   Purpose -- Cal. Term of Rankin's source potential.

C   Notes   -- Single Precision

C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      SUBROUTINE GRTRM1(XL,YL,ZL,XG,YG,ZG,TR,EL,SS,ENXYZ,XXI,XYI,YYI

     *                 ,CXI,CETA,FAI,FAIN,FAIDX)

C----------------------------------------------------------------

C     TR   - Transform Matrix ( direction cosine of Local axis )

C     SS   - Area of Element

C     FAI  - Potential of Source

C     FAIN - Partial derivative of potential with respect to Normal

C     ENXYZ- direction cosine of Normal with respect to (x,y,z)-coord.

C----------------------------------------------------------

      COMMON /FILE/ KEY(72),ITI,ITO

      COMMON /IWAVE/ RK,V,H

      DIMENSION  TR(9),ENXYZ(3),ZZ(2),CXI(5),CETA(5)

C-*-Parameter Initialize

      FAI  = 0.0

      FAXX = 0.0

      FAYY = 0.0

      FAZZ = 0.0

C-*-Parameter Set

      P4=-0.07957747

      CXI(5)  = CXI(1)

      CETA(5) = CETA(1)

      ZZ(1) =  ZL

      ZZ(2) = -ZL

C---Cal. Process

      DO 2000 I=1,2

C-*-Cal. Potential & derivative

        CALL RANKIN (XL,YL,ZZ(I),XG,YG,ZG,EL,SS,XXI,XYI,YYI,CXI,CETA

     *               ,TR,FAII,FAX,FAY,FAZ )

C-*-Transform (gzai,eya,zeta)-derivative to (x,y,z)-derivative

        FX = TR(1)*FAX + TR(4)*FAY + TR(7)*FAZ

        FY = TR(2)*FAX + TR(5)*FAY + TR(8)*FAZ

        FZ = TR(3)*FAX + TR(6)*FAY + TR(9)*FAZ

        IF (I.EQ.2) FZ = -FZ

C-*-Sum up Result

        IF (I.EQ.1) SIGN=+1.0

        IF (I.EQ.2) SIGN=-1.0

          FAI  = FAI  + SIGN*FAII

          FAXX = FAXX + SIGN*FX

          FAYY = FAYY + SIGN*FY

          FAZZ = FAZZ + SIGN*FZ

C------------

 2000 CONTINUE

        FAI = FAI*P4

C-*-Cal. Normal-derivative

        FAIDX= FAXX*P4

        FAIN = FAXX*ENXYZ(1) + FAYY*ENXYZ(2) + FAZZ*ENXYZ(3)

        FAIN = FAIN*P4

C---Cal. check No.66

      IF(KEY(66).LE.0) GO TO 9000

c        WRITE(ITO,100) FAI,FAIN,FAXX,FAYY,FAZZ

  100 FORMAT(1H0,'FAI= ',E12.5,3X,'FAIN=',E12.5,3X,'FAXX=',E12.5

     *      ,3X,'FAYY=',E12.5,3X,'FAZZ=',E12.5                     )

C---Normal Return

 9000 RETURN

      END

C

C

C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

C   Name    -- Subroutine GRTRM2

C   Purpose -- Cal. Wave Term Potential

C   Notes   -- Single precision

C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      SUBROUTINE GRTRM2(m,yleng,XGZAI,YETA,ZZ1,ENXYZ,SS,GRR,GRI,GRDNR,	

     *          GRDNI,GRDXR,GRDXI)

              

      REAL*8 GRAV,DOMEGA,SPEED,DX,DY,DZ,GRENR,GRENI,DUM

     *      ,SGR,SGI,SGNR,SGNI

      COMMON/FILE  / KEY(72),ITI,ITO

      COMMON/PARAMT/ GRAV,DOMEGA,SPEED,DX,DY,DZ,DUM

      COMMON/RESULT/ GRENR(4),GRENI(4)

      COMMON/IWAVE / DUMMY(3),WD(2)

      COMMON/WAVE  / NTW,TWAVE(5),IWW,FN,FORWRD

      DIMENSION  ENXYZ(3)

	
CC	NX= NUMBER OF ELEMENTS AROUND THE FIRST WATER LINE
CC	nd= number of water lines
cc	fn= froude number
cc	alg= length of the ship


	domega=.001
	grav=9.8
	fn=0.2670
	alg=2.0
	nx=50
	nd=5

	forwrd=fn*sqrt(grav*alg)


	
	
	 
C---Parameter set

      P4=-0.07957747

      DX =DBLE(XGZAI)

      DY =DBLE(YETA)

      DZ =DBLE(ZZ1)

c      GRAV =DBLE(WD(2))

c      DOMEGA=DBLE(WD(1))

      SPEED=DBLE(FORWRD)

      IF (SPEED.LE.0.005D0) SPEED=0.005D0

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

C     GRDNR - real part of Normal-derivative

C     GRDNI - imagnary part of Normal-derivative

C----------------------------------------------------------------

C---Cal. Potential & Derivative of unit Source

      CALL GREEN

      SGR = GRENR(1)

      SGI = GRENI(1)

      SGNR= GRENR(2)*ENXYZ(1) + GRENR(3)*ENXYZ(2) + GRENR(4)*ENXYZ(3)

      SGNI= GRENI(2)*ENXYZ(1) + GRENI(3)*ENXYZ(2) + GRENI(4)*ENXYZ(3)

C---Cal. Potential & Derivative of Panel

czam      GRR  = SNGL(SGR) *SS*P4

czam      GRI  = SNGL(SGI) *SS*P4

czam      GRDNR= SNGL(SGNR)*SS*P4

czam      GRDNI= SNGL(SGNI)*SS*P4

czam      GRDXR= SNGL(GRENR(2))*SS*P4

czam      GRDXI= SNGL(GRENI(2))*SS*P4



C---Cal. Potential & Derivative of Panel

c	dndphi(m,1)=ENXYZ(1)

c	N1=m+(nx*nd)
c	N2=m+(nx*nd)

C	if((m.eq.elwl).or.(m.eq.elwl1))

	if(m.gt.nx)THEN
	go to 100
		else
			go to 50
	ENDIF
	



c	if((m.le.nx))then

50	grr1=p4*sngl(sgr)*ss
	gri1=p4*sngl(sgi)*ss
	grr2=p4*forwrd**2/grav*sngl(sgr)*(enxyz(1))*yleng
	gri2=p4*forwrd**2/grav*sngl(sgi)*(enxyz(1))*yleng
c	write(*,*)m,grr2,gri2,yleng
	GRR  = grr1+grr2
      GRI  = gri1+gri2
c	write(*,*)m,grr2,gri2,yleng
	


	grdnr1=p4*sngl(sgnr)*ss
	grdni1=p4*sngl(sgni)*ss
	grdnr2=p4*forwrd**2/grav*sngl(sgnr)*(enxyz(1))*yleng
	grdni2=p4*forwrd**2/grav*sngl(sgni)*(enxyz(1))*yleng

	GRDNR= grdnr1+grdnr2
      GRDNI= grdni1+grdni2



	grdxr1=p4*sngl(GRENR(2))*ss
	grdxi1=p4*sngl(GRENI(2))*ss
	grdxr2=p4*forwrd**2/grav*sngl(GRENR(2))*(enxyz(1))*yleng
	grdxi2=p4*forwrd**2/grav*sngl(GRENI(2))*(enxyz(1))*yleng


	GRDXR=GRDXR1+GRDXR2 
      GRDXI=GRDXI1+GRDXI2 


	GO TO 9000






C	else
100	grr1=p4*sngl(sgr)*ss
	gri1=p4*sngl(sgi)*ss
	grr2=0.0
	gri2=0.0
C	endif
c	write(*,*)m,grr2,gri2

	GRR  = grr1+grr2
      GRI  = gri1+gri2


	grdnr1=p4*sngl(sgnr)*ss
	grdni1=p4*sngl(sgni)*ss
	grdnr2=0.0
	grdni2=0.0
C	endif

	GRDNR= grdnr1+grdnr2
      GRDNI= grdni1+grdni2



C	else
	grdxr1=p4*sngl(GRENR(2))*ss
	grdxi1=p4*sngl(GRENI(2))*ss
	grdxr2=0.0
	grdxi2=0.0
C	endif

	GRDXR=GRDXR1+GRDXR2 
      GRDXI=GRDXI1+GRDXI2 




 9000 RETURN

      END

C

C

C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

C   Name    -- Subroutine RANKIN

C   Purpose -- Cal. Potential & derivative of Rankin source.

C   Notes   -- Single Precision

C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      SUBROUTINE RANKIN(XL,YL,ZL,XG,YG,ZG,EL,SS,XXI,XYI,YYI,CXI,CETA

     *                 ,TR,FAI,FAX,FAY,FAZ)

C----------------------------------------------------------------

C     TR   - Transform Matrix ( direction cosine of Local axis )

C     EL   - Length of Element

C     SS   - Area of Element

C     XXI  - Second moment of Area

C     XYI  - Product of Inertia

C     YYI  - Second moment of Area

C     CXI  - Local Coordinate | gzai

C     CETA -  of Nodal Points | eta

C     FAI  - Potential of Source

C     FAIX - Partial        | with respect to gzai-axis

C     FAIY -  derivative    | with respect to  eta-axis

C     FAIZ -   of potential | with respect to zeta-axis

C----------------------------------------------------------

      COMMON /FILE/ KEY(72),ITI,ITO

      DIMENSION TR(9),CXI(5),CETA(5)

      DATA PAI / 3.14159265358979 /

C---Cal. Process

      XD = TR(1)*(XL-XG) + TR(2)*(YL-YG) + TR(3)*(ZL-ZG)

      YD = TR(4)*(XL-XG) + TR(5)*(YL-YG) + TR(6)*(ZL-ZG)

      ZD = TR(7)*(XL-XG) + TR(8)*(YL-YG) + TR(9)*(ZL-ZG)

      XX = XD*XD

      YY = YD*YD

      ZZ = ZD*ZD

C-*-Cal. check No.32

czzeded    IF(KEY(32).LE.0) GO TO 10

czzeked      WRITE(ITO,1000) XD,YD,ZD

czeked   10 WRITE(8,1001) XD,YD,ZD

	R0=SQRT(XX+YY+ZZ)

      IF(R0.LE.1.E-5) GO TO  20

      W=1.0/R0

      R03=R0*R0*R0

      R05=R03*R0*R0

      R03=1.0/R03

      R05=1.0/R05

      WX=-XD*R03

      WY=-YD*R03

      WZ=-ZD*R03

   20 RL=R0/EL

C-1)-Process for ( Length ratio RL > 4.0 )

      IF(RL.LT.4.0) GO TO 30

      FAI = SS*W

      FAX = SS*WX

      FAY = SS*WY

      FAZ = SS*WZ

      GO TO 160

C-2)-Process for ( 2.45 < RL < 4.0 )

   30 IF(RL.LT.2.45) GO TO 40

      RW=3.0*R05

      WXX=RW*XX-R03

      WXY=RW*XD*YD

      WYY=RW*YY-R03

      R07=R05*R03*R0

      P=YY+ZZ-4.0*XX

      Q=ZZ+XX-4.0*YY

C-*

      RW=3.0*R07

      WXXX=RW*XD*(3.0*P+10.0*XX)

      WXXY=RW*YD*P

      WXYY=RW*XD*Q

      WYYY=RW*YD*(3.0*Q+10.0*YY)

      WXXZ=RW*ZD*P

      WXYZ=-5.0*RW*XD*YD*ZD

      WYYZ=RW*ZD*Q

C-*

      FAI=SS*W+0.5*XXI*WXX+XYI*WXY+0.5*YYI*WYY

      FAX=SS*WX+0.5*XXI*WXXX+XYI*WXXY+0.5*YYI*WXYY

      FAY=SS*WY+0.5*(XXI*WXXY+2.0*XYI*WXYY+YYI*WYYY)

      FAZ=SS*WZ+0.5*(XXI*WXXZ+2.0*XYI*WXYZ+YYI*WYYZ)

      GO TO 160

C-3)-Process for (RL < 2.45 )

   40 IF(ABS(XD).LT.1.E-5.AND.ABS(YD).LT.1.E-5.AND.ABS(ZD).LT.1.E-5)

     *  GO TO 60

C-*-Cal. Potential

      FAI =0.0

        N=4

c      IF(RL.GT.1.0) GO TO 50

c        N=6

c      IF(RL.GT.0.5) GO TO 50

c       N=8

c      IF(RL.GT.0.1) GO TO 50

c        N=10

c   50 CALL HESS1(CXI,CETA,N,XD,YD,ZD,EL,FAI,ITI,ITO,KEY(34))
      CALL HESS1(CXI,CETA,N,XD,YD,ZD,EL,FAI)

C-*-Cal. derivative

      FAX=0.0

      FAY=0.0

      FAZ=0.0

c      CALL HESS2(CXI,CETA,XD,YD,ZD,EL,FAX,FAY,FAZ,ITI,ITO,KEY(34))
	CALL HESS2(CXI,CETA,XD,YD,ZD,EL,FAX,FAY,FAZ)

      GO TO 160

C-4)-Process for ( RL is nearly equal zero )

   60 P4 = PAI/4.0

      FAI = 0.0

      IF (ABS(CXI(1)).LT.1.E-6) GO TO 70

      CTAI = ATAN2(CETA(1),CXI(1))

      GO TO 90

   70 IF (CETA(1).LT.0.0) GO TO 80

      CTAI = PAI / 2.0

      GO TO 90

   80 CTAI  = -PAI / 2.0

   90 DO 150 I=1,4

        Y1 = CXI(I+1) - CXI(I)

        Y2 = CETA(I+1) - CETA(I)

        Y3 = CXI(I)*CETA(I+1) - CETA(I)*CXI(I+1)

        ELI = SQRT(Y1*Y1 + Y2*Y2)

        P = ABS(Y3) / ELI

        IF (ABS(CXI(I+1)).LT.1.E-6)  GO TO 100

        CTAI1 = ATAN2 (CETA(I+1),CXI(I+1))

        GO TO 120

  100   IF(CETA(I+1).LT.0.0) GO  TO 110

          CTAI1 = PAI / 2.0

          GO TO 120

  110   CTAI1 = -PAI / 2.0

  120   CTAI0 = ATAN2(-Y1,Y2)

C---Cal check No.33

        IF(KEY(32).LE.0) GO TO 130

          WRITE(ITO,2000) CTAI0

  130   A1 = (CTAI1 - CTAI0) /2.0 + P4

        A2 = (CTAI  - CTAI0) /2.0 + P4

        A1=SIN(A1) / COS(A1)

        A2=SIN(A2) / COS(A2)

        A3 = ABS(A1/A2)

        FAI = FAI + P*ALOG(A3)

C-*-Cal. check No.33

        IF(KEY(32).LE.0) GO TO 140

          WRITE(ITO,4000) Y1,Y2,Y3,ELI,P,CTAI,CTAI1,CTAI0,A1,A2,A3,FAI

  140 CTAI =CTAI1

  150 CONTINUE

      FAX = 0.0

      FAY = 0.0

      FAZ = 0.0

C---Cal. check No.32

c  160 IF(KEY(33).LE.0) GO TO 900
160	GO TO 900

      WRITE(ITO,3000) FAI,FAX,FAY,FAZ

C---Output format

 1000 FORMAT(1H ,5X,3HXD=,E12.5,3X,3HYD=,E12.5,3X,3HZD=,E12.5)
 1001 format(2x,3(2x,f12.6))
 2000 FORMAT(1H ,5X,10H*CTAI0* = ,E12.5)

 3000 FORMAT(1H0,10X,4HFAI=,F12.8,3X,4HFAX=,E12.5,3X,4HFAY=,E12.5,3X

     *      ,4HFAZ=,E12.5    )

 4000 FORMAT(1H ,5X,3HY1=,E12.5,3X,3HY2=,E12.5,3X,3HY3=,E12.5,3X,4HELI=

     *      ,E12.5,3X,2HP=,E12.5 / 6X

     *      ,5HCTAI=,E12.5,3X,6HCTAI1=,E12.5,3X,6HCTAI0=,E12.5,3X

     *      ,3HA1=,E12.5,3X,3HA2=,E12.5,3X,3HA3=,E12.5 / 6X,4HFAI=

     *      ,E12.5   )

C---Normal Return

  900 RETURN

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

      COMMON/CONST / ICONST,PI

      COMMON/INTGRL/ GAMMA,POLE1,POLE2

      COMMON/AUTPMT/ SRBASE(4,32),SIBASE(4,32),BOTTOM,PERIOD,EPS,IPHASE

     *              ,NSECT

      COMMON/COUNT/JDIV

      COMMON/SELECT/ISEL1

      COMMON/RESULT/ GRENR(4),GRENI(4)

      DATA PI2 /1.570796326794896619D0/

C---Initial Value Set

C-*-Definition of Values ------------------------------

C     DOMEGA    - Frequency of Oscillation (rad/s)

C     SPEED    - Steady Forward Speed of Source (m/s)

C------------------------------------------------------

      PI    = 3.141592653589793238D0

      ICONST= (0.D0,1.D0)

C---Information of Cal.Condition

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

C

C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

C   Name    -- Subroutine AUTDIV

C   Purpose -- Divide into subintervals by automatic processing.

C              Quadrate Integrand by Gauss-Legendre Method.

C   History -- Developed by Y.Makino    1988.10.1

C   Notes   -- Double Precision

C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      SUBROUTINE AUTDIV( ISECT,SECSR,SECSI )

      IMPLICIT REAL*8(A-H,O-Z)

      COMMON/AUTPMT/ SRBASE(4,32),SIBASE(4,32),BOTTOM,PERIOD,EPS,IPHASE

     *              ,NSECT

      COMMON/SELECT/ ISEL1

      DIMENSION  SECSR(4),SECSI(4)

     *          ,DIVSR(4,32),DIVSI(4,32),SRERR(4),SIERR(4)

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

   90 CONTINUE

      GO TO 10

C---Result of Integral don't satisfy the prescribed precision.

C    Preparing Process for Phase2

  100 IPHASE=2

      DO 110 IDIV=1,NSECT

        ISECT=IDIV

        DO 110 IDERIV=1,4

          SRBASE(IDERIV,ISECT)=DIVSR(IDERIV,IDIV)

          SIBASE(IDERIV,ISECT)=DIVSI(IDERIV,IDIV)

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

C

C

C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

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

        VAL=U(J)

        CALL GREENF (VAL,GRID)

        DO 3 IDERIV=1,4

          GRIDC=GRID(IDERIV)

          GRIDR=DREAL(GRIDC)

          GRIDI=DIMAG(GRIDC)

          SUMR(IDERIV,IDIV)=SUMR(IDERIV,IDIV)+AA(J)*GRIDR

          SUMI(IDERIV,IDIV)=SUMI(IDERIV,IDIV)+AA(J)*GRIDI

    3   CONTINUE

      DO 4 IDERIV=1,4

        SUMR(IDERIV,IDIV)=SUMR(IDERIV,IDIV)*C2

        SUMI(IDERIV,IDIV)=SUMI(IDERIV,IDIV)*C2

    4 CONTINUE

C---Normal Return

      JDIV=JDIV+1

      RETURN

      END

C

C

C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

C   Name    -- Subroutine GREENF

C   Purpose -- Cal. Integrand of Theta range Quadrature

C   Notes   -- Double Precision

C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      SUBROUTINE GREENF (THETA,GRID)

      IMPLICIT REAL*8 (A-H,O-Z)

      COMPLEX*16 GRID,GRIDPV,GRIDDT,QWBAS1,QWBAS2,QWNUM1,QWNUM2,QWNUM3,

     *           QWDEN,HBAS1,HBAS2,HBAS3,K1,K2,FKJ,ICONST,PART3,

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

          GRIDDT(4)= 2.D0*ICONST      *(HBAS2(1)+HBAS2(2))/PART3

            GO TO 2

    1     GRIDDT(1)= 2.D0*ICONST      *(HBAS1(1)-HBAS1(2))/PART3

          GRIDDT(2)=-2.D0       *PART1*(HBAS2(1)-HBAS2(2))/PART3

          GRIDDT(3)=-2.D0*ICONST*PART2*(HBAS3(1)-HBAS3(2))/PART3

          GRIDDT(4)= 2.D0*ICONST      *(HBAS2(1)-HBAS2(2))/PART3

    2   CONTINUE

C-<3>-Integrand Green Function

  100   DO 3 IDERIV=1,4

          GRID(IDERIV)=GRIDPV(IDERIV)+GRIDDT(IDERIV)

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
c      Name    -- Function FKJ
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

      COMPLEX*16 PART3,PART4

      COMMON/PARAMT/ GRAV,DOMEGA,SPEED,DX,DY,DZ,GRBATA

C---Cal.Process

      PART1=DCOS(THETA)

      PART2=(DOMEGA**2)/GRAV

      PART3=1.D0-4.D0*GRBATA*PART1

      PART3=CDSQRT(PART3)

      PART4=1.D0+PART3

      PART5=2.D0*GRBATA*PART1

      GO TO (1,2),J

    1 FKJ=4.D0*PART2/(PART4**2)

      GO TO 3

    2 FKJ=PART2*((PART4/PART5)**2)

C---Normal Return

    3 CONTINUE

      RETURN

      END

C

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

   60   DO 70 IDERIV=1,4

          GRENR(IDERIV)=GRENR(IDERIV)+RESR(IDERIV)

          GRENI(IDERIV)=GRENI(IDERIV)+RESI(IDERIV)

   70   CONTINUE

      RETURN

      END

C

C
	SUBROUTINE HESS1(X,Y,N1,XL,YL,ZL,EL,ANS)

      REAL INT
	
C

C

C

      DIMENSION  X(4),Y(4),CX(52),CW(52)

C

C

C*****    ABSCISSAE COEFFICIENTS

C

      DATA CX /-.7745966  ,0.,.7745966  ,-.8611363  ,-.3399810  ,

     A         .3399810  ,.8611363  ,-.9061798  ,-.5384693  ,0.,

     B         .5384963  ,.9061798  ,-.9324695  ,-.6612093  ,

     C         -.2386191  ,.2386191  ,.6612093  ,.9324695  ,

     D        -.9491079  ,-.7415311  ,-.4058451  ,0.,.4058451  ,

     E         .7415311  ,.9491079  ,-.9602898  ,-.7966664  ,

     F        -.5255324  ,-.1834346  ,.1834346  ,.5255324  ,

     G         .7966664  ,.9602898  ,-.9681602  ,-.8360311  ,

     H         -.6133714  ,-.3242534  ,0.,.3242534    ,.6133714  ,

     I         .8360311  ,.9681602  ,-.9739065  ,-.8650633  ,

     J        -.6794095  ,-.4333953  ,-.1488743  ,.1488743  ,

     K         .4333953  ,.6794095  ,.8650633  ,.9739065  /

C

C*****  WEIGHT COEFFICIENTS

C

      DATA CW / .5555556,.8888889,.5555556,.3478548,

     A         .6521452,.6521452,.3478543,.2369269,

     B         .4786286,.5688889,.4786287,.2369269,

     C         .1713244,.3607618,.4679139,.4679139,

     D         .3607615,.1713244,.1294850,.2797054,

     E         .3818300,.4179591,.3818300,.2797054,

     F         .1294850,.1012285,.2223810,.3137066,

     G		 .3626837,.3626837,.3137066,.2223810,	

     H         .1012285,.0812743,.1806482,.2606107,

     I         .3123471,.3302393,.3123471,.2606107,

     J         .1806481,.0812743,.0666713,.1494513,

     K         .2190863,.2692667,.2955242,.2955242,

     L         .2692667,.2190863,.1494513,.0666713/

C

C

      NN =  N1*(N1-1) /2-3

      ANS=0.

C

      DO 5000 I=1,4

      J=MOD(I,4)+1

      NN=N1*(N1-1)/2-3

      INT=0.

      IF(ABS(X(J)-X(I)).LT..1E-4*EL) GO TO 5000

      IF(ABS(ZL).GT..1E-4*EL) GO TO 1000

      IF((XL-X(I))*(XL-X(J)).GE..1E-4*EL) GO TO 1000

      IF(AMIN1(Y(I),Y(J))-YL.LT..1E-4*EL) GO TO 3000

 1000 DO 2000 N=1,N1

      NN=NN+1

      UMU=CX(NN)

      XI=(X(J)-X(I))*UMU/2.+(X(I)+X(J))/2.

      ET=(Y(J)-Y(I))*UMU/2.+(Y(I)+Y(J))/2.

      RR=SQRT((XL-XI)**2+(YL-ET)**2+ZL*ZL)

      ARG=ET-YL+RR

      INT=INT-ALOG(ARG)*CW(NN)

 2000 CONTINUE

      ANS=ANS+INT*(X(J)-X(I))/2.

      GO TO 5000

 3000 DO 4000 N=1,N1

      NN=NN+1

      UMU=CX(NN)

      XI=(X(J)-X(I))*UMU/2.+(X(I)+X(J))/2.

      ET=(Y(J)-Y(I))*UMU/2.+(Y(J)+Y(I))/2.

      RR=SQRT((XL-XI)**2+(YL-ET)**2+ZL*ZL)

      ARG=YL-ET+RR

      INT=INT+ALOG(ARG)*CW(NN)

 4000 CONTINUE

      ANS=ANS+INT*(X(J)-X(I))/2.

      ANS=ANS+2.*(X(J)-X(I))

      IF(ABS(X(J)-XL).GT..1E-4*EL) ANS=ANS-(X(J)-XL)*2.*

     *                             ALOG(ABS(X(J)-XL))

      IF(ABS(X(I)-XL).GT..1E-4*EL) ANS=ANS+(X(I)-XL)*2.*

     $                           ALOG(ABS(X(I)-XL))

C

 5000 CONTINUE

c      IF(ICK.LE.0) RETURN

c      WRITE(*,10000) ANS

10000 FORMAT(1H ,'ANS IN HESS1---',E12.5)

      RETURN

      END

C

C

      SUBROUTINE HESS2(X,Y,XL,YL,ZL,EL,VX,VY,VZ)

C

       DIMENSION X(4),Y(4)

      REAL MIJ

      VX=0.

      VY=0.

      VZ=0.

C

      DO 1000 I=1,4

      J=MOD(I,4)+1

      RI=SQRT((XL-X(I))**2+(YL-Y(I))**2+ZL*ZL)

      RJ=SQRT((XL-X(J))**2+(YL-Y(J))**2+ZL*ZL)

      DIJ=SQRT((X(J)-X(I))**2+(Y(J)-Y(I))**2)

      EI=ZL*ZL+(XL-X(I))**2

      HI=(YL-Y(I))*(XL-X(I))

      EJ=ZL*ZL+(XL-X(J))**2

      HJ=(YL-Y(J))*(XL-X(J))

      AL=(RI+RJ-DIJ)/(RI+RJ+DIJ)

      AL=ALOG(AL)

      IF(ABS(Y(J)-Y(I)).LT..1E-4*EL) GO TO 100

      VX=VX+(Y(J)-Y(I))*AL/DIJ

  100 IF(ABS(X(J)-X(I)).LT..1E-4*EL) GO TO 200

      VY=VY-(X(J)-X(I))*AL/DIJ

  200 IF(ABS(ZL).LT..1E-4*EL) GO TO 1000

      IF(ABS(X(J)-X(I)).LT..1E-4*EL) GO TO 1000

      MIJ=(Y(J)-Y(I))/(X(J)-X(I))

      TI=(MIJ*EI-HI)/(ZL*RI)

      TJ=(MIJ*EJ-HJ)/(ZL*RJ)

      TI=ATAN(TI)

      TJ=ATAN(TJ)

      VZ=VZ+TI-TJ

 1000 CONTINUE

C      IF(ICK.LE.0) RETURN

C      WRITE(ITO,2000) VX,VY,VZ

 2000 FORMAT(1H ,'VX,VY,VZ IN HESS2---',3E12.5)

      RETURN

      END

C

C
		SUBROUTINE MATINV(A,N)


	 COMPLEX A(480,480),
C	COMPLEX A,AMAX,SWAP,T,PIVOT
	
	
     &        AMAX,SWAP,T,
     &        PIVOT(480)
	integer Nelem,index1(480,2),
     &        j,i,k,irow,icolum,l,l1,jrow,jcolum,
     &        ipivot(480)
C      DIMENSION A(N,N),INDEX1(N,2),PIVOT(N),IPIVOT(N)                   MAT00150
C
C	WRITE(1,*)((A(I,J),I=1,N),J=1,N)
C      num=nelem
	NUM=N
	NELEM=N
C		WRITE(1,*)((A(I,J),I=1,NELEM),J=1,NELEM)
	DO 1000 J=1,num
      IPIVOT(J)=0                                                       MAT00250
 1000 CONTINUE                                                          MAT00300
      DO 2200 I=1,Nelem                                                 MAT00350
      AMAX=cmplx(0.0,0.0)                                               MAT00400
      DO 1500 J=1,Nelem                                                 MAT00450
      IF(IPIVOT(J)-1) 1100,1500,1100                                    MAT00500
 1100 DO 1400 K=1,Nelem                                                 MAT00550
      IF(IPIVOT(K)-1) 1200,1400,1400
 1200 IF(CABS(AMAX)-CABS(A(J,K))) 1300,1400,1400                        MAT00600
 1300 IROW=J                                                            MAT00650
      ICOLUM=K                                                          MAT00700
      AMAX=A(J,K)                                                       MAT00750
 1400 CONTINUE                                                          MAT00800
 1500 CONTINUE                                                          MAT00850
      IPIVOT(ICOLUM)=IPIVOT(ICOLUM)+1                                   MAT00900
      IF(IROW-ICOLUM) 1600,1800,1600                                    MAT00950
 1600 DO 1700 L=1,Nelem                                                 MAT01000
      SWAP=A(IROW,L)                                                    MAT01050
      A(IROW,L)=A(ICOLUM,L)                                             MAT01100
      A(ICOLUM,L)=SWAP                                                  MAT01150
 1700 CONTINUE                                                          MAT01200
 1800 INDEX1(I,1)=IROW                                                  MAT01250
      INDEX1(I,2)=ICOLUM                                                MAT01300
      PIVOT(I)=A(ICOLUM,ICOLUM)                                         MAT01350
      A(ICOLUM,ICOLUM)=(1.0,0.0)                                        MAT01400
      DO 1900 L=1,Nelem                                                 MAT01450
      A(ICOLUM,L)=A(ICOLUM,L)/PIVOT(I)                                  MAT01500
 1900 CONTINUE                                                          MAT01550
      DO 2200 L1=1,Nelem                                                MAT01600
      IF(L1-ICOLUM) 2000,2200,2000                                      MAT01650
 2000 T=A(L1,ICOLUM)                                                    MAT01700
      A(L1,ICOLUM)=(0.0,0.0)                                            MAT01750
      DO 2100 L=1,Nelem                                                 MAT01800
      A(L1,L)=A(L1,L)-A(ICOLUM,L)*T                                     MAT01850
 2100 CONTINUE                                                          MAT01900
 2200 CONTINUE                                                          MAT01950
      DO 2500 I=1,Nelem                                                 MAT02000
      L=Nelem+1-I                                                       MAT02050
      IF(INDEX1(L,1)-INDEX1(L,2)) 2300,2500,2300                        MAT02100
 2300 JROW  =INDEX1(L,1)                                                MAT02150
      JCOLUM=INDEX1(L,2)                                                MAT02200
      DO 2400 K=1,Nelem                                                 MAT02250
      SWAP=A(K,JROW)                                                    MAT02300
      A(K,JROW)=A(K,JCOLUM)                                             MAT02350
      A(K,JCOLUM)=SWAP                                                  MAT02400
 2400 CONTINUE                                                          MAT02450
 2500 CONTINUE                                                          MAT02500
      RETURN                                                            MAT02550
      END                                                               MAT02600 