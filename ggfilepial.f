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
c	THIS IS THE FIRST PROGRAM FOR THE CALCULATION  OF STEADY FREE SURFACE FLOW PATTERN     *
C																						   *	
C					THIS PROGRAM WILL GENERATE INPUT FOR OTHERS PROGRAMS				   *	
C																						   *	
c*********************************************************************************************
c

C---Input File to Green Function generator

      PARAMETER(PAR=1000)
	
	COMPLEX GAMMA,ABC,F,dndphi
	

      REAL N1,N2,N3,N4,N5,N6,IXX,IXY,IYY,LMAX

	INTEGER INF


      COMMON/STRG/RDATA(100)/CENTR/XYZG(3)

      COMMON/WAVE/NTW,TWAVE(5),IWW,FN,FORWRD,KLK

      COMMON /MUGEN/INF,rk,v

      COMMON/NONDIM/ALG,ABD,DVV

      COMMON/GRV/ DGX(480),DGY(480),DGZ(480),DSS(480)

      COMMON/PANEL/XVALUE(4,480),YVALUE(4,480),ZVALUE(4,480)

	DIMENSION XGG(PAR),YGG(PAR),ZGG(PAR),SS(PAR)

      DIMENSION NT(4,480),NV(PAR),NE(PAR),X(1000),Y(1000),

     1          Z(1000),CR(6,6),PG(3),W2(3,3)


      DIMENSION N1(PAR),N2(PAR),N3(PAR),N4(PAR),

     * N5(PAR),N6(PAR),dndphi(PAR,6),XLENG(PAR),YLENG(PAR),zleng(PAR)	



	DIMENSION FAIN1(PAR),FAIN2(PAR),FAIN3(PAR)

      DIMENSION FAI2(PAR),FAIN4(PAR),FAI1(PAR)

      DIMENSION FAIX1(PAR),FAIX2(PAR)

      DIMENSION XE(4),YE(4),ZE(4),WW(3),XC(4),YC(4),ROT(9)

      DIMENSION GAMMA(PAR),SA(2),xxg(3)

      DIMENSION SG(3),XEE(4),YEE(4),ZEE(4)

C	INPUT DATA FOR BODY 
	OPEN(1,FILE='data/data_5',STATUS='UNKNOWN')


c	OUTPUT DATA
	OPEN(2,FILE='output/ggfile',STATUS='UNKNOWN')
	OPEN(3,FILE='output/centroid',STATUS='UNKNOWN')
	OPEN(4,FILE='output/enxyz',STATUS='UNKNOWN')
	OPEN(5,FILE='output/area',STATUS='UNKNOWN')
	open(10,FILE='output/dndphi',STATUS='UNKNOWN')
	OPEN(11,FILE='output/fain',STATUS='UNKNOWN')
c	OPEN(12,FILE='svgren',STATUS='UNKNOWN')	
	OPEN(13,FILE='output/Length',STATUS='UNKNOWN')

	WRITE(*,*)'WHAT IS THE ANGLE i,e, BETA ?'
	READ(*,*)ANGL


	AMP=0.00
	GRAV=9.8
	omega=.001
	fn=.2
c	omega=2.155
c	fn=.50
	RK=omega*omega/grav

cc	inf=1 for infinite depth
cc	inf=0 for finite depth in that case 
cc	 Depth of water H is nedded

	INF=1
	H=0.0
	RS=1.0
	NDIV=NELEM
	xxg(1)=0.0
	xxg(2)=0.0
	xxg(3)=0.0
	rc=1.0
c	rho=1000 kg/m3,so
	rho=1000.0

	pi=acos(-1.0)

	BATA=0.0174532925199*ANGL
	BETA=0.0174532925199*ANGL

C
C----->0.01745329251994329576e0=PI/180
C
		 
      CBATA=RK*COS(beta)
      SBATA=RK*SIN(beta)

	B=GRAV*AMP
      C=OMEGA
      CONST=-B/C

	READ(1,*)NELEM
	WRITE(*,*)nelem

	DO 101 J=1,NELEM
	read(1,1110) NE(J),(NT(IP,J),IP=1,4)
101	continue

c	READ(1,110)NNODE

	READ(1,*)NNODE
	write(*,*) nnode

	 DO 102 J=1,NNODE
 
      read(1,1120) NV(J),X(J),Y(J),Z(J)

102	continue	  


  110	FORMAT(8X,I4)
 1109 FORMAT(A4,4X,I4,4X,4(I4,1X))
 1119 FORMAT(A4,4X,I4,3E12.5)



 1110 FORMAT(8X,I4,4X,4(I4,1X))
 1120 FORMAT(8X,I4,4x,3(E12.5,2x))

	
	do 700 j=1,6
	DO 700 I=1,NELEM
700	dndphi(i,j)=0.0

	do 800 i=1,nelem
	fai1(i)=0.0
	fai2(i)=0.0
	fain1(i)=0.0
	fain2(i)=0.0
	faix1(i)=0.0
	faix2(i)=0.0
  
  800 CONTINUE



C---Parameter Initialize

      DV=0.0

      ST=0.

      XB=0.

      YB=0.

      ZB=0.

      WA=0.

      BA=0.

      FX=0.

      FY=0.

      BX=0.

      BY=0.

      WXX=0.

      WYY=0.

      WXY=0.

      XMAX=-100000.

      YMAX=-100000.

	ZMAX=-100000.

      XMIN=100000.

      YMIN=100000.

	ZMIN=100000.

      DO 950 I=1,6

      DO 950 J=1,6

      CR(J,I)=0.0
	
  950 CONTINUE

C---Repeat Element Number Times

      DO 1600 I=1,NELEM

ct      ISP=0

      R0=0.

      DO 1200 J=1,4

      NODE=NT(J,I)

      DO 1000 L=1,NNODE

      IF(NV(L)-NODE) 1000,1100,1000

c	1000 f0r -ve, 1100 for 0, 1000 for +ve

 1000 CONTINUE

      IER=IER+1

      GO TO 1200

1100  XE(J)=X(L)

      YE(J)=Y(L)

      ZE(J)=Z(L)


      XMAX=AMAX1(XMAX,XE(J))
      YMAX=AMAX1(YMAX,YE(J))
	ZMAX=AMAX1(ZMAX,ZE(J))

      XMIN=AMIN1(XMIN,XE(J))
      YMIN=AMIN1(YMIN,YE(J))
	ZMIN=AMIN1(ZMIN,ZE(J))

 1200 CONTINUE


	CALL COMOD(XE,YE,ZE)



 1299 DO 1300 J=1,2

      K=J+1

      L=J+2

      SA(J)=((YE(K)-YE(1))*(ZE(L)-ZE(1))-(ZE(K)-ZE(1))*(YE(L)-YE(1)))**2

     1     +((ZE(K)-ZE(1))*(XE(L)-XE(1))-(XE(K)-XE(1))*(ZE(L)-ZE(1)))**2

     2     +((XE(K)-XE(1))*(YE(L)-YE(1))-(YE(K)-YE(1))*(XE(L)-XE(1)))**2

      SA(J)=SQRT(SA(J))

 1300 CONTINUE

C

	asa=abs((xe(2)-xe(1))**2+(ye(2)-ye(1))**2)
	asa=sqrt(asa)
	yleng(i)=asa

	write(13,*)YLENG(I)		



      XG1=(XE(1)+XE(2)+XE(3))/3.

      XG2=(XE(1)+XE(3)+XE(4))/3.

      YG1=(YE(1)+YE(2)+YE(3))/3.

      YG2=(YE(1)+YE(3)+YE(4))/3.

      ZG1=(ZE(1)+ZE(2)+ZE(3))/3.

      ZG2=(ZE(1)+ZE(3)+ZE(4))/3.

C

      XG=(SA(1)*XG1+SA(2)*XG2)/(SA(1)+SA(2))

      YG=(SA(1)*YG1+SA(2)*YG2)/(SA(1)+SA(2))

      ZG=(SA(1)*ZG1+SA(2)*ZG2)/(SA(1)+SA(2))

c	

	write(3,*)xg,yg,zg
c
	sg(1)=xg
	sg(2)=yg
	sg(3)=zg


      XGG(I)=XG

      YGG(I)=YG

      ZGG(I)=ZG

      DGX(I)=XG

      DGY(I)=YG

      DGZ(I)=ZG

      CALL TRAN(ROT,XE,YE,ZE)

      DO 1500 J=1,4

      WW(1)=XE(J)-XG

      WW(2)=YE(J)-YG

      WW(3)=ZE(J)-ZG

      SUM=0.0

      TUM=0.0

      DO 1400 L=1,3

	SUM=SUM+ROT(L)*WW(L)

      TUM=TUM+ROT(L+3)*WW(L)

 1400 CONTINUE

 
c   xc --> partial coordinate of X-axis
c   yc --> partial coordinate of Y-axis
c

      XC(J)=SUM

      YC(J)=TUM

 1500 CONTINUE

c1530 CONTINUE

      CA=XC(3)-XC(1)

      CB=YC(4)-YC(2)

      CC=YC(1)-YC(2)

      CD=YC(1)+YC(2)

      CE=YC(1)-YC(4)

      CF=YC(1)+YC(4)

      CG=XC(3)+XC(1)

      S=0.5*CA*CB

c     
	write(5,*) s
c

      DSS(I)=S

      IXX=S*(CA*CA+3.0*XC(1)*XC(3)-XC(2)*XC(4))/6.0

      IXY=CA*(2.0*XC(2)*CC*CD-2.0*XC(4)*CE*CF+CG*(CF*CF-CD*CD))/24.0

      IYY=S*(YC(1)*CC-YC(4)*CD)/6.0

      SS(I)=S

ct      IF(KLK.EQ.2.AND.KL.EQ.1) GO TO 1510

      PS=S*(ROT(7)*XG+ROT(8)*YG+ROT(9)*ZG)

      DV=DV+PS

      XB=XB+XG*PS

      YB=YB+YG*PS

      ZB=ZB+ZG*PS

      ST=ST+S

      IF(ABS(ZG).LE..1E-4) GO TO 1510

      DO 1501 J=1,4

      J1=MOD(J,4)+1

      IF(ABS(ZE(J)+ZE(J1)).LT..1E-4) GO TO 1502

      IF(ABS(ZE(J)+ZE(J1)+2.*H).LT..1E-4) GO TO 1503

 1501 CONTINUE

      GO TO 1510
 1502 DS=XE(J1)*YE(J)-XE(J)*YE(J1)

      WA=WA+DS

      FX=FX+DS*(XE(J)+XE(J1))

      FY=FY+DS*(YE(J)+YE(J1))

      WXX=WXX+DS*(YE(J)*YE(J)+YE(J1)*YE(J1)+YE(J)*YE(J1))

      WYY=WYY+DS*(XE(J)*XE(J)+XE(J1)*XE(J1)+XE(J)*XE(J1))

      WXY=WXY+DS*((XE(J)+XE(J1))*(YE(J)+YE(J1))+XE(J)*YE(J)+

     *      XE(J1)*YE(J1))

ct      ISP=1

      R0=((XE(J)-XE(J1))**2+(YE(J)-YE(J1))**2)/2.

      GO TO 1510

 1503 DS=XE(J)*YE(J1)-XE(J1)*YE(J)

      BA=BA+DS

      BX=BX+DS*(XE(J)+XE(J1))

      BY=BY+DS*(YE(J)+YE(J1))

ct	e.q. (56) Vector r_13

 1510 XX=XE(3)-XE(1)

      YY=YE(3)-YE(1)

      ZZ=ZE(3)-ZE(1)

      SUM=XX**2+YY**2

      SUM1=SUM+ZZ*ZZ

c	e.q. (56) Vector r_24

      XX=XE(4)-XE(2)

      YY=YE(4)-YE(2)

      ZZ=ZE(4)-ZE(2)

      TUM=XX**2+YY**2

      TUM1=TUM+ZZ*ZZ

      LMAX=SQRT(AMAX1(SUM1,TUM1))

      HMAX=SQRT(AMAX1(SUM,TUM))


	WRITE(2,*) NE(I),(XC(J),J=1,4),(YC(J),J=1,4),(ROT(J),J=1,9),S,
     1				IXX,IXY,IYY,LMAX,HMAX



	CALL GAUSS1(XE,YE,ZE,CBATA,SBATA,GAMMA(I))


      TT=EXP(RK*ZG)



       A=XG*CBATA

       B=YG*SBATA

      IF(INF.GE.1) GO TO 2900



 2900 CONTINUE



	FAI1(I) =-CONST*TT*SIN(A+B)

      FAI2(I) =+CONST*TT*COS(A+B)

      FAIX1(I)=-CONST*TT*COS(A+B)*CBATA

      FAIX2(I)=-CONST*TT*SIN(A+B)*CBATA

	write(11,*)fai1(i),fai2(i)


C

	C1=COS(A)

      C2=COS(B)

      S1=SIN(A)

      S2=SIN(B)

      IF(INF.GE.1) GO TO 3110


 3110 T1=CONST*CBATA*TT*ROT(7)

      T2=CONST*SBATA*TT*ROT(8)

      T3=CONST*RK*TT*ROT(9)


c*****************************************
c   Fain1 - Normal Differential of Fai1
c   Fain2 - Normal Differential of Fai2
******************************************

C 3300 FAIN1(I)=(T1+T2)*COS(A+B)+T3*SIN(A+B)+FAIN1(I)
	FAIN1(I)=(T1+T2)*COS(A+B)+T3*SIN(A+B)+FAIN1(I)

      FAIN2(I)=(T1+T2)*SIN(A+B)-T3*COS(A+B)+FAIN2(I)

c	write(11,*)fain1(i),fain2(i)

C 3400 CONTINUE


C

C 1280 CONTINUE

 1580 N1(I)=ROT(7)

      N2(I)=ROT(8)

      N3(I)=ROT(9)

      N4(I)=(YG-XYZG(2))*ROT(9)-(ZG-XYZG(3))*ROT(8)

      N5(I)=(ZG-XYZG(3))*ROT(7)-(XG-XYZG(1))*ROT(9)

      N6(I)=(XG-XYZG(1))*ROT(8)-(YG-XYZG(2))*ROT(7)

      DO 1585 J=1,3

      CR(J,3)=CR(J,3)+ROT(J+6)*S

	 J1=MOD(J,3)+1

      J2=MOD(J+1,3)+1

      RN=(SG(J1)-XXG(J1))*ROT(J2+6)-(SG(J2)-XXG(J2))*ROT(J1+6)

      CR(J+3,3)=CR(J+3,3)+RN*S

      CR(J,4)=CR(J,4)+ROT(J+6)*(SG(2)-XXG(2))*S

      CR(J,5)=CR(J,5)-ROT(J+6)*(SG(1)-XXG(1))*S

      CR(J,J2+3)=CR(J,J2+3)-ROT(J1+6)*SG(3)*S

      CR(J,J1+3)=CR(J,J1+3)+ROT(J2+6)*SG(3)*S

C************

      CR(J+3,4)=CR(J+3,4)+RN*(SG(2)-XXG(2))*S

      CR(J+3,5)=CR(J+3,5)-RN*(SG(1)-XXG(1))*S

      RN1=(SG(J2)-XXG(J2))*ROT(J+6)-(SG(J)-XXG(J))*ROT(J2+6)

      CR(J2+3,J1+3)=CR(J2+3,J1+3)-RN*SG(3)*S

      CR(J2+3,J+3)=CR(J2+3,J +3)+RN1*SG(3)*S

 1585 CONTINUE

C     IF(KEY(22)-1)1600,1590,1590

C 1590 WRITE(ITO,1049) N1(I),N2(I),N3(I),N4(I),N5(I),N6(I)
c	WRITE(2,*) N1(I),N2(I),N3(I),N4(I),N5(I),N6(I)


 1600  CONTINUE

	BA =BA/2.0
      DV =(DV+BA*H)/3.0
      DVV=DV
      XB =(XB+BX*H/6.0)/DV/4.0
      YB =(YB+BY*H/6.0)/DV/4.0
      ZB =(ZB-BA*H**2.0)/DV/4.0
      WA =WA/2.0
      FX =FX/WA/6.0
      FY =FY/WA/6.0
      WXX=WXX/12.0
      WYY=WYY/12.0
      WXY=WXY/24.0
      WXX=WXX-FY*FY*WA
      WYY=WYY-FX*FX*WA
      WXY=WXY-FX*FY*WA
      ALG=XMAX-XMIN 
      ABD=YMAX-YMIN
C*********************
C   Real Forward Speed
C*********************
      FORWRD=FN*SQRT(GRAV*ALG)
C**********************************************
C   Encounter Circular Frequency of Oscillation
C**********************************************
c      OMEGA0=OMEGA
	OMEGAE =OMEGA-(OMEGA**2/GRAV)*FORWRD*COS(beta)
C*********
c      DUMY(1)=DV*GRAV*rho
c      DUMY(2)=WA
	GL=OMEGAe*OMEGAe/GRAV*ALG
	DO 1605 J=1,6
      DO 1605 K=1,6
 1605 CR(K,J)=-CR(K,J)*GRAV*rho      


 1701 FORMAT(1H0,5X,'********  GZAI-L(OMEGA-E**2/G*L) :',F10.5
     &          ,' ( WAVE NUMBER :',F9.5
     &          ,', OMEGA-O :'     ,F9.5
     &          ,', OMEGA-E :'     ,F9.5,' ) ************')

C**************************
C   Hull Boundary Condition
C**************************
C	write(*,*)'forward speed=',forwrd,'alg=',alg,fn,abd
C	write(*,*)'omegae=',omegae
      DO 2500 ii=1,nelem

C**********************************************
C   Real & Imag Value Set for Radiation Problem
C**********************************************
c        HBC1R = N1(IPANEL)
c	  HBC2R = N2(IPANEL)
c	  HBC3R = N3(IPANEL)
c	  HBC4R = N4(IPANEL)
c	  HBC5R = N5(IPANEL)
c	  HBC6R = N6(IPANEL)
c	  HBC1I = 0.0
c	  HBC2I = 0.0
c	  HBC3I = 0.0
c	  HBC4I = 0.0
c	  HBC5I = +FORWRD*N3(IPANEL)/OMEGAe
c	  HBC6I = -FORWRD*N2(IPANEL)/OMEGAe
       dndphi(ii,1)=CMPLX(n1(ii),0.0)
       dndphi(ii,2)=CMPLX(n2(ii),0.0)
	 dndphi(ii,3)=CMPLX(n3(ii),0.0)
	 dndphi(ii,4)=CMPLX(n4(ii),0.0)
	 dndphi(ii,5)=CMPLX(n5(ii), forwrd*n3(ii)/omegae)
	 dndphi(ii,6)=CMPLX(n6(ii),-forwrd*n2(ii)/omegae)

	write(4,1)dndphi(ii,1),dndphi(ii,2),dndphi(ii,3)
	write(10,*)(dndphi(ii,j),j=1,6)
	
	
	
2500	CONTINUE


    1 format(3(2f12.8))
    2 format(6(2f12.8))	
	STOP
	END

     
	SUBROUTINE TRAN(TOR,X,Y,Z)                                        

C*****  O DIMENSIONAL ROTATION MATRIX                                   

      DIMENSION ROT(9),X(4),Y(4),Z(4),TOR(9)

      DOUBLE PRECISION A,B,C,ROT

      ROT(1)=X(3)-X(1)                                                  TRA00150

      ROT(2)=Y(3)-Y(1)                                                  TRA00200

      ROT(3)=Z(3)-Z(1)                                                  TRA00250

      ROT(4)=X(4)-X(2)                                                  TRA00300

      ROT(5)=Y(4)-Y(2)                                                  TRA00350

      ROT(6)=Z(4)-Z(2)                                                  TRA00400

      ROT(7)=ROT(2)*ROT(6)-ROT(3)*ROT(5)                                TRA00450

      ROT(8)=ROT(3)*ROT(4)-ROT(1)*ROT(6)                                TRA00500

      ROT(9)=ROT(1)*ROT(5)-ROT(2)*ROT(4)                                TRA00550

      ROT(4)=ROT(3)*ROT(8)-ROT(2)*ROT(9)                                TRA00600

      ROT(5)=ROT(1)*ROT(9)-ROT(3)*ROT(7)                                TRA00650

      ROT(6)=ROT(2)*ROT(7)-ROT(1)*ROT(8)                                TRA00700

      A=ROT(1)*ROT(1)+ROT(2)*ROT(2)+ROT(3)*ROT(3)

      A=1.0/SQRT(A)                                                     TRA00800

      B=ROT(7)*ROT(7)+ROT(8)*ROT(8)+ROT(9)*ROT(9)

      B=1.0/SQRT(B)                                                     TRA00900

      C=A*B                                                             TRA00950

      DO 1000 I=1,3                                                     TRA01000

      TOR(I)=ROT(I)*A                                                   TRA01050

      TOR(I+3)=ROT(I+3)*C                                               TRA01100

      TOR(I+6)=ROT(I+6)*B                                               TRA01150

 1000 CONTINUE                                                          TRA01200

      RETURN                                                            TRA01250

      END                                                               TRA01300

C

C

	SUBROUTINE COMOD(X,Y,Z)                                           

C*****  REDEFINE CO-ORDINATES DATA IF NECESSARY                          CMD00050

      COMMON/FILE/KEY(72),ITI,ITO                                       CMD00100

      DIMENSION X(4),Y(4),Z(4)                                          CMD00150

c    IF(KEY(15).EQ.0) GO TO 1000                                       CMD00200

c      WRITE(*,*) X,Y,Z                                               

c1009 FORMAT(1H0,29H*** COMOD --- DEGUG *** X,Y,Z/(1X,4E12.4))          CMD00300

c 1000 X21=X(2)-X(1)                                                     CMD00350

      X21=X(2)-X(1) 

      Y21=Y(2)-Y(1)                                                     CMD00400

      Z21=Z(2)-Z(1)                                                     CMD00450

      X41=X(4)-X(1)                                                     CMD00500

      Y41=Y(4)-Y(1)                                                     CMD00550

      Z41=Z(4)-Z(1)                                                     CMD00600

      XA=Y21*Z41-Y41*Z21                                                CMD00650

      YA=Z21*X41-Z41*X21                                                CMD00700

      ZA=X21*Y41-X41*Y21                                                CMD00750

      X43=X(4)-X(3)                                                     CMD00800

      Y43=Y(4)-Y(3)                                                     CMD00850

      Z43=Z(4)-Z(3)                                                     CMD00900

      X23=X(2)-X(3)                                                     CMD00950

      Y23=Y(2)-Y(3)                                                     CND01000

      Z23=Z(2)-Z(3)                                                     CMD01050

      XB=Y43*Z23-Y23*Z43                                                CMD01100

      YB=Z43*X23-Z23*X43                                                CMD01150

      ZB=X43*Y23-X23*Y43                                                CMD01200

      XA=XA-X(1)                                                        CMD01250

      YA=YA-Y(1)                                                        CMD01300

      ZA=ZA-Z(1)                                                        CMD01350

      XB=XB-X(3)                                                        CMD01400

      YB=YB-Y(3)                                                        CMD01450

      ZB=ZB-Z(3)                                                        CMD01500

      CX=XA-XB                                                          CMD01550

      CY=YA-YB                                                          CMD01600

      CZ=ZA-ZB                                                          CMD01650

      CXY=CX-CY                                                         CMD01700

      CYZ=CY-CZ                                                         CMD01750

      IF(CXY.GT.1.0E-4)  GO  TO 1100

      IF(CYZ.GT.1.0E-4) GO TO 1100

      IF(KEY(15).EQ.0) GO TO 1400

      WRITE(2,1019) CXY,CYZ                                             CMD01950

 1019 FORMAT(1H0,35H*** COMOD --- DEBUG *** CX-CY,CY-CZ/1X,2E12.4)      CM:02000

      GO TO 1400                                                        CMD02050

 1100 X21=X(3)-X(1)                                                     CMD02100

      Y21=Y(3)-Y(1)                                                     CMD02150

      Z21=Z(3)-Z(1)                                                     CMD02200

      X41=X(4)-X(2)                                                     CMD02250

      Z41=Z(4)-Z(2)                                                     CMD02350

      Y41=Y(4)-Y(2)                                                     CMD02300

      XA=Y21*Z41-Y41*Z21                                                CMD02400

      YA=Z21*X41-Z41*X21                                                CMD02450

      ZA=X21*Y41-X41*Y21                                                CMD02500

      SUM=XA**2+YA**2+ZA**2                                             CMD02550

      SUM=SQRT(SUM)                                                     CMD02600

      XA=XA/SUM                                                         CMD02650

      YA=YA/SUM                                                         CMD02700

      ZA=ZA/SUM                                                         CMD02750

      XB=0.0                                                            CMD02800

      YB=0.0                                                            CMD02850

      ZB=0.0                                                            CMD02900

      DO 1200 I=1,4                                                     CMD02950

      XB=XB+X(I)                                                        CMD03000

      YB=YB+Y(I)                                                        CMD03050

      ZB=ZB+Z(I)                                                        CMD03100

 1200 CONTINUE                                                          CMD03150

      XB=0.25*XB                                                        CMD03200

      YB=0.25*YB                                                        CMD03250

      ZB=0.25*ZB                                                        CMD03300

      CC=XB*XA+YB*YA+ZB*ZA                                              CMD03350

      DO 1300 I=1,4                                                     CMD03400

      SUM=X(I)*XA+Y(I)*YA+Z(I)*ZA                                       CMD03450

      SUM=SUM-CC                                                        CMD03500

      X(I)=X(I)-SUM*XA                                                  CMD03550

      Y(I)=Y(I)-SUM*YA                                                  CMD03600

      Z(I)=Z(I)-SUM*ZA                                                  CMD03650

 1300 CONTINUE                                                          CMD03700

c      IF(KEY(15).EQ.0) GO TO 1400                                       CMD03750

c     WRITE(2,*) X, Y, Z                                                    

 1400 RETURN                                                            CMD03850

      END                                                               CMD03900

C


C
c	SUBROUTINE GAUSS1 (X,Y,Z,COBK,SIBK,INTEGR,COCO)
	SUBROUTINE GAUSS1 (X,Y,Z,COBK,SIBK,INTEGR)

      COMPLEX ABC,F,INTGR0,INTEGR

      DOUBLE PRECISION CX,CW

      COMMON /FILE/KEY(72),ITI,ITO

      COMMON/IWAVE/RK,V,H

      COMMON/MUGEN/INF

      DIMENSION  X(4), Y(4),Z(4), RA(4),UM(4),CX(52),CW(52)

       DATA RA/-1.0,1.0,1.0,-1.0/,UM/-1.0,-1.0,1.0,1.0 /

C*****  WEIGHT COEFFICIENTS

      DATA CX /-.7745966  ,0.,.7745966  ,-.8611363  ,-.3399810  ,

     A         .3399810,.8611363,-.9061798,-.5384693,0.,

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

      DATA CW / .5555556,.8888889,.5555556,.3478548,

     A         .6521452,.6521452,.3478543,.2369269,

     B         .4786286,.5688889,.4786287,.2369269,

     C         .1713244,.3607618,.4679139,.4679139,

     D         .3607615,.1713244,.1294850,.2797054,

     E         .3818300,.4179591,.3818300,.2797054,

     F         .1294850,.1012285,.2223810,.3137066,

     G         .3626837,.3626837,.3137066,.2223810,

     H         .1012285,.0812743,.1806482,.2606107,

     I         .3123471,.3302393,.3123471,.2606107,

     J         .1806481,.0812743,.0666713,.1494513,

     K         .2190863,.2692667,.2955242,.2955242,

     L         .2692667,.2190863,.1494513,.0666713/

c	write(*,*)'entered in gauss and assumed coco=1.0'
	

       N1=4

       N2=4

      NN =  N1*(N1-1) /2 - 3

      MM = N2*(N2-1) / 2 - 3

      INTEGR = (0.0,0.0)

      DO 4000 I=1,N1

      NN = NN + 1

      UMU = CX( NN)

      INTGR0 = (0.0,0.0)

      DO 3000 J=1,N2

      MM = MM + 1

      RAM = CX( MM)

      COI = 0.0

      COJ = 0.0

      COK = 0.0

      DO 1000 K=1,3

      K1 = K+1

       DO 1000  L=K1,4

       CCC = RA(K)*UM(L)-RA(L)*UM(K) + RA(K)*RA(L)*(UM(L)-UM(K))*RAM

     1        +UM(K)*UM(L)*(RA(K)-RA(L))*UMU

      COI=COI+CCC*(Y(K)*Z(L)-Z(K)*Y(L))

      COJ=COJ+CCC*(Z(K)*X(L)-X(K)*Z(L))

        COK = COK + CCC*(  X(K)* Y(L)- Y(K)* X(L))

 1000  CONTINUE

      E3=SQRT(COI*COI+COJ*COJ+COK*COK)/16.0

      GZAI = 0.0

      ETA = 0.0

      ZETA = 0.0

       DO 2000 M=1,4

       COO = 0.25*(1.+RA(M)*RAM)*(1.+UM(M)*UMU)

       GZAI = GZAI + COO* X(M)

       ETA  = ETA  + COO* Y(M)

      ZETA=ZETA+COO*Z(M)

 2000  CONTINUE

      GE = GZAI*COBK + ETA*SIBK

      ABC = CMPLX(0.0,GE)

czaman      F  = CEXP(ABC)*COSH(RK*(ZETA+H))/COCO

c      IF(INF.GE.1) then
	 F=CEXP(ABC)*EXP(RK*ZETA)
c	else
c	F  = CEXP(ABC)*COSH(RK*(ZETA+H))/COCO
c	end if

      INTGR0=INTGR0+F*E3*CW(MM)

 3000 CONTINUE

      MM = MM-N2

      INTEGR = INTEGR + INTGR0*CW(NN)

 4000 CONTINUE

 5009 FORMAT(1H0,10X,10HINTEGR =  ,2(E12.5,1X))

 9000 RETURN

      END
