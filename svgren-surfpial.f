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
c	THIS IS THE THIRD PROGRAM FOR THE CALCULATION  OF STEADY FREE SURFACE FLOW PATTERN     *
C																						   *	
C	     THIS PROGRAM GENERATES SOURCE DENSITIES BY USING THE INPUT FROM MGGREN.F		   *	
C																						   *	
c********************************************************************************************



	integer nelem,ind(480,2),ipivot(480)
	real fain1(480),fain2(480),x(1000),y(1000),z(1000)
c    &     dndphi(480,6)
      complex alfa(480,480),dndphi(480,6),beta(480,480),
     &        pivot(480),phi(480,7),
     &        sigma(480,7)
c	
	open(3,file='output/dndphi',status='unknown')
	open(4,file='output/alfa_sig',status='unknown')
	open(7,file='output/sigma_body',status='unknown')






	write(*,*)'nelem= ?'
	read(*,*) nelem

	num=nelem
	ndiv=nelem

14	format(3(2f20.12,x))	
	

	read(4,*)((alfa(i,j),i=1,nelem),j=1,nelem)
	

	do 12 i=1,nelem	

	read(3,*)(dndphi(i,jj),jj=1,6)
12	continue


113	format(4(2f18.10))

	pi=acos(-1.0)
	pi4=-1.0/(4.0*pi)
	grav=9.80000065

      call matinv(alfa,nelem,ind,ipivot,pivot)
      
	
C	SIGMA= SOURCE DENSITIES AROUND THE BODY SURFACE
	
	do 100 jv=1,num                                                               
      do 100 iv=1,7                                                           
      sigma(jv,iv)=cmplx(0.0,0.0)                                                          
  100 continue                                                                  
c
      do 200 iv=1,nelem
      do 200 jv=1,nelem
	do 300 kv=1,1                                                             
	sigma(iv,kv)=alfa(jv,iv)*dndphi(jv,kv)+sigma(iv,kv)
                      
  300 continue                                                           

  200 continue 
	write(7,*)(sigma(jv, 1),jv=1,nelem)

 	stop  
c
      end 
	
	
	
	SUBROUTINE MATINV(A,Nelem,INDEX1,IPIVOT,PIVOT)                    


	 COMPLEX A(480,480),
     &        AMAX,SWAP,T,
     &        PIVOT(480)
	integer Nelem,index1(480,2),
     &        j,i,k,irow,icolum,l,l1,jrow,jcolum,
     &        ipivot(480)
c      DIMENSION A(N,N),INDEX(N,2),PIVOT(N),IPIVOT(N)                   MAT00150
C
      num=nelem
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