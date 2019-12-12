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
c	THIS IS THE FIFTH PROGRAM FOR THE CALCULATION  OF STEADY FREE SURFACE FLOW PATTERN     *
C																						   *	
C	       THIS PROGRAM TAKES INPUT FROM  SURFACE_GRN.F, SVGREN_SURF AND GGFILE.F          *
C																						   *		
C		            	IT GENERATES  FREE SURFACE ELEVATION             				   *	
c*********************************************************************************************
	integer nelem
	real x(9480),y(9480),z(9480),ds(500),eliv(9480)     		

      complex beta(9480,1480),dxbeta(9480,1480),dxdphi(9480,7),
     &        phi(9480,7),graph(9480),
     &        sigma(480,7)



	open(7,file='output/sigma_body',status='unknown')
	open(2,file='output/area',status='unknown')

	open(6,file='data/fine_surf2',status='unknown')
	open(4,file='output/grn_surf',status='unknown')
	open(44,file='output/dxgrn1',status='unknown')

	open(1,file='output/elevtn1',status='unknown')


	Write(*,*)'how many points in surface'
	
	Write(*,*)
	Write(*,*)'FOR PRESENT CASE, PLEASE WRITE  2772'
	read(*,*)nsurp
	Write(*,*)'how many elements in body'
	
	Write(*,*)
	Write(*,*)'FOR PRESENT CASE, PLEASE WRITE  250'
	read(*,*)nelem




cc	fn= froude number
cc	alg= length of the ship
CC	DOMEGA= FREQUENCY OF OSCILLATION




	
	fn=0.2670
	alg=2.0
	pi=3.14
	pi4=1/(4.0*pi)
	grav=9.8
	OMEGA=.001
	num=nelem

	speed=fn*sqrt(grav*alg)
	
	Write(*,*)SPEED

	do 101 i=1,nsurp
	read(6,*)x(i),y(i),z(i)
101	continue



	read(7,*)(sigma(jv, 1),jv=1,nelem)
	read(2,*)(ds(i),i=1,nelem)
	read(4,*)((beta(i,j),j=1,nelem),i=1,nsurp)
	read(44,*)((dxbeta(i,j),j=1,nelem),i=1,nsurp)



	do 1000 iv=1,1                                                                                                                               
	do 1000 jv=1,nsurp
      phi(jv,iv) =cmplx(0.0,0.0)
C	eliv(jv)  =	cmplx(0.0,0.0)
	eliv(jv)  =	0.0


      do 1000 kv=1,nelem
      phi(jv,iv)   =(pi4*beta(jv,kv)*sigma(kv,iv)*ds(kv))+phi(jv,iv)
	dxdphi(jv,iv)=(pi4*dxbeta(jv,kv)*sigma(kv,iv)*ds(kv))
	1+dxdphi(jv,iv)

	eliv(jv)=-(speed/grav)*real(dxdphi(jv,iv))


 1000 continue
	
	do 1001 jv=1,nsurp 
	write(1,11) x(jv),y(jv),z(jv),eliv(jv)
1001	continue

11	format(3f15.9,4x,2f20.9)

	
 	stop  
c
      end 

