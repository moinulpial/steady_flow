c*********************************************************************************************    
C																						   *	
c	THIS IS THE PROGRAM FOR THE CALCULATION  OF STEADY FREE SURFACE FLOW PATTERN           *
C																						   *	
C	     THIS PROGRAM GENERATES ELEVATIONS BY USING THE INPUT FROM PHI.F	    	       *	
C																						   *	
c********************************************************************************************

	
	
	dimension x(11000),dx(11000),dy(11000),dz(11000),pp1(11000)
	dimension y(11000),z(11000),p1(11000),p2(11000),pp2(11000)
	1,p11(11000),p12(11000)
	

      open(3,file='output/elevtn1',status='unknown')
	open(4,file='output/gnu_fine',status='unknown')
	


	pi=3.14
	pi4=1/(4.0*pi)

	write(*,*) ' what is the number of surface element?'
	write(*,*)
	Write(*,*)'FOR PRESENT CASE, PLEASE WRITE  2772'
	write(*,*)
	read (*,*)n
	do 121 i=1,n
	
	read (3,*) x(i),y(i),z(i),p1(i)
	write(4,1)x(i),y(i),z(i),p1(i)
121	continue
c	write(4,*)


	do 125  ii=1,n	
	dx(ii)=x(ii)
	dy(ii)= -y(ii)
	dz(ii)=0.0
	pp1(ii)=p1(ii)
	pp2(ii)=p2(ii)
	
	write(4,1)dx(ii),dy(ii),dz(ii),pp1(ii)

	
125   CONTINUE

	     
1	format(6(f14.6,2x))
	stop
	end
	

