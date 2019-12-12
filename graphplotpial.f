c*********************************************************************************************    
C																						   *	
c	THIS IS THE PROGRAM FOR  CALCULATION  OF STEADY FREE SURFACE FLOW PATTERN              *
C																						   *	
C	     THIS PROGRAM IS USED FOR PLOTTING GRAPH                                	       *	
C																						   *	
c********************************************************************************************



	open(3,file='output/gnu_fine',status='unknown')
      open(4,file='graph/graph_fine',status='unknown')
				
	

	write(*,*) ' what is the  number of  surface points in one side ?'
	Write(*,*)
	Write(*,*)'FOR PRESENT CASE, PLEASE WRITE  2772'
	read (*,*)nn
	n1=2*nn

	write(*,*) ' what is the total number of points in one row ?'
	Write(*,*)
	Write(*,*)'FOR PRESENT CASE, PLEASE WRITE  21'
	read (*,*)n2
	n=n1/n2

	do 122 j=1,n
	do 121 i=1,n2
	read (3,*)x,y,z,p1
	write(4,1)x,y,z,p1

121	continue
	write(4,*)
122   continue
   
1	format(3(f10.6,2x),3f15.8)
	stop
	end
	

