subroutine sinerp(abd,ld4,nk,p1ip,p2ip,ldnk,flag)
implicit double precision(a-h,o-z) 
integer	flag,ld4,nk,ldnk,i,j,k
double precision abd(ld4,nk),p1ip(ld4,nk),p2ip(ldnk,nk),
	wjm3(3),wjm2(2),wjm1(1),c0,c1,c2,c3





	# Purpose :  Computes Inner Products between columns of L(-1)
	#	     where L = abd is a Banded Matrix with 3 subdiagonals

	# The algorithm works in two passes:
	#
	#		Pass 1 computes (cj,ck) k=j,j-1,j-2,j-3 ,j=nk, .. 1
	#		Pass 2 computes (cj,ck) k<=j-4  (If flag == 1 ).
	#
	#		A refinement of Elden's trick is used.
	#




			# Pass 1


			wjm3(1)=0e0; wjm3(2)=0e0; wjm3(1)=0e0
			wjm2(1)=0e0; wjm2(2)=0e0
			wjm1(1)=0e0

		do i=1,nk { j=nk-i+1

			             c0 = 1e0/abd(4,j)
			if(j<=nk-3) {c1 = abd(1,j+3)*c0
			             c2 = abd(2,j+2)*c0
				     c3 = abd(3,j+1)*c0 }
		   else if(j==nk-2) {c1 = 0e0
			             c2 = abd(2,j+2)*c0
				     c3 = abd(3,j+1)*c0 }
		   else if(j==nk-1) {c1 = 0e0
			             c2 = 0e0
				     c3 = abd(3,j+1)*c0 }
		   else if(j==nk)   {c1 = 0e0
			             c2 = 0e0
				     c3 = 0e0}


		p1ip(1,j) = 0e0- (c1*wjm3(1)+c2*wjm3(2)+c3*wjm3(3))
		p1ip(2,j) = 0e0- (c1*wjm3(2)+c2*wjm2(1)+c3*wjm2(2))
		p1ip(3,j) = 0e0- (c1*wjm3(3)+c2*wjm2(2)+c3*wjm1(1))

		p1ip(4,j) = c0**2 + 
		c1**2*wjm3(1)+2.*c1*c2*wjm3(2)+2.*c1*c3*wjm3(3) +
		c2**2*wjm2(1)+2.*c2*c3*wjm2(2) +
		c3**2*wjm1(1)

			wjm3(1)=wjm2(1) ; wjm3(2)=wjm2(2) ; wjm3(3)=p1ip(2,j)
			wjm2(1)=wjm1(1) ; wjm2(2)=p1ip(3,j)
			wjm1(1)=p1ip(4,j)

				}


	if(flag==0) {return}


		   # Pass 2


	else	    { # Compute p2ip

			do i=1,nk { j=nk-i+1
			for(k=1;k<=4 & j+k-1<=nk;k=k+1)
			{ p2ip(j,j+k-1) = p1ip(5-k,j) }}
		

			do i=1,nk { j=nk-i+1
			for(k=j-4;k>=1;k=k-1){

			c0 = 1./abd(4,k) ; c1 = abd(1,k+3)*c0
			c2 = abd(2,k+2)*c0 ; c3 = abd(3,k+1)*c0
				    
			p2ip(k,j) = 0e0- ( c1*p2ip(k+3,j)  + 
					     c2*p2ip(k+2,j)  + 
					     c3*p2ip(k+1,j)  )  }
				  }

			return

			}


end
