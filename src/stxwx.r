subroutine stxwx(x,z,w,k,xknot,n,y,hs0,hs1,hs2,hs3)
implicit double precision(a-h,o-z) 
integer          k,n,
		 j,i,ilo,ileft,mflag       # local
double precision             z(k),w(k),x(k),
		 xknot(n+4),
		 y(n),
		 hs0(n),hs1(n),hs2(n),hs3(n),
                
		 eps,vnikx(4,1),work(16)     # local



lenxk=n+4

     # Initialise the output vectors

     do i=1,n { y(i)=0e0 ; hs0(i)=0e0 ; hs1(i)=0e0
			  hs2(i)=0e0 ; hs3(i)=0e0 }


     # Compute X'WX -> hs0,hs1,hs2,hs3  and X'WZ -> y

	ilo=1  ; eps = .1e-9

     do i=1,k {

	      #call intrv(xknot(1),(n+1),x(i),ilo,ileft,mflag)
	      call interv(xknot(1),(n+1),x(i),ileft,mflag)


#	      if(mflag==-1) {write(6,'("Error in hess ",i2)')mflag;stop}
#	      if(mflag==-1) {return}
	      if(mflag== 1) {if(x(i)<=(xknot(ileft)+eps)){ileft=ileft-1}
#			   else{write(6,'("Error in hess ",i2)')mflag;stop}}
			   else{return}}


              #call bspvd (xknot,4,1,x(i),ileft,4,vnikx,work)
              call bsplvd (xknot,lenxk,4,x(i),ileft,work,vnikx,1)
	# check 2
	# print 2999;2999 format(" got through check2 stxwx ")



	      j= ileft-4+1
	       y(j) =  y(j)+w(i)**2*z(i)*vnikx(1,1)
	      hs0(j)=hs0(j)+w(i)**2*vnikx(1,1)**2
	      hs1(j)=hs1(j)+w(i)**2*vnikx(1,1)*vnikx(2,1)
	      hs2(j)=hs2(j)+w(i)**2*vnikx(1,1)*vnikx(3,1)
	      hs3(j)=hs3(j)+w(i)**2*vnikx(1,1)*vnikx(4,1)

	      j= ileft-4+2
	       y(j) =  y(j)+w(i)**2*z(i)*vnikx(2,1)
	      hs0(j)=hs0(j)+w(i)**2*vnikx(2,1)**2
	      hs1(j)=hs1(j)+w(i)**2*vnikx(2,1)*vnikx(3,1)
	      hs2(j)=hs2(j)+w(i)**2*vnikx(2,1)*vnikx(4,1)

	      j= ileft-4+3
	       y(j) =  y(j)+w(i)**2*z(i)*vnikx(3,1)
	      hs0(j)=hs0(j)+w(i)**2*vnikx(3,1)**2
	      hs1(j)=hs1(j)+w(i)**2*vnikx(3,1)*vnikx(4,1)

	      j= ileft-4+4
	       y(j) =  y(j)+w(i)**2*z(i)*vnikx(4,1)
	      hs0(j)=hs0(j)+w(i)**2*vnikx(4,1)**2  }

return
end
