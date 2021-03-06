subroutine sslvrg(penalt,dofoff,x,y,w,n,knot,nk,
		  coef,sz,lev,
		  crit,icrit,
		  spar,ratio,
		  xwy,
		  hs0,hs1,hs2,hs3,
		  sg0,sg1,sg2,sg3,
		  abd,p1ip,p2ip,ld4,ldnk,info)

implicit double precision(a-h,o-z) 
integer  n,nk,icrit,ld4,ldnk,i,icoef,ileft,ilo,info,j,mflag
double precision  penalt,dofoff,x(n),y(n),w(n),
	 knot(nk+4),
	 coef(nk),sz(n),lev(n),
	 crit,
	 ratio,spar,
	 xwy(nk),
	 hs0(nk),hs1(nk),hs2(nk),hs3(nk),
         sg0(nk),sg1(nk),sg2(nk),sg3(nk),
	 abd(ld4,nk),p1ip(ld4,nk),p2ip(ldnk,nk),
	 lambda,
	 b0,b1,b2,b3,eps,
	 vnikx(4,1),work(16),
	# xv,bvalu,rss,df
	 xv,bvalue,rss,df


lenkno=nk+4

         ilo = 1 ; eps = .1e-10

    # Purpose : Solves the smoothing problem and computes the
    #           criterion function (OCV or GCV).


	# The coeficients of estimated smooth

	                      lambda = ratio*16.**(-2. + spar*(6.))

          do i=1,nk     { coef(i)    = xwy(i) }
          do i=1,nk     { abd(4,i)   = hs0(i)+lambda*sg0(i) }
          do i=1,(nk-1) { abd(3,i+1) = hs1(i)+lambda*sg1(i) }
          do i=1,(nk-2) { abd(2,i+2) = hs2(i)+lambda*sg2(i) }
          do i=1,(nk-3) { abd(1,i+3) = hs3(i)+lambda*sg3(i) }


	  call dpbfa(abd,ld4,nk,3,info)
	  if(info.ne.0) { return } 
	  call dpbsl(abd,ld4,nk,3,coef)


	# Value of smooth at the data points

                     icoef = 1
          do i=1,n {    xv = x(i)
			 #sz(i) = bvalu(knot,coef,
			 #               nk,4,0,xv,icoef,work(1)) }
			 sz(i) = bvalue(knot,lenkno,coef,
			                nk,4,xv,0)

						 }


	# Compute the criterion function if requested

	  if(icrit==0) { return}

	  else         { # Ordinary or Generalized CV

	                 # Get Leverages First

	      call sinerp(abd,ld4,nk,p1ip,p2ip,ldnk,0)

     do i=1,n { xv = x(i)
	        #call intrv(knot(1),(nk+1),xv,ilo,ileft,mflag)
	        call interv(knot(1),(nk+1),xv,ileft,mflag)

	        if(mflag==-1) { ileft = 4   ; xv = knot(4)+eps }
	        if(mflag==1)  { ileft = nk  ; xv = knot(nk+1)-eps }
 	        j=ileft-3
                #call bspvd(knot,4,1,xv,ileft,4,vnikx,work)
                call bsplvd(knot,lenkno,4,xv,ileft,work,vnikx,1)

	        b0=vnikx(1,1);b1=vnikx(2,1);b2=vnikx(3,1);b3=vnikx(4,1)

	lev(i) = (p1ip(4,j)*b0**2   + 2.*p1ip(3,j)*b0*b1 +
				   2.*p1ip(2,j)*b0*b2 + 2.*p1ip(1,j)*b0*b3 +
		  p1ip(4,j+1)*b1**2 + 2.*p1ip(3,j+1)*b1*b2 +
				   2.*p1ip(2,j+1)*b1*b3 +
		  p1ip(4,j+2)*b2**2 + 2.*p1ip(3,j+2)*b2*b3 +
		  p1ip(4,j+3)*b3**2 )*w(i)**2      }


			  


	# Evaluate Criterion

	      if(icrit==1) { # Generalized CV
	       
			                rss = 0e0 ; df = 0e0 ; sumw=0e0
			     do i=1,n { rss = rss + ((y(i)-sz(i))*w(i))**2}
			     do i=1,n { df  = df  + lev(i)}
#			     do i=1,n { sumw  = sumw  + w(i)}
		
			     crit = (rss/n)/((1e0-(dofoff+penalt*df)/n)**2)
#			     crit = (rss/n)/((1e0-(dofoff+penalt*df)/sumw)**2)
				}

	      else if(icrit==2)        { # Ordinary CV
	       
			                crit = 0e0
			     do i=1,n { crit = crit +
			                (((y(i)-sz(i))*w(i))/(1-lev(i)))**2 }
					crit=crit/n}
		else		{ #df matching
				crit=0e0
				do i=1,n {crit=crit+lev(i)}
				crit=3+(dofoff-crit)**2}


                        return }

end
