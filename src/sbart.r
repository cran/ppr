subroutine sbart(penalt,dofoff,xs,ys,ws,n,knot,nk,
		  coef,sz,lev,
		  crit,icrit,spar,ispar,lspar,uspar,tol,
		  isetup,
		  xwy,
		  hs0,hs1,hs2,hs3,sg0,sg1,sg2,sg3,
		  abd,p1ip,p2ip,ld4,ldnk,ier)
implicit double precision(a-h,o-z) 




       # A Cubic B-spline Smoothing routine.

#
#          The algorithm minimises:
#
#      (1/n) * sum ws(i)**2 * (ys(i)-sz(i))**2 + lambda* int ( sz"(xs) )**2 dxs
#
#        lambda is a function of the spar which is assumed to be between
#        0 and 1


	      # Input
#penalt 		A penalty > 1 to be used in the gcv criterion

#   n               number of data points
#  ys(n)	    vector of length n containing the observations
#  ws(n)            vector containing the weights given to each data point
#  xs(n)            vector containing the ordinates of the observations
 

#  nk               number of b-spline coefficients to be estimated
#                   nk <= n+2
#  knot(nk+4)       vector of knot points defining the cubic b-spline basis.
# 		    To obtain full cubic smoothing splines one might
#		    have (provided the xs-values are strictly increasing)


#  spar             penalised likelihood smoothing parameter
#  ispar            indicator saying if spar is supplied or to be estimated
#  lspar, uspar     lower and upper values for spar 0.,1. are good values
#  tol              used in Golden Search routine

#  isetup           setup indicator

#  icrit            indicator saying which cross validation score
#		    is to be computed

#  ld4              the leading dimension of abd (ie ld4=4)
#  ldnk             the leading dimension of p2ip (not referenced)


                # Output

#   coef(nk)       vector of spline coefficients
#   sz(n)          vector of smoothed z-values
#   lev(n)         vector of leverages
#   crit           either ordinary of generalized CV score
#   ier            error indicator
#                  ier = 0 ___  everything fine
#                  ier = 1 ___  spar too small or too big
#                               problem in cholesky decomposition



         # Working arrays/matrix
#   xwy		      X'Wy
#   hs0,hs1,hs2,hs3   the diagonals of the X'WX matrix
#   sg0,sg1,sg2,sg3   the diagonals of the Gram matrix
#   abd(ld4,nk)       [ X'WX+lambda*SIGMA] in diagonal form
#   p1ip(ld4,nk)       inner products between columns of L inverse
#   p2ip(ldnk,nk)      all inner products between columns of L inverse
#                          L'L = [X'WX+lambdaSIGMA]  NOT REFERENCED

integer         n,nk,isetup,icrit,ispar,ld4,ldnk,ier
double precision penalt,dofoff,xs(n),ys(n),ws(n),
	 	knot(nk+4),
	 	coef(nk),sz(n),lev(n),
	 	crit,spar,lspar,uspar,tol,
		xwy(nk),
	 	hs0(nk),hs1(nk),hs2(nk),hs3(nk),
         	sg0(nk),sg1(nk),sg2(nk),sg3(nk),
	 	abd(ld4,nk),p1ip(ld4,nk),p2ip(ldnk,nk)






		  # Local variables

double precision		  t1,t2,ratio,
                  a,b,c,d,e,eps,xm,p,q,r,tol1,tol2,u,v,w,
                  fu,fv,fw,fx,x,
		  ax,bx

integer           i





     #  Compute SIGMA, X' W**2 X, X' W**2 z, trace ratio, s0, s1.

                     # SIGMA     -> sg0,sg1,sg2,sg3
                     # X' W**2 X -> hs0,hs1,hs2,hs3
     	             # X' W**2 Z -> xwy

# trevor fixed this 4/19/88
# Note: sbart uses the square of the weights
# the following rectifies that
	for(i=1;i<=n;i=i+1)
		{ if (ws(i)>0) ws(i)=sqrt(ws(i)) }
# 
#
	   if(isetup==0){
           call sgram(sg0,sg1,sg2,sg3,knot,nk)
           call stxwx(xs,ys,ws,n,knot,nk,xwy,hs0,hs1,hs2,hs3)


		    t1=0. ; t2=0.
		    do i=3,nk-3 { t1 = t1 + hs0(i) }
		    do i=3,nk-3 { t2 = t2 + sg0(i) }
		    ratio = t1/t2

		    isetup = 1 }



     # Compute estimate


     if(ispar==1) { # Value of spar supplied

		call sslvrg(penalt,dofoff,xs,ys,ws,n,knot,nk,
		  		coef,sz,lev,crit,icrit,
		  		spar,ratio,
				xwy,
		  		hs0,hs1,hs2,hs3,
		  		sg0,sg1,sg2,sg3,
		  		abd,p1ip,p2ip,ld4,ldnk,ier)
	# check 2
	# print 2999;2999 format(" got through check 2")

		    return }

     else         {
     
     # Use Forsythe Malcom and Moler routine to minimise criterion
     
      ax=lspar ; bx=uspar  # f denotes the value of the criterion

#
#      an approximation  x  to the point where	f  attains a minimum  on
#  the interval  (ax,bx)  is determined.
#
#
#  input..
#
#  ax	 left endpoint of initial interval
#  bx	 right endpoint of initial interval
#  f	 function subprogram which evaluates  f(x)  for any  x
#	 in the interval  (ax,bx)
#  tol	 desired length of the interval of uncertainty of the final
#	 result ( .ge. 0.0)
#
#
#  output..
#
#  fmin  abcissa approximating the point where	f  attains a minimum
#
#
#      the method used is a combination of  golden  section  search  and
#  successive parabolic interpolation.	convergence is never much slower
#  than  that  for  a  fibonacci search.  if  f  has a continuous second
#  derivative which is positive at the minimum (which is not  at  ax  or
#  bx),  then  convergence  is	superlinear, and usually of the order of
#  about  1.324....
#      the function  f	is never evaluated at two points closer together
#  than  eps*abs(fmin) + (tol/3), where eps is	approximately the square
#  root  of  the  relative  machine  precision.   if   f   is a unimodal
#  function and the computed values of	 f   are  always  unimodal  when
#  separated by at least  eps*abs(x) + (tol/3), then  fmin  approximates
#  the abcissa of the global minimum of  f  on the interval  ax,bx  with
#  an error less than  3*eps*abs(fmin) + tol.  if   f	is not unimodal,
#  then fmin may approximate a local, but perhaps non-global, minimum to
#  the same accuracy.
#      this function subprogram is a slightly modified	version  of  the
#  algol  60 procedure	localmin  given in richard brent, algorithms for
#  minimization without derivatives, prentice - hall, inc. (1973).
#
#
#      double precision  a,b,c,d,e,eps,xm,p,q,r,tol1,tol2,u,v,w
#      double precision  fu,fv,fw,fx,x
#
#  c is the squared inverse of the golden ratio
#
      c = 0.5*(3. - sqrt(5e0))
#
#  eps is approximately the square root of the relative machine
#  precision.
#
      eps = 1e0
   10 eps = eps/2e0
      tol1 = 1e0 + eps
      if (tol1 .gt. 1e0) go to 10
      eps = sqrt(eps)
      eps=.000244

#
#  initialization
#
      a = ax
      b = bx
      v = a + c*(b - a)
      w = v
      x = v
      e = 0.0

		spar = x
		call sslvrg(penalt,dofoff,xs,ys,ws,n,knot,nk,
		  		coef,sz,lev,crit,icrit,
		  		spar,ratio,
				xwy,
		  		hs0,hs1,hs2,hs3,
		  		sg0,sg1,sg2,sg3,
		  		abd,p1ip,p2ip,ld4,ldnk,ier)

      fx = crit
      fv = fx
      fw = fx
#
#  main loop starts here
#
   20 xm = 0.5*(a + b)
      tol1 = eps*abs(x) + tol/3e0
      tol2 = 2e0*tol1
#
#  check stopping criterion
#
      if (abs(x - xm) .le. (tol2 - 0.5*(b - a))) go to 90
#
# is golden-section necessary
#
      if (abs(e) .le. tol1) go to 40
#
#  fit parabola
#
      r = (x - w)*(fx - fv)
      q = (x - v)*(fx - fw)
      p = (x - v)*q - (x - w)*r
      q = 2.00*(q - r)
      if (q .gt. 0.0) p = -p
      q =  abs(q)
      r = e
      e = d
#
#  is parabola acceptable
#
   30 if (abs(p) .ge. abs(0.5*q*r)) go to 40
      if (p .le. q*(a - x)) go to 40
      if (p .ge. q*(b - x)) go to 40
#
#  a parabolic interpolation step
#
      d = p/q
      u = x + d
#
#  f must not be evaluated too close to ax or bx
#
      if ((u - a) .lt. tol2) d = sign(tol1, xm - x)
      if ((b - u) .lt. tol2) d = sign(tol1, xm - x)
      go to 50
#
#  a golden-section step
#
   40 if (x .ge. xm) e = a - x
      if (x .lt. xm) e = b - x
      d = c*e
#
#  f must not be evaluated too close to x
#
   50 if (abs(d) .ge. tol1) u = x + d
      if (abs(d) .lt. tol1) u = x + sign(tol1, d)

		spar = u
		call sslvrg(penalt,dofoff,xs,ys,ws,n,knot,nk,
		  		coef,sz,lev,crit,icrit,
		  		spar,ratio,
				xwy,
		  		hs0,hs1,hs2,hs3,
		  		sg0,sg1,sg2,sg3,
		  		abd,p1ip,p2ip,ld4,ldnk,ier)

      fu = crit
#
#  update  a, b, v, w, and x
#
      if (fu .gt. fx) go to 60
      if (u .ge. x) a = x
      if (u .lt. x) b = x
      v = w
      fv = fw
      w = x
      fw = fx
      x = u
      fx = fu
      go to 20
   60 if (u .lt. x) a = u
      if (u .ge. x) b = u
      if (fu .le. fw) go to 70
      if (w .eq. x) go to 70
      if (fu .le. fv) go to 80
      if (v .eq. x) go to 80
      if (v .eq. w) go to 80
      go to 20
   70 v = w
      fv = fw
      w = u
      fw = fu
      go to 20
   80 v = u
      fv = fu
      go to 20
#
#  end of main loop
#
   90 continue ; spar = x ; crit = fx
      return         }





return
end
