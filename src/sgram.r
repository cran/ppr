subroutine sgram(sg0,sg1,sg2,sg3,tb,nb)
implicit double precision(a-h,o-z) 
integer   nb,ileft,ilo,mflag,
	  i,ii,jj		# indices
double precision   sg0(nb),sg1(nb),sg2(nb),sg3(nb),tb(nb+4),
		 vnikx(4,3),work(16),yw1(4),yw2(4),
		 wpt



 		#PURPOSE

# 	Calculation of the cubic B-spline smoothness prior
# 	for "usual" interior knot setup.


#	Uses BSPVD and INTRV in the CMLIB




#	sgm[0-3](nb)    Symmetric matrix
#                       whose (i,j)'th element contains the integral of
#                       B''(i,.) B''(j,.) , i=1,2 ... nb and j=i,...nb.
#                       Only the upper four diagonals are computed.


lentb=nb+4
	#Initialise the sigma vectors
		do i=1,nb{ sg0(i)=0.;sg1(i)=0.;sg2(i)=0.;sg3(i)=0.}
	
		ilo = 1

	do i=1,nb {

                # Calculate a linear approximation to the
		# second derivative of the non-zero B-splines
		# over the interval [tb(i),tb(i+1)].

			#call intrv(tb(1),(nb+1),tb(i),ilo,ileft,mflag)
			call interv(tb(1),(nb+1),tb(i),ileft,mflag)


		        #Left end second derivatives
			#call bspvd (tb,4,3,tb(i),ileft,4,vnikx,work)
			call bsplvd (tb,lentb,4,tb(i),ileft,work,vnikx,3)

		        # Put values into yw1
			do ii=1,4 { yw1(ii) = vnikx(ii,3) }

	
		        #Right end second derivatives
			#call bspvd (tb,4,3,tb(i+1),ileft,4,vnikx,work)
			call bsplvd (tb,lentb,4,tb(i+1),ileft,work,vnikx,3)


		# Slope*(length of interval) in Linear Approximation to B''
			do ii=1,4 { yw2(ii) = vnikx(ii,3) - yw1(ii)  }
		        wpt = tb(i+1) - tb(i)


# Calculate Contributions to the simga vectors

if(ileft>=4){
do ii=1,4 {

jj=ii
sg0(ileft-4+ii) = sg0(ileft-4+ii) +
wpt* (yw1(ii)*yw1(jj) + (yw2(ii)*yw1(jj) + yw2(jj)*yw1(ii))*.50  +
yw2(ii)*yw2(jj)*.3330 )
jj=ii+1
if(jj<=4) {sg1(ileft+ii-4) = sg1(ileft+ii-4) +
wpt* (yw1(ii)*yw1(jj) + (yw2(ii)*yw1(jj) + yw2(jj)*yw1(ii))*.50  +
yw2(ii)*yw2(jj)*.3330 )}
jj=ii+2
if(jj<=4) {sg2(ileft+ii-4) = sg2(ileft+ii-4) +
wpt* (yw1(ii)*yw1(jj) + (yw2(ii)*yw1(jj) + yw2(jj)*yw1(ii))*.50  +
yw2(ii)*yw2(jj)*.3330 )}
jj=ii+3
if(jj<=4) {sg3(ileft+ii-4) = sg3(ileft+ii-4) +
wpt* (yw1(ii)*yw1(jj) + (yw2(ii)*yw1(jj) + yw2(jj)*yw1(ii))*.50  +
yw2(ii)*yw2(jj)*.3330 )}

          }
	  }

else if(ileft==3){
do ii=1,3 {

jj=ii
sg0(ileft-3+ii) = sg0(ileft-3+ii) +
wpt* (yw1(ii)*yw1(jj) + (yw2(ii)*yw1(jj) + yw2(jj)*yw1(ii))*.50  +
yw2(ii)*yw2(jj)*.3330 )
jj=ii+1
if(jj<=3) {sg1(ileft+ii-3) = sg1(ileft+ii-3) +
wpt* (yw1(ii)*yw1(jj) + (yw2(ii)*yw1(jj) + yw2(jj)*yw1(ii))*.50  +
yw2(ii)*yw2(jj)*.3330 )}
jj=ii+2
if(jj<=3) {sg2(ileft+ii-3) = sg2(ileft+ii-3) +
wpt* (yw1(ii)*yw1(jj) + (yw2(ii)*yw1(jj) + yw2(jj)*yw1(ii))*.50  +
yw2(ii)*yw2(jj)*.3330 )}

          }
	  }

else if(ileft==2){
do ii=1,2 {

jj=ii
sg0(ileft-2+ii) = sg0(ileft-2+ii) +
wpt* (yw1(ii)*yw1(jj) + (yw2(ii)*yw1(jj) + yw2(jj)*yw1(ii))*.50  +
yw2(ii)*yw2(jj)*.3330 )
jj=ii+1
if(jj<=2) {sg1(ileft+ii-2) = sg1(ileft+ii-2) +
wpt* (yw1(ii)*yw1(jj) + (yw2(ii)*yw1(jj) + yw2(jj)*yw1(ii))*.50  +
yw2(ii)*yw2(jj)*.3330 )}

          }
	  }

else if(ileft==1){
do ii=1,1 {

jj=ii
sg0(ileft-1+ii) = sg0(ileft-1+ii) +
wpt* (yw1(ii)*yw1(jj) + (yw2(ii)*yw1(jj) + yw2(jj)*yw1(ii))*.50  +
yw2(ii)*yw2(jj)*.3330 )

          }
	  }}


return
end
