subroutine bvalus(n,knot,coef,nk,x,s,order)
double precision knot(1),coef(1),x(1),s(1)
double precision bvalue
integer order,nk
do i=1,n{s(i)=bvalue(knot,n+4,coef,nk,4,x(i),order)}
return
end
