real(8) function PARTF(i, iso, T)
 use PART_FUNCTION
 implicit none
 integer, intent(in) :: i, iso
 real(8), intent(in) :: T

 integer :: lo, hi, mid
 real(8) :: t1, t2, p1, p2, w

 if(T <= tpf(1)) then
	PARTF = pf(i,iso,1)
	return
 endif
 if(T >= tpf(301)) then
	PARTF = pf(i,iso,301)
	return
 endif
 lo = 2; hi = 301
 do while(lo < hi)
	mid = (lo + hi)/2
	if(T < tpf(mid)) then
		hi = mid
	else
		lo = mid + 1
	endif
 enddo
 t1 = tpf(hi-1); t2 = tpf(hi)
 p1 = pf(i,iso,hi-1); p2 = pf(i,iso,hi)
 if(t2 == t1) then
	PARTF = p1
 else
	w = (T - t1)/(t2 - t1)
	PARTF = (1.d0 - w)*p1 + w*p2
 endif
end function PARTF


integer(4) function INDWAV(wc,v)
 use OBSERVED; use CONSTANTS
 implicit none
 integer    :: j
 real(8)    :: wc, v
 INDWAV = 0
 wc = wc*(1.d0+v/cvel)
 do j = 1,nobs-1
	if(wc > wobs(j) .and. wc <= wobs(j+1)) then
		if(wc < wobs(j)+(wobs(j+1)-wobs(j))/2) INDWAV = j
		if(wc >= wobs(j)+(wobs(j+1)-wobs(j))/2) INDWAV = j+1
		exit
	endif
 enddo
end function INDWAV


subroutine QINTER(xold, fold, nold, xnew, fnew, nnew)
 implicit none
 integer, intent(in) :: nold, nnew
 real(8), intent(in) :: xold(nold), fold(nold), xnew(nnew)
 real(8), intent(out) :: fnew(nnew)
 integer :: k, l, i0, i1, i2
 real(8) :: x, x0, x1, x2, f0, f1, f2
 real(8) :: d01, d02, d12
 real(8) :: w0, w1, w2
 if(nold < 2) then
	fnew = 0.d0
	return
 endif
 l = 2
 do k = 1, nnew
	x = xnew(k)
	if(k > 1) then
		if(xnew(k) < xnew(k-1)) l = 2
	endif
	if(x <= xold(1)) then
		if(xold(2) /= xold(1)) then
			fnew(k) = fold(1) + (fold(2)-fold(1))*(x-xold(1))/(xold(2)-xold(1))
		else
			fnew(k) = fold(1)
		endif
		cycle
	endif
	if(x >= xold(nold)) then
		if(xold(nold) /= xold(nold-1)) then
			fnew(k) = fold(nold-1) + (fold(nold)-fold(nold-1))*(x-xold(nold-1))/(xold(nold)-xold(nold-1))
		else
			fnew(k) = fold(nold)
		endif
		cycle
	endif
	do while (l < nold .and. x >= xold(l))
		l = l + 1
	enddo
	i1 = l-1
	if(i1 <= 1 .or. i1 >= nold) then
		if(xold(i1+1) /= xold(i1)) then
			fnew(k) = fold(i1) + (fold(i1+1)-fold(i1))*(x-xold(i1))/(xold(i1+1)-xold(i1))
		else
			fnew(k) = fold(i1)
		endif
		cycle
	endif
	i0 = i1 - 1
	i2 = i1 + 1
	x0 = xold(i0); f0 = fold(i0)
	x1 = xold(i1); f1 = fold(i1)
	x2 = xold(i2); f2 = fold(i2)
	d01 = x0 - x1
	d02 = x0 - x2
	d12 = x1 - x2
	if(d01 == 0.d0 .or. d02 == 0.d0 .or. d12 == 0.d0) then
		if(xold(i1+1) /= xold(i1)) then
			fnew(k) = fold(i1) + (fold(i1+1)-fold(i1))*(x-xold(i1))/(xold(i1+1)-xold(i1))
		else
			fnew(k) = fold(i1)
		endif
		cycle
	endif
	w0 = (x - x1)*(x - x2) / (d01*d02)
	w1 = (x - x0)*(x - x2) / ((-d01)*d12)
	w2 = (x - x0)*(x - x1) / ((-d02)*(-d12))
	fnew(k) = w0*f0 + w1*f1 + w2*f2
 enddo
end subroutine QINTER

subroutine LININT_UNIFORM(x, x0, dx, n, f, y, ok)
 implicit none
 integer :: n, i0
 real(8) :: x, x0, dx, f(n), y, t
 logical :: ok

 if(dx <= 0.d0) then
	ok = .false.
	y = 0.d0
	return
 endif
 t = (x - x0)/dx
 i0 = int(t) + 1
 if(i0 < 1 .or. i0 >= n) then
	ok = .false.
	y = 0.d0
	return
 endif
 t = t - dble(i0-1)
 y = (1.d0 - t)*f(i0) + t*f(i0+1)
 ok = .true.
end subroutine LININT_UNIFORM

subroutine LINTER ( xold, yold, nold, xnew, ynew, nnew )
  real(8) xold(*), yold(*), xnew(*), ynew(*)
! XOLD AND XNEW INCREASING
  iold = 2
  do inew = 1, nnew
     do while( iold < nold .and. xold(iold) < xnew(inew) )
        iold = iold + 1
     enddo
     ynew(inew) = yold(iold-1) + (yold(iold) - yold(iold-1)) / &
	              (xold(iold) - xold(iold-1)) * (xnew(inew) - xold(iold-1))
  enddo
end subroutine LINTER



subroutine MINI(xopt,obs,weight,n,nwl,fun,sw,nonneg_flag)
! nonneg_flag = 1 enables non-negativity x(i) >= 0
real(8) xopt(n), obs(nwl), weight(nwl)
real(8) ftoler, freduc, fvalue, fopt, dev, oldfun, d, regul, cormax, corlim, corscl, di, qlambda, qscale
logical failure, done
integer info, wrksiz, iterf, maxitr, iterc, nspec, sw
parameter (ftoler = 0.0001d0, freduc = 0.9999d0, corlim = 3.d0)
real(8), allocatable :: x(:), prf(:), work(:), sigc(:,:), diag(:), rhsold(:), deriva(:,:), derivc(:,:), rhs(:)
integer, allocatable :: ipiv(:)
integer nonneg_flag

allocate (x(n), prf(nwl), ipiv(n), work(45*n), sigc(n,n), diag(n), rhsold(n), deriva(nwl,n), derivc(n,n), rhs(n), stat=ios)
if(ios /= 0)stop 'Cannot allocate arrays in MINI subroutine'
work=0.d0
iterc=0; deriva=0.d0; wrksiz=45*n; maxitr=50; ncalls=0; failure=.false.; qlambda=10.0; qscale=10.0; fopt=1.d30; iterf=iterc

open(6,carriagecontrol = 'fortran')
do 
 iterc = iterc + 1

 call DCOPY(n,xopt,1,x,1)
 call APPLY_NNORM(x,n,nonneg_flag)
 call APPLY_NNORM(xopt,n,nonneg_flag)

 oldfun = fopt

 call FUN(x,prf,deriva,n,nwl,'grad')

 do i = 1, n
  do iwl = 1, nwl
   deriva(iwl,i) = deriva(iwl,i) * sqrt(weight(iwl) / nwl)
  enddo
 enddo  

 call DEVIAT(obs,prf,weight,nwl,d,dev)
 fvalue = d

 if(fvalue < fopt) then
  call DCOPY(n,x,1,xopt,1)
  fopt = fvalue
 endif

 call DSYRK('u','t',n,nwl,1.d0,deriva,nwl,0.d0,derivc,n)
 do i1 = 2, n
  do i2 = 1, i1-1
   derivc(i1,i2) = derivc(i2,i1)
  enddo
 enddo

 do i = 1, n
  di = 0.d0
  do iwl = 1, nwl
   di = di + sqrt(weight(iwl)) * (obs(iwl) - prf(iwl)) * deriva(iwl,i)
  enddo
  rhs(i) = di / sqrt(dble(nwl))
 enddo

 call DCOPY(n,derivc,n+1,diag,1)
 call DCOPY(n,rhs,1,rhsold,1)

 done = .false.
 ifail = 0
 itry = 0
 oldfun = fvalue
  
 do
  call DSCAL(n,1.d0+qlambda,derivc,n+1)
  ifail = ifail + 1
  itry = itry + 1 

  call DSYSV('u',n,1,derivc,n,ipiv,rhs,n,work,wrksiz,info)

  if(info /= 0) failure = .true.
  if(failure) then
   write(*,'('' Degenerate matrix DERIVC.'','' Cannot solve SLE.'',T95,''*** (MINI  ) ***'')')
   deallocate (x, prf, ipiv, work, sigc, diag, rhsold, deriva, derivc, rhs)
   stop
  endif
  
  cormax = 0.d0
  do i = 1, n
   if(cormax < abs(rhs(i))) cormax = abs(rhs(i))
  enddo

  corscl = 1.d0
  if(cormax > corlim) corscl = corlim / cormax
  do i = 1, n
   x(i) = x(i) + rhs(i) * corscl
  enddo
  call APPLY_NNORM(x,n,nonneg_flag)

  call FUN(x,prf,deriva,n,nwl,'func')
  call DEVIAT(obs,prf,weight,nwl,d,dev)
  fvalue = d; dev = sqrt(dev)

  if(itry == 2 .and. fvalue > oldfun*freduc) then
   if(qscale > 1.) then
    qscale = 0.1
   else
    qscale = 10.
   endif
   qlambda = qlambda * qscale * qscale
   call DCOPY(n,xopt,1,x,1)
   call APPLY_NNORM(x,n,nonneg_flag)
   fvalue = fopt
   done = .false.
  elseif(itry > 2 .and. fvalue > oldfun*freduc) then
   qlambda = qlambda / qscale
   done = .true.   
  else
   qlambda = qlambda * qscale
   if(qlambda > 1.e10) then
    qlambda = 1.e9
	done = .true.
   elseif(qlambda < 1.e-10) then
    qlambda = 1.e-9
	done = .true.
   else
    done = .false.
   endif
  endif

  if(fvalue < fopt*(1.d0-ftoler)) then
   call DEVIAT(obs,prf,weight,nwl,di,dev)
   dev = sqrt(dev)
   call DCOPY(n,x,1,xopt,1)
   call APPLY_NNORM(xopt,n,nonneg_flag)
   fopt = fvalue
   ifail = 0
  endif

  do i1 = 2, n
   do i2 = 1, i1-1
    derivc(i2,i1) = derivc(i1,i2)
   enddo
  enddo
  call DCOPY(n,diag,1,derivc,n+1)
  call DCOPY(n,rhsold,1,rhs,1)
  
  if(.not. done .or. itry <= 2 .or. fvalue < oldfun*freduc) then
   oldfun = fvalue
   cycle
  else
   exit ! From internal loop
  endif
 
 enddo

 if(sw == 0) write(6,"(1h+'       Correct SLD profiles:            Iter:',I3,'    Mean DEV.=',F7.3,'%')") iterc, dev*100
 if(sw == 1) write(6,"(1h+'      Correct INSTRUMENTAL profiles:   Iter:',I3,'    Mean DEV.=',F7.3,'%')") iterc, dev*100
 if(sw == 2) write(6,"(1h+'      Correct LINE profiles:           Iter:',I3,'    Mean DEV.=',F7.3,'%')") iterc, dev*100
 if(sw == 3) write(6,"(1h+'      Correct ABUNDANCE molecule:           Iter:',I3,'    Mean DEV.=',F7.3,'%')") iterc, dev*100

 if(ifail < itry .or. (iterc - iterf) < 2) then
  call DCOPY(n,xopt,1,x,1)
  call APPLY_NNORM(x,n,nonneg_flag)
  fvalue = fopt
  oldfun = fopt 
  if(iterc >= maxitr) then
 !  write(*,'(''      The limit of iteration number has been reached'',t95,''*** (mini  ) ***'')')
   call DEVIAT(obs,prf,weight,nwl,di,dev)     
   dev = sqrt(dev)
   exit ! From main loop
  endif          
 else
  exit ! From main loop
 endif 
enddo ! End of main loop 
    
call DCOPY(n,xopt,1,x,1)
call APPLY_NNORM(x,n,nonneg_flag)
call DCOPY(n,x,1,xopt,1)

if(sw == 0) write(6,"(1h+'       Correct SLD profiles:            OK          Mean DEV.=',F7.3,'%',/)") dev*100
if(sw == 1) write(6,"(1h+'      Correct INSTRUMENTAL profiles:   OK          Mean DEV.=',F7.3,'%',/)") dev*100
if(sw == 2) write(6,"(1h+'      Correct LINE profiles:           OK          Mean DEV.=',F7.3,'%',/)") dev*100
if(sw == 3) write(6,"(1h+'      Correct ABUNDANCE molecule:           OK          Mean DEV.=',F7.3,'%',/)") dev*100

deallocate (x, prf, ipiv, work, sigc, diag, rhsold, deriva, derivc, rhs)
close(6)

end subroutine MINI


subroutine DEVIAT(obs,prf,weight,nwl,dev1,dev2)
real(8) prf(nwl), obs(nwl), weight(nwl), dev1, dev2
dev1 = 0.d0
dev2 = 0.d0
do iwl = 1, nwl
 dev1 = dev1 + weight(iwl) * (prf(iwl) - obs(iwl))**2
 dev2 = dev2 + (prf(iwl) - obs(iwl))**2
enddo
dev1 = dev1 / nwl
dev2 = dev2 / nwl
end subroutine DEVIAT



! Ax,Bx - min&max interval for x
! F -external function f(x)
! Tol - tolerantnost (1.e-5)
real(8) function FMIN(AX,BX,F,TOL,Fout,number)
implicit real(8) (a-h,o-z)
real(8) Fout
    C=0.5d0*(3.d0-sqrt(5.d0))
    EPS=1.d0
10  EPS=EPS/2.d0
    TOL1=1.d0+EPS
    if(TOL1.GT.1.d0)go to 10
    EPS=sqrt(EPS)
    A=AX; B=BX
    V=A+C*(B-A)
    W=V; X=V; E=0.
    FX=F(X,number)
	FV=FX; FW=FX
20  XM=.5d0*(A+B)
    TOL1=EPS*abs(X)+TOL/3.d0
    TOL2=2.*TOL1
    if(abs(X-XM).LE.(TOL2-.5d0*(B-A)))go to 90
    if(abs(E).LE.TOL1)go to 40
    R=(X-W)*(FX-FV)
    Q=(X-V)*(FX-FW)
    P=(X-V)*Q-(X-W)*R
    Q=2.*(Q-R)
    if(Q.GT.0.)P=-P
    Q=abs(Q)
    R=E; E=D
30  if(abs(P).GE.abs(.5d0*Q*R))go to 40
    if(P.LE.Q*(A-X))go to 40
    if(P.GE.Q*(B-X))go to 40
    D=P/Q; U=X+D
    if((U-A).LT.TOL2)D=sign(TOL1,XM-X)
    if((B-U).LT.TOL2)D=sign(TOL1,XM-X)
    go to 50
40  if(X.GE.XM)E=A-X
    if(X.LT.XM)E=B-X
    D=C*E
50  if(abs(D).GE.TOL1)U=X+D
    if(abs(D).LT.TOL1)U=X+sign(TOL1,D)
    FU=F(U,number)
    if(FU.GT.FX)go to 60
    if(U.GE.X)A=X
    if(U.LT.X)B=X
    V=W; FV=FW; W=X
    FW=FX; X=U; FX=FU
    go to 20
60  if(U.LT.X)A=U
    if(U.GE.X)B=U
    if(FU.LE.FW)go to 70
    if(W.EQ.X)go to 70
    if(FU.LE.FV)go to 80
    if(V.EQ.X)go to 80
    if(V.EQ.W)go to 80
    go to 20
70  V=W; FV=FW
    W=U; FW=FU
    go to 20
80  V=U; FV=FU
    go to 20
90  FMIN=X
	Fout=FU
    return
end function FMIN


subroutine FUNCTION_GAUSS(nfunc,nunknowns,xi,x,f) ! calculate approximation function
 real(8) f(nfunc),xi(nunknowns),x(nfunc)
 do i=1,nfunc
   f(i)=xi(1)*exp(-((x(i)-xi(2))/(xi(3)))**2)
 enddo
end subroutine FUNCTION_GAUSS



subroutine MINI_REG(xopt,absc,obs,weight,regpar,n,nwl,fun,sw,nonneg_flag) ! Tichonov's regularization included
! nonneg_flag = 1 enables non-negativity x(i) >= 0
real(8) xopt(n),obs(nwl),weight(nwl),absc(n),regpar(2),admin,admax,regul 
real(8) ftoler, freduc, fvalue, fopt, dev, oldfun, d, cormax, corlim, corscl, di, qlambda, qscale
logical failure, done
integer info, wrksiz, iterf, maxitr, iterc, nspec, sw
parameter (ftoler = 0.0002d0, freduc = 0.9998d0, corlim = 3.d0)
real(8), allocatable :: x(:), prf(:), work(:), sigc(:,:), diag(:), rhsold(:), deriva(:,:), derivc(:,:), rhs(:), aderiv(:)
integer, allocatable :: ipiv(:)
integer nonneg_flag

allocate (x(n), prf(nwl), ipiv(n), work(45*n), sigc(n,n), diag(n), rhsold(n), deriva(nwl,n), derivc(n,n), rhs(n),aderiv(n),stat=ios)
if(ios /= 0)stop 'Cannot allocate arrays in MINI_REG subroutine'
work=0.d0

iterc=0; deriva=0.d0; wrksiz=45*n; maxitr=50; ncalls=0; failure=.false.; qlambda=10.0; qscale=10.0; fopt=1.d30; iterf=iterc

open(6,carriagecontrol = 'fortran')

do 
 iterc = iterc + 1
 call DCOPY(n,xopt,1,x,1)
 call APPLY_NNORM(x,n,nonneg_flag)
 call APPLY_NNORM(xopt,n,nonneg_flag)
 oldfun = fopt
 
 call FUN(x,prf,deriva,n,nwl,'grad')

 do i = 1, n
  do iwl = 1, nwl
   deriva(iwl,i) = deriva(iwl,i) * sqrt(weight(iwl) / nwl)
  enddo
 enddo  

 call DEVIAT(obs,prf,weight,nwl,d,dev)

  do i=1,n
   aderiv(i)=0.d0
   do iwl=1,nwl    
    aderiv(i)=aderiv(i)+deriva(iwl,i)**2
   enddo
  enddo 
 	 
  aderiv=sqrt(aderiv)

  admin=minval(aderiv)  
  admax=maxval(aderiv)  
  do i=1,n
   aderiv(i)=(aderiv(i)-admin)/(admax-admin) 
  enddo

  call REGFUNQ(x,absc,aderiv,regpar,n,regul)
  fvalue=d+regul 

 if(fvalue < fopt) then
  call DCOPY(n,x,1,xopt,1)
  fopt = fvalue
 endif
 
 call DSYRK('u','t',n,nwl,1.d0,deriva,nwl,0.d0,derivc,n)
 do i1 = 2, n
  do i2 = 1, i1-1
   derivc(i1,i2) = derivc(i2,i1)
  enddo
 enddo

 do i = 1, n
  di = 0.d0
  do iwl = 1, nwl
   di = di + sqrt(weight(iwl)) * (obs(iwl) - prf(iwl)) * deriva(iwl,i)
  enddo
  rhs(i) = di / sqrt(dble(nwl))
 enddo

 call REGDER(x,absc,aderiv,regpar,n,derivc,rhs)
 call DCOPY(n,derivc,n+1,diag,1)
 call DCOPY(n,rhs,1,rhsold,1)

 done = .false.
 ifail = 0
 itry = 0
 oldfun = fvalue
  
 do
  call DSCAL(n,1.d0+qlambda,derivc,n+1)
  ifail = ifail + 1
  itry = itry + 1 
  
  call DSYSV('u',n,1,derivc,n,ipiv,rhs,n,work,wrksiz,info)
  if(info /= 0) failure = .true.

  if(failure) then
   write(*,'('' Degenerate matrix DERIVC.'','' Cannot solve SLE.'',T95,''*** (MINI  ) ***'')')
   write(*,*) 'INFO=', info  
   deallocate (x, prf, ipiv, work, sigc, diag, rhsold, deriva, derivc, rhs)
   stop
  endif
 
  cormax = 0.d0
  do i = 1, n
   if(cormax < abs(rhs(i))) cormax = abs(rhs(i))
  enddo

  corscl = 1.d0
  if(cormax > corlim) corscl = corlim / cormax
  do i = 1, n
   x(i) = x(i) + rhs(i) * corscl
  enddo
  call APPLY_NNORM(x,n,nonneg_flag)

  call FUN(x,prf,deriva,n,nwl,'func')
  call DEVIAT(obs,prf,weight,nwl,d,dev)
  call REGFUNQ(x,absc,aderiv,regpar,n,regul)
  fvalue = d+regul; dev = sqrt(dev)

  if(itry == 2 .and. fvalue > oldfun*freduc) then
   if(qscale > 1.) then
    qscale = 0.1
   else
    qscale = 10.
   endif
   qlambda = qlambda * qscale * qscale
   call DCOPY(n,xopt,1,x,1)
   call APPLY_NNORM(x,n,nonneg_flag)
   fvalue = fopt
   done = .false.
  elseif(itry > 2 .and. fvalue > oldfun*freduc) then
   qlambda = qlambda / qscale
   done = .true.   
  else
   qlambda = qlambda * qscale
   if(qlambda > 1.e10) then
    qlambda = 1.e9
	done = .true.
   elseif(qlambda < 1.e-10) then
    qlambda = 1.e-9
	done = .true.
   else
    done = .false.
   endif
  endif

  if(fvalue < fopt*(1.d0-ftoler)) then
   call DEVIAT(obs,prf,weight,nwl,di,dev)
   dev = sqrt(dev)
   call DCOPY(n,x,1,xopt,1)
   call APPLY_NNORM(xopt,n,nonneg_flag)
   fopt = fvalue
   ifail = 0
  endif

  do i1 = 2, n
   do i2 = 1, i1-1
    derivc(i2,i1) = derivc(i1,i2)
   enddo
  enddo
  call DCOPY(n,diag,1,derivc,n+1)
  call DCOPY(n,rhsold,1,rhs,1)
 
  if(.not. done .or. itry <= 2 .or. fvalue < oldfun*freduc) then
   oldfun = fvalue
   cycle
  else
   exit ! From internal loop
  endif
 enddo

 if(sw == 0) write(6,"(1h+'      Correct SLD profiles:            Iter:',I3,'    Mean DEV.=',F7.3,'%')") iterc, dev*100
 if(sw == 1) write(6,"(1h+'     Correct INSTRUMENTAL profiles:   Iter:',I3,'    Mean DEV.=',F7.3,'%')") iterc, dev*100
 if(sw == 2) write(6,"(1h+'     Correct LINE profiles:           Iter:',I3,'    Mean DEV.=',F7.3,'%')") iterc, dev*100
 if(ifail < itry .or. (iterc - iterf) < 2) then
  call DCOPY(n,xopt,1,x,1)
  call APPLY_NNORM(x,n,nonneg_flag)
  fvalue = fopt
  oldfun = fopt 
 
  if(iterc >= maxitr) then
 !  write(*,'(''      The limit of iteration number has been reached'',t95,''*** (mini  ) ***'')')
   call DEVIAT(obs,prf,weight,nwl,di,dev)     
   dev = sqrt(dev)
   exit ! From main loop
  endif          
 else
  exit ! From main loop
 endif 
enddo ! End of main loop     

call DCOPY(n,xopt,1,x,1)
call APPLY_NNORM(x,n,nonneg_flag)
call DCOPY(n,x,1,xopt,1)

if(sw == 0) write(6,"(1h+'      Correct SLD profiles:            OK          Mean DEV.=',F7.3,'%',/)") dev*100
if(sw == 1) write(6,"(1h+'     Correct INSTRUMENTAL profiles:   OK          Mean DEV.=',F7.3,'%',/)") dev*100
if(sw == 2) write(6,"(1h+'     Correct LINE profiles:           OK          Mean DEV.=',F7.3,'%',/)") dev*100 

deallocate (x, prf, ipiv, work, sigc, diag, rhsold, deriva, derivc, rhs,aderiv)
close(6)

end subroutine MINI_REG

subroutine REGFUNQ(x,absc,aderiv,regpar,n,regul)
real(8) x(n),absc(n),aderiv(n),regpar(2),regul,wt,da,df
 regul=0.d0
  do i=2,n
   i1=i
   i2=i-1
   wt=1.d0-regpar(2)*0.5d0*(aderiv(i1)+aderiv(i2))
   da=absc(i1)-absc(i2)
   df=x(i1)-x(i2)
   regul=regul+wt*(df/da)**2
  enddo 
  regul=regpar(1)*regul/n  
end subroutine REGFUNQ

subroutine REGDER(x,absc,aderiv,regpar,n,derivc,rhs)
real(8) x(n),absc(n),aderiv(n),regpar(2),rgn,wt,da,df,derivc(n,n),rhs(n)
  do i=2,n
   i1=i
   i2=i-1
   wt=1.d0-regpar(2)*0.5d0*(aderiv(i1)+aderiv(i2))
   da=absc(i1)-absc(i2)
   df=x(i1)-x(i2)
   rgn=regpar(1)*wt/n/(da*da)
   rhs(i1)=rhs(i1)-rgn*df
   rhs(i2)=rhs(i2)+rgn*df
   derivc(i1,i1)=derivc(i1,i1)+rgn
   derivc(i2,i2)=derivc(i2,i2)+rgn
   derivc(i1,i2)=derivc(i1,i2)-rgn
   derivc(i2,i1)=derivc(i2,i1)-rgn
  enddo
end subroutine REGDER

subroutine APPLY_NNORM(x,n, nonneg_flag)
implicit none
integer n, nonneg_flag, i
real(8) x(n)
  if(nonneg_flag /= 0) then
	do i = 1,n
		if(x(i) < 0.d0) x(i) = 0.d0
	enddo
  endif
end subroutine APPLY_NNORM


real(8) function VOIGT(nu, nu0, alpha, gamma)
 implicit none
 real(8), intent(in) :: nu, nu0, alpha, gamma
 real(8) :: wn, wn0, VOIGT_WN
 wn = 1.d8/nu
 wn0 = 1.d8/nu0
 VOIGT = VOIGT_WN(wn, wn0, alpha, gamma)
end function VOIGT

real(8) function VOIGT_WN(wn, wn0, alpha, gamma)
 implicit none
 real(8), intent(in) :: wn, wn0, alpha, gamma
 real(8) :: sigma, x, y
 real(8), parameter :: RT2LN2 = 1.1774100225154747d0 ! sqrt(2*ln2)
 real(8), parameter :: RT2 = 1.4142135623730951d0 ! sqrt(2)
 real(8), parameter :: RT2PI = 2.5066282746310002d0 ! sqrt(2*pi)
 real(8) :: HUMLIC
 sigma = alpha/RT2LN2
 x = (wn0 - wn) / (sigma*RT2)
 y = gamma / (sigma*RT2)
 VOIGT_WN = HUMLIC(x,y) / (sigma*RT2PI)
end function VOIGT_WN

real(8) function HUMLIC(x,y)
 implicit none
 real(8), intent(in) :: x, y
 complex(8) :: t, u, hum
 real(8) :: s, ax
 t = cmplx(y, -x, kind=8)
 ax = abs(x)
 s = ax + y
 if (s >= 15.d0) then
	hum = t*0.5641896d0/(0.5d0 + t*t)
	HUMLIC = real(hum,8)
	return
 endif
 if (s >= 5.5d0) then
	u = t*t
	hum = t*(1.410474d0 + u*0.5641896d0)/(0.75d0 + u*(3.d0 + u))
	HUMLIC = real(hum,8)
	return
 endif
 if (y >= 0.195d0*ax - 0.176d0) then
	hum = (16.4955d0 + t*(20.20933d0 + t*(11.96482d0 + t*(3.778987d0 + t*0.5642236d0)))) / &
	(16.4955d0 + t*(38.82363d0 + t*(39.27121d0 + t*(21.69274d0 + t*(6.699398d0 + t)))))
	HUMLIC = real(hum,8)
	return
 endif
 u = t*t
 hum = exp(u) - t*(36183.31d0 - u*(3321.9905d0 - u*(1540.787d0 - u*(219.0313d0 - u*(35.76683d0 - u*(1.320522d0 - u*0.56419d0)))))) / &
 (32066.6d0 - u*(24322.84d0 - u*(9022.228d0 - u*(2186.181d0 - u*(364.2191d0 - u*(61.57037d0 - u*(1.841439d0 - u)))))))
 HUMLIC = real(hum,8)
end function HUMLIC


subroutine GAUS_FOR_RESOLUTION
 use INSTRUMENT
 use OBSERVED
 implicit none

 real(8) :: z(2000)
 real(8) :: lam0, sigma_lam, delta, sigma_pix, invsig2, dd
 integer :: i
 lam0 = 0.5d0*(wlast + wfirst)
! sigma ??? exp(-(x/sigma)^2) ??? ???????? Resolution (FWHM = lam0/Resolution)
 sigma_lam = (lam0/Resolution) / (2.d0*sqrt(log(2.d0)))
 delta = wobs(2) - wobs(1)
 sigma_pix = sigma_lam / delta
 nip = 0
 invsig2 = 1.d0/(sigma_pix*sigma_pix)
 do i = 1, 2000
	dd = dble(i-1)*dble(i-1) * invsig2
	if(dd > 20.d0) exit
	nip = nip + 1
	z(nip) = exp(-dd)
 enddo
 if(allocated(inst_prof)) deallocate(inst_prof)
 allocate(inst_prof(nip))
 inst_prof = z(1:nip)
end subroutine GAUS_FOR_RESOLUTION


integer function LOCATE_LEFT(x, n, val)
 implicit none
 integer :: n, lo, hi, mid
 real(8) :: x(n), val
 lo = 1
 hi = n
 if(val <= x(1)) then
	LOCATE_LEFT = 1
	return
 endif
 if(val > x(n)) then
	LOCATE_LEFT = n+1
	return
 endif
 do while (lo < hi)
	mid = (lo + hi)/2
	if(x(mid) < val) then
		lo = mid + 1
	else
		hi = mid
	endif
 enddo
 LOCATE_LEFT = lo
end function LOCATE_LEFT

integer function LOCATE_RIGHT(x, n, val)
 implicit none
 integer :: n, lo, hi, mid
 real(8) :: x(n), val
 lo = 1
 hi = n
 if(val < x(1)) then
	LOCATE_RIGHT = 0
	return
 endif
 if(val >= x(n)) then
	LOCATE_RIGHT = n
	return
 endif
 do while (lo < hi)
	mid = (lo + hi + 1)/2
	if(x(mid) <= val) then
		lo = mid
	else
		hi = mid - 1
	endif
 enddo
 LOCATE_RIGHT = lo
end function LOCATE_RIGHT


subroutine START_TIMER(time_start)
 use DFPORT
 implicit none
 character(8), intent(out) :: time_start
 time_start = clock()
end subroutine START_TIMER

subroutine SHOW_DURATION(time_start)
 use DFPORT
 implicit none
 character(*), intent(in) :: time_start
 character(8) :: time_finish
 integer :: ih1, im1, is1, ih2, im2, is2
 integer :: t1, t2, dt
 integer :: hh, mm, ss

 time_finish = clock()

 read(time_start,'(i2,1x,i2,1x,i2)') ih1, im1, is1
 read(time_finish,'(i2,1x,i2,1x,i2)') ih2, im2, is2

 t1 = ih1*3600 + im1*60 + is1
 t2 = ih2*3600 + im2*60 + is2
 dt = t2 - t1
 if(dt < 0) dt = dt + 24*3600

 hh = dt / 3600
 mm = (dt - hh*3600) / 60
 ss = dt - hh*3600 - mm*60

 if(hh > 0) then
	write(*,'(/,a13,i3,a1,i2.2,a1,i2.2,/)') '  Duration:  ', hh, ':', mm, ':', ss
 else
	write(*,'(/,a13,i2.2,a1,i2.2,/)') '  Duration:  ', mm, ':', ss
 endif

end subroutine SHOW_DURATION

 subroutine QUICK_SORT(a,n)
  use OBSERVED
  implicit none
  integer    :: i, j, n, k(1)
  integer    ::  a(n)
  integer, allocatable :: c(:), d(:), b(:)
  
  allocate(b(n),d(n),c(n))
  d = a
  do i = 1,n
	k = maxloc(d)
	b(i) = k(1)
	c = -1.d0
	do j = 1,n
		if(j/=b(i))c(j)=d(j)
	enddo
	d = c
  enddo
  c = b
  do i = 1,n
	b(i) = c(n+1-i)
  enddo
  allocate(bi(n))
  do i = 1,n
	bi(i) = b((n+1)-i)
  enddo
  deallocate(d,c,b)
 end subroutine QUICK_SORT

