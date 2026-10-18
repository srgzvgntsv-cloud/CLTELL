! Mathematical subroutines

pure subroutine MakeGaussMatrix(x1,x2,x,w,n)
implicit none
intent(in) :: x1,x2,n
intent(out) :: x,w
real(8), parameter :: pi=3.141525d0
real(8),parameter :: eps = 3.d-14 ! eps is the relative precision
integer(4) :: n,i,j,m
real(8)    :: x1,x2,x(n),w(n)
real(8)    :: p1,p2,p3,pp,xl,xm,z,z1

! The roots are symmetric in the interval, so we only have to find half of them
m = (n + 1)/2
xm = 0.5d0*(x2 + x1)
xl = 0.5d0*(x2 - x1)

! Loop over the desired roots
do i=1,m
  z = cos(pi*(i-0.25d0)/(n+0.5d0))
  ! Starting with the above approximation to the ith root, we enter the main loop of
  ! refinement by Newton's method
  do while( .true. )
    p1 = 1.d0; p2 = 0.d0
    ! Loop up the recurrence relation to get the legendre polynomial evaluated at z
    do j=1,n
      p3 = p2; p2 = p1; p1 = ((2.d0*j-1.d0)*z*p2-(j-1.d0)*p3)/j
    enddo
    ! p1 is now the desired legendre polynomial. We next compute pp, its derivative, by
    ! a standard relation involving also p2, the polynomial of one lower order
    pp = n*(z*p1-p2)/(z*z-1.d0)
    z1 = z; z = z1-p1/pp            ! Newton's method.
    if( abs(z-z1)<=eps ) exit
  enddo
  x(i) = xm - xl*z                  ! Scale the root to the desired interval,
  x(n+1-i) = xm + xl*z              ! and put in its symmetric counterpart.
  w(i) = 2.d0*xl/((1.d0-z*z)*pp*pp) ! Compute the weight
  w(n+1-i) = w(i)                   ! and its symmetric counterpart.
enddo

end subroutine MakeGaussMatrix
subroutine integ ( x, f, fint, n, start)
  real(8) x(*), f(*), fint(*)
  real(8) A(n), B(n), C(n)
  real(8) start
  call parcoe ( f, x, A, B, C, n )
  fint(1) = start
  do i = 1, n - 1
     fint(i+1) = fint(i) + (A(i) + B(i) / 2.d0 * (x(i+1) + x(i)) + &
                 C(i) / 3.d0 * ((x(i+1) + x(i)) * x(i+1) + x(i) *  &
				 x(i))) * (x(i+1) - x(i))
  enddo
end subroutine integ
subroutine parcoe ( f, x, A, B, C, n )
  real(8) f(*), x(*), A(*), B(*), C(*)
  real(8) wt, d
  C(1) = 0.d0
  B(1) = (f(2) - f(1)) / (x(2) - x(1))
  A(1) = f(1) - x(1) * B(1)
  C(n) = 0.d0
  B(n) = (f(n) - f(n-1)) / (x(n) - x(n-1))
  A(n) = f(n) - x(n) * B(n)
  if ( n > 2 ) then
	  do j = 2, n - 1
	     d    = (f(j) - f(j-1)) / (x(j) - x(j-1))
	     C(j) = f(j+1) / ((x(j+1) - x(j)) * (x(j+1) - x(j-1))) - f(j) / &
		        ((x(j) - x(j-1)) * (x(j+1) - x(j))) + f(j-1) /          &
				((x(j) - x(j-1)) * (x(j+1) - x(j-1)))
	     B(j) = d - (x(j) + x(j-1)) * C(j)
	     A(j) = f(j-1) - x(j-1) * d + x(j) * x(j-1) * C(j)
	  enddo
	  C(2) = 0.d0
	  B(2) = (f(3) - f(2)) / (x(3) - x(2))
	  A(2) = f(2) - x(2) * B(2)
	  C(3) = 0.d0
	  B(3) = (f(4) - f(3)) / (x(4) - x(3))
	  A(3) = f(3) - x(3) * B(3)
	  do j = 2, n - 1
	     if ( C(j) /= 0.d0 ) then
			wt   = abs(C(j+1)) / (abs(C(j+1)) + abs(C(j)))
			A(j) = A(j+1) + wt * (A(j) - A(j+1))
			B(j) = B(j+1) + wt * (B(j) - B(j+1))
			C(j) = C(j+1) + wt * (C(j) - C(j+1))
	     endif
	  enddo
	  A(n-1) = A(n)
	  B(n-1) = B(n)
	  C(n-1) = C(n)
  endif
end subroutine parcoe
subroutine linter ( xold, yold, nold, xnew, ynew, nnew )
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
end subroutine linter
integer function map1 ( xold, fold, nold, xnew, fnew, nnew )
 implicit real(8) (a-h,o-z)
 real(8) xold(*), fold(*), xnew(*), fnew(*)
 l  = 2
 ll = 0
 do k = 1, nnew
	do while ( xnew(k) >= xold(l) .and. l <= nold )
   	 l = l + 1
	enddo
	if ( l > nold ) then
		if ( l /= ll ) then
		
		   l  = amin0(nold,l)
		   c  = 0.d0
		   b  = (fold(l) - fold(l-1)) / (xold(l) - xold(l-1))
		   a  = fold(l) - xold(l) * b
		   ll = l
		
		endif
		
		fnew(k) = a + (b + c * xnew(k)) * xnew(k)
		cycle
	endif
	if ( l /= ll ) then
	    if ( l /= 2 .and. l /= 3 ) then
		   l1 = l - 1
		   if ( .not. ( l > ll+1 .or. ( l == 3 .or. l == 4 ) ) ) then
		
		      cbac = cfor
		      bbac = bfor
		      abac = afor
		    if ( l == nold ) then
		
			   c  = cbac
			   b  = bbac
			   a  = abac
			   ll = l
			   fnew(k) = a + (b + c * xnew(k)) * xnew(k)
			   cycle
		
		    endif
		else
			
		   l2 = l - 2
		   d  = (fold(l1) - fold(l2)) / (xold(l1) - xold(l2))
		   cbac = fold(l) / ((xold(l) - xold(l1)) * (xold(l) - xold(l2))) + &
				  (fold(l2) / (xold(l) - xold(l2)) - fold(l1) /             &
				  (xold(l) - xold(l1))) / (xold(l1) - xold(l2))
		   bbac = d - (xold(l1) + xold(l2)) * cbac
		   abac = fold(l2) - xold(l2) * d + xold(l1) * xold(l2) * cbac
		   if ( l >= nold ) then
		
			  c  = cbac
			  b  = bbac
			  a  = abac
			  ll = l
			  fnew(k) = a + (b + c * xnew(k)) * xnew(k)
			  cycle
		
		   endif
		endif
		d = (fold(l) - fold(l1)) / (xold(l) - xold(l1))
		cfor = fold(l+1) / ((xold(l+1) - xold(l)) * (xold(l+1) - xold(l1))) + &
			   (fold(l1) / (xold(l+1) - xold(l1)) - fold(l) /                 &
			   (xold(l+1) - xold(l))) / (xold(l) - xold(l1))
		bfor = d - (xold(l) + xold(l1)) * cfor
		afor = fold(l1) - xold(l1) * d + xold(l) * xold(l1) * cfor
		wt = 0.d0
		if ( abs(cfor) /= 0.d0 ) wt = abs(cfor) / ( abs(cfor) + abs(cbac) )
		
		   a  = afor + wt * (abac - afor)
		   b  = bfor + wt * (bbac - bfor)
		   c  = cfor + wt * (cbac - cfor)
		   ll = l
		   fnew(k) = a + (b + c * xnew(k)) * xnew(k)
		   cycle
		endif
		if ( l /= ll ) then
		
		   l  = amin0(nold,l)
		   c  = 0.d0
		   b  = (fold(l) - fold(l-1)) / (xold(l) - xold(l-1))
		   a  = fold(l) - xold(l) * b
		   ll = l
		
		endif
	endif
	fnew(k) = a + (b + c * xnew(k)) * xnew(k)
 enddo
 map1    = ll - 1
end function map1
subroutine deriv ( x, f, dfdx, n )
! ASSUMES THAT ANY ZERO IN X OCCURS AT A ENDPOINT
  real(8) x(*), f(*), dfdx(*)
  real(8) s, scale, d1, d, Tan1, Tan
  dfdx(1) = (f(2) - F(1)) / (x(2) - x(1))
  dfdx(n) = (f(n) - f(n-1)) / (x(n) - x(n-1))
  if ( n > 2) then
	  s = abs( x(2) - x(1)) / (x(2) - x(1) )
	  do j = 2, n - 1
	     scale = dmax1(abs(f(j-1)),abs(f(j)),abs(f(j+1))) / abs(x(j))
	     if ( scale == 0.d0 ) scale = 1.d0
	     d1      = (f(j+1) - f(j)) / (x(j+1) - x(j)) / scale
	     d       = (f(j) - f(j-1)) / (x(j) - x(j-1)) / scale
	     Tan1    = d1 / (s * sqrt(1.d0 + d1**2) + 1.d0)
	     Tan     = d / (s * sqrt(1.d0 + d**2) + 1.d0)
	     dfdx(j) = (Tan1 + Tan) / (1.d0 - Tan1 * Tan) * scale
	  enddo
  endif
end subroutine deriv
real(8) function Voigt(v,a) ! Width5 subroutine, modified by T.Kipper
 implicit real(8) (a-h,o-z)
 dimension H0(41),H1(81),H2(41),Gaus20(20),Gaus10(10),wht10(10),&
           Gaus3(3),wht3(3),wht20(20)
 Data H0/ &
  1.0000000d0,0.9900500d0,0.9607890d0,0.9139310d0,0.8521440d0,0.7788010d0,&
  0.6976760d0,0.6126260d0,0.5272920d0,0.4448580d0,0.3678790d0,0.2981970d0,&
  0.2369280d0,0.1845200d0,0.1408580d0,0.1053990d0,0.0773050d0,0.0555760d0,&
  0.0391640d0,0.0270520d0,0.0183156d0,0.0121552d0,0.0079071d0,0.0050418d0,&
  0.0031511d0,0.0019305d0,0.0011592d0,0.0006823d0,0.0003937d0,0.0002226d0,&
  0.0001234d0,0.0000671d0,0.0000357d0,0.0000186d0,0.0000095d0,0.0000048d0,&
  0.0000024d0,0.0000011d0,0.0000005d0,0.0000002d0,0.0000001d0/
 Data H1/ &
  -1.1283800d0,-1.1059600d0,-1.0404800d0,-0.9370300d0,-0.8034600d0,&
  -0.6494500d0,-0.4855200d0,-0.3219200d0,-0.1677200d0,-0.0301200d0,&
   0.0859400d0, 0.1778900d0, 0.2453700d0, 0.2898100d0, 0.3139400d0,&
   0.3213000d0, 0.3157300d0, 0.3009400d0, 0.2802700d0, 0.2564800d0,&
   0.2317260d0, 0.2075280d0, 0.1848820d0, 0.1643410d0, 0.1461280d0,&
   0.1302360d0, 0.1165150d0, 0.1047390d0, 0.0946530d0, 0.0860050d0,&
   0.0785650d0, 0.0721290d0, 0.0665260d0, 0.0616150d0, 0.0572810d0,&
   0.0534300d0, 0.0499880d0, 0.0468940d0, 0.0440980d0, 0.0415610d0,&
   0.0392500d0, 0.0351950d0, 0.0317620d0, 0.0288240d0, 0.0262880d0,&
   0.0240810d0, 0.0221460d0, 0.0204410d0, 0.0189290d0, 0.0175820d0,&
   0.0163750d0, 0.0152910d0, 0.0143120d0, 0.0134260d0, 0.0126200d0,&
   0.0118860d0, 0.0112145d0, 0.0105990d0, 0.0100332d0, 0.0095119d0,&
   0.0090306d0, 0.0085852d0, 0.0081722d0, 0.0077885d0, 0.0074314d0,&
   0.0070985d0, 0.0067875d0, 0.0064967d0, 0.0062243d0, 0.0059688d0,&
   0.0057287d0, 0.0055030d0, 0.0052903d0, 0.0050898d0, 0.0049006d0,&
   0.0047217d0, 0.0045526d0, 0.0043924d0, 0.0042405d0, 0.0040964d0,&
   0.0039595d0/
  Data H2/ &
    1.0000000d0, 0.9702000d0, 0.8839000d0, 0.7494000d0, 0.5795000d0,&
    0.3894000d0, 0.1953000d0, 0.0123000d0,-0.1476000d0,-0.2758000d0,&
   -0.3679000d0,-0.4234000d0,-0.4454000d0,-0.4392000d0,-0.4113000d0,&
   -0.3689000d0,-0.3185000d0,-0.2657000d0,-0.2146000d0,-0.1683000d0,&
   -0.1282100d0,-0.0950500d0,-0.0686300d0,-0.0483000d0,-0.0331500d0,&
   -0.0222000d0,-0.0145100d0,-0.0092700d0,-0.0057800d0,-0.0035200d0,&
   -0.0021000d0,-0.0012200d0,-0.0007000d0,-0.0003900d0,-0.0002100d0,&
   -0.0001100d0,-0.0000600d0,-0.0000300d0,-0.0000100d0,-0.0000100d0,&
    0.0000000d0/
  Data Gaus20/ &
    0.05d0,0.15d0,0.25d0,0.35d0,0.45d0,0.6d0,0.75d0,0.9d0,1.05d0,1.2d0,&
    1.35d0,1.5d0,1.65d0,1.8d0,1.95d0,2.1d0,2.25d0,2.4d0,2.55d0,2.7d0/
  Data  wht20/ &
    0.0996d0,0.097634d0,0.0938184d0,0.0883723d0,0.100389d0,0.104549d0,&
    0.0855171d0,0.0668929d0,0.0500384d0,0.035796d0,0.0244891d0,0.0160216d0,&
    0.010024d0,0.00599769d0,0.00343183d0,0.00187789d0,0.000982679d0,&
    0.000491757d0,0.000235342d0,0.0000681124d0/
  Data Gaus10/ &
    0.245341d0,0.737474d0,1.234076d0,1.738538d0,2.254974d0,2.788806d0,&
    3.347855d0,3.944764d0,4.603682d0,5.38748d0/
  Data wht10/  &
    0.462244d0,0.2866755d0,0.1090172d0,0.02481052d0,0.003243773d0,&
    0.0002283386d0,0.7802556d-5,0.1086069d-6,0.4399341d-9,0.2229394d-12/
  Data Gaus3/ 0.436077d0,1.335849d0,2.35605d0/
  Data wht3/0.7246296d0,0.1570673d0,0.00453001d0/
 if (a <= 0.175)then
  v0=v*10.d0
  n=v0
   if(n >= 120)then
     vv=v*v
     Voigt=(0.56419+0.846/vv)/vv*a
    return
   endif
   if(n < 40)then
     v1=n
     n=n+1
     v2=v0-v1
     n1=n+1
     Voigt=v2*(H0(n1)-H0(n)+a*(H1(n1)-H1(n)+a*(H2(n1)-H2(n))))+&
           H0(n)+a*(H1(n)+a*H2(n))
    return
   endif
     n=n/2+20
     v1=(n-20)*2
     n=n+1
     v2=(v0-v1)/2.d0
     n1=n+1
     Voigt=A*((H1(n1)-H1(n))*v2+H1(n))
    return
 endif
 if(a <= 0.5d0.and.v <= 2.7d0)then ! Super 20 point quadrature
  Voigt=0.d0
   do i=1,20
    Voigt=Voigt+wht20(i)*A/3.1415925d0*(1.d0/(A*A+(v-Gaus20(i))**2)+&
                1.d0/(a*a+(v+Gaus20(i))**2))
   enddo
  return
 endif
 if((a <= 0.5d0.and.v <= 7.d0).or.(a <= 1.d0.and.v <= 2.7d0))then ! 10 point
  Voigt=0.d0
   do i=1,10
    Voigt=Voigt+wht10(i)*a/3.1415925d0*(1.d0/(a*a+(v-Gaus10(i))**2)+&
                1.d0/(A*A+(v+Gaus10(i))**2))
   enddo
  return
 endif
 if((a <= 0.5d0.and.v <=100.d0).or.(a <= 1.d0.and.v <= 40.d0).or.&
    (a <= 15.d0.and.v <= 20.d0))then ! 3 point Gaussian quadrature
  Voigt=0.d0
   do i=1,3
    Voigt=Voigt+wht3(i)*A/3.1415925d0*(1.d0/(a*a+(v-Gaus3(i))**2)+&
                1.d0/(A*A+(v+Gaus3(i))**2))
   enddo
  return
 endif
   Voigt=a/1.77245d0/(a*a+v*v)  ! Lorentzian
end

real(8) function expi ( n, x )

 implicit none

! EXPONENTIAL INTEGRAL FOR POSITIVE ARGUMENTS AFTER CODY AND
! THACHER, MATH. OF COMP.,22,641(1968)
  real(8), parameter :: A0 = -44178.5471728217d0
  real(8), parameter :: A1 =  57721.7247139444d0
  real(8), parameter :: A2 =  9938.31388962037d0
  real(8), parameter :: A3 =  1842.11088668000d0
  real(8), parameter :: A4 =  101.093806161906d0
  real(8), parameter :: A5 =  5.03416184097568d0
  real(8), parameter :: B0 =  76537.3323337614d0
  real(8), parameter :: B1 =  32597.1881290275d0
  real(8), parameter :: B2 =  6106.10794245759d0
  real(8), parameter :: B3 =  635.419418378382d0
  real(8), parameter :: B4 =  37.2298352833327d0
  real(8), parameter :: C0 =  4.65627107975096d-7
  real(8), parameter :: C1 =  0.999979577051595d0
  real(8), parameter :: C2 =  9.04161556946329d0
  real(8), parameter :: C3 =  24.3784088791317d0
  real(8), parameter :: C4 =  23.0192559391333d0
  real(8), parameter :: C5 =  6.90522522784444d0
  real(8), parameter :: C6 =  0.430967839469389d0
  real(8), parameter :: D1 =  10.0411643829054d0
  real(8), parameter :: D2 =  32.4264210695138d0
  real(8), parameter :: D3 =  41.2807841891424d0
  real(8), parameter :: D4 =  20.4494785013794d0
  real(8), parameter :: D5 =  3.31909213593302d0
  real(8), parameter :: D6 =  0.103400130404874d0
  real(8), parameter :: E0 = -0.999999999998447d0
  real(8), parameter :: E1 = -26.6271060431811d0
  real(8), parameter :: E2 = -241.055827097015d0
  real(8), parameter :: E3 = -895.927957772937d0
  real(8), parameter :: E4 = -1298.85688746484d0
  real(8), parameter :: E5 = -545.374158883133d0
  real(8), parameter :: E6 = -5.66575206533869d0
  real(8), parameter :: F1 = 28.6271060422192d0
  real(8), parameter :: F2 = 292.310039388533d0
  real(8), parameter :: F3 = 1332.78537748257d0
  real(8), parameter :: F4 = 2777.61949509163d0
  real(8), parameter :: F5 = 2404.01713225909d0
  real(8), parameter :: F6 = 631.657483280800d0
  real(8) :: x1 = -1.0d20
  real(8)    x, ex, ex1
  integer n, i


  if ( x /= x1 ) then

    ex = exp( - x )
    x1 = x
    if ( x > 4.d0 ) then
      ex1 = (ex + ex * (E0 + (E1 + (E2 + (E3 + (E4 + (E5 + E6 / x) / x) / x) /&
             x) / x) / x) / &
             (x + F1 + (F2 + (F3 + (F4 + (F5 + F6 / x) / x) / x) / x) / x)) / x
    else if ( x > 1.d0 ) then
      ex1 = ex * (C6 + (C5 + (C4 + (C3 + (C2 + (C1 + C0 * x) * x) * x) * x) * &
            x) * x) / (D6 + (D5 + (D4 + (D3 + (D2 + (D1 + x) * x) * x) * x) * &
            x) * x)
    else if ( x > 0.d0 ) then
      ex1 = (A0 + (A1 + (A2 + (A3 + (A4 + A5 * x) * x) * x) * x) * x) / &
            (B0 + (B1 + (B2 + (B3 + (B4 + x) * x) * x) * x) * x) - log(x)
    else
      ex1 = 0.d0
    endif

  endif

  expi = ex1

  if ( n > 1 ) then

    do i = 1, n - 1
      expi = (ex - x * expi) / dble( i )	
    enddo

  endif


end function expi

subroutine solvit(A,n,b,ipivot)
 implicit real*8 (a-h,o-z)
!     solves linear equations
!     a is a completely filled n by n array that is destroyed
!     b is the right side vector of length n and returns as the solution
!     ipivot is a scratch area of length n
 dimension A(n,n),b(n),ipivot(n)
! equivalence (t,pivot,c)
  n1=n-1
   do i=1,n1
    m=i
    i1=i+1
    do  k=i1,n
     if(abs(A(k,i)) > abs(A(m,i)))m=k
    enddo
   ipivot(i)=m
   if(m /= i)then
    do k=i1,n
     t=A(i,k)
     A(i,k)=A(m,k)
     A(m,k)=t
    enddo
   endif
   pivot=1.d0/A(m,i)
   A(m,i)=A(i,i)
   a(i,i)=pivot
    do k=i1,n
     a(k,i)=a(k,i)*pivot
    enddo
    do j=i1,n
     c=A(i,j)
     if(c /= 0.d0)then
      do  k=i1,n
       A(k,j)=A(k,j)-A(k,i)*c
      enddo
     endif
    enddo
  enddo
   a(n,n)=1.d0/a(n,n)
  do i=1,n1
   m=ipivot(i)
   if(m /= i)then
    t=b(m)
    b(m)=b(i)
    b(i)=t
   endif
   c=b(i)
   i1=i+1
    do k=i1,n
     b(k)=b(k)-A(k,i)*c
    enddo
  enddo
  j1=n
  do i=1,n1
   j=j1
   j1=j1-1
   b(j)=b(j)*A(j,j)
   c=b(j)
    do k=1,j1
     b(k)=b(k)-A(k,j)*c
    enddo
  enddo
   b(1)=b(1)*A(1,1)
end 