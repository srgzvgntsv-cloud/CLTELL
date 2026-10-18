subroutine FORM_TELLURIC
 use OBSERVED; use LINES; use SLD
 implicit none
 integer                    :: i, j, nX, k
 real(8), allocatable       :: x(:), weights(:), absc(:)
 external                   :: SLD_FUN

 nX = SLD_npoint*nsld ! number of unknowns

 allocate (weights(nobs),x(nX),absc(nX))
 weights = 1.d0; k = 0

 do i = 1,nsld
	do j = 1,SLD_npoint
		k = k+1
		x(k) = SLD_intensity(i,j)
		absc(k) = j
	enddo
 enddo

 call MINI_REG(x,absc,robs,weights,regpar,nX,nobs,SLD_FUN,0,1)

 deallocate(weights,x,absc)

end subroutine FORM_TELLURIC




subroutine SLD_FUN(x,f,Jacobian,n,nwl,mode) ! Calculate model data for solver mini
 use MODELS; use SLD; use LINES
 implicit none
 integer               :: n, nwl, i, j, k
 real(8)              :: x(n),f(nwl),Jacobian(nwl,n)
 character(4)         :: mode

 k = 0
 do i = 1,nsld
	do j = 1,SLD_npoint
		k = k+1
		SLD_intensity(i,j) = x(k)
	enddo
 enddo
 
 call LINE_ABSORBTION(mode,'fsld',Jacobian)

 f = model

end subroutine SLD_FUN