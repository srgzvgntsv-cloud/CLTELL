subroutine FIND_INSTRUMENTAL
 use OBSERVED; use INSTRUMENT; use LINES 
 implicit none
 integer                    :: nX
 real(8), allocatable       :: x(:), weights(:)
 external                   :: INSTRUMENTAL_FUN

 nX = nip ! number of unknowns

 allocate (weights(nobs),x(nX))

 weights = 1.d0
 
 x = inst_prof

 call MINI(x,robs,weights,nX,nobs,INSTRUMENTAL_FUN,1,0)

 deallocate(weights,x)

 call LINE_ABSORBTION('skip','skip')

end subroutine FIND_INSTRUMENTAL


subroutine INSTRUMENTAL_FUN(x,f,Jacobian2,n,nwl,mode)
 use INSTRUMENT; use MODELS
 implicit none
 integer              :: n, nwl
 real(8)              :: x(n),f(nwl),Jacobian2(nwl,n)
 character(4)         :: mode

 inst_prof = x

 call INSTRUMENTAL(mode,Jacobian2)

 f = model

end subroutine INSTRUMENTAL_FUN


subroutine INSTRUMENTAL(mode,Jacobian2)
 use INSTRUMENT
 use OBSERVED
 use MODELS
 implicit none
 character(4), optional :: mode
 real(8), optional :: Jacobian2(:,:)

 real(8) :: f1, f2
 real(8), allocatable :: work(:)
 real(8), allocatable :: sm(:)
 integer, allocatable :: cnt(:)
 integer :: i, k, m, k1, k2
 logical :: do_grad

 do_grad = .false.
 if(present(mode)) then
	if(mode == 'grad') do_grad = .true.
 endif

 if(do_grad) then
	if(.not. present(Jacobian2)) stop 'INSTRUMENTAL: Jacobian2 not present in grad mode'
	Jacobian2 = 0.d0
 endif

 model = model_tmp

 allocate(work(nobs))
 if(do_grad) then
	allocate(sm(nip), cnt(nip))
 endif

 do i = 1, nobs
	k1 = i - nip
	if(k1 < 1) k1 = 1
	k2 = i + nip
	if(k2 > nobs) k2 = nobs
	f1 = 0.d0; f2 = 0.d0
	if(do_grad) then
		sm = 0.d0
		cnt = 0
	endif
	do k = k1, k2
		m = abs(k - i) + 1
		if(m > nip) cycle
		f2 = f2 + inst_prof(m)
		f1 = f1 + model(k) * inst_prof(m)
		if(do_grad) then
			sm(m) = sm(m) + model(k)
			cnt(m) = cnt(m) + 1
		endif
	enddo
	if(f2 <= 0.d0) then
		work(i) = 0.d0
		if(do_grad) then
			do m = 1, nip
				Jacobian2(i,m) = 0.d0
			enddo
		endif
	else
		work(i) = f1 / f2
		if(do_grad) then
			do m = 1, nip
				if(cnt(m) > 0) then
					Jacobian2(i,m) = (sm(m) - work(i)*dble(cnt(m))) / f2
				else
					Jacobian2(i,m) = 0.d0
				endif
			enddo
		endif
	endif
 enddo

 model(1:nobs) = work(1:nobs)

 deallocate(work)
 if(do_grad) then
	deallocate(sm, cnt)
 endif

end subroutine INSTRUMENTAL


subroutine INSTRUMENTAL_SLD_JACOBIAN(Jacobian2,npar)
 use INSTRUMENT
 use OBSERVED
 implicit none

 integer, intent(in) :: npar
 real(8), intent(inout) :: Jacobian2(nobs,npar)

 real(8), allocatable :: source(:), work(:)
 real(8) :: f1, f2
 integer :: i, k, m, k1, k2, p

 allocate(source(nobs),work(nobs))

 do p = 1, npar

    source(1:nobs) = Jacobian2(1:nobs,p)
    work = 0.d0

    do i = 1, nobs

       k1 = i-nip
       if(k1 < 1) k1 = 1

       k2 = i+nip
       if(k2 > nobs) k2 = nobs

       f1 = 0.d0
       f2 = 0.d0

       do k = k1, k2

          m = abs(k-i)+1

          if(m > nip) cycle

          f1 = f1+source(k)*inst_prof(m)
          f2 = f2+inst_prof(m)

       enddo

       if(f2 > 0.d0) then
          work(i) = f1/f2
       else
          work(i) = 0.d0
       endif

    enddo

    Jacobian2(1:nobs,p) = work(1:nobs)

 enddo

 deallocate(source,work)

end subroutine INSTRUMENTAL_SLD_JACOBIAN