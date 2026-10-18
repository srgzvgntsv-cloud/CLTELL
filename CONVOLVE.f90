subroutine CONVOLVE(stp,n1,v1,r1,n2,r2)
 use LINES; use SLD
 implicit none
 integer :: n1, n2, half, i, j
 real(8) :: v1(n1), r1(n1), r2(n2), stp

 np = (n2-1) + (int((v1(n1)-v1(1))/stp)/2)*2 + 1

 if(allocated(vp)) deallocate(vp)
 if(allocated(rp)) deallocate(rp)
 if(allocated(r3p)) deallocate(r3p)

 allocate(vp(np), rp(np), r3p(np))

 vp((np-1)/2 + 1) = 0.d0
 r3p = 0.d0

 do i = ((np-1)/2), 1, -1
	vp((np-1)/2 + 1 - i) = - stp * dble(i)
	vp((np-1)/2 + 1 + i) = stp * dble(i)
 enddo

 call LINTER(v1,r1,n1,vp,r3p,np)

 do i = 1, np
	if(vp(i) < v1(1) .or. vp(i) > v1(n1)) r3p(i) = 0.d0
	if(r3p(i) < 0.d0) r3p(i) = 0.d0
 enddo

 rp = 0.d0
 half = (n2-1)/2

 do i = half+1, np+half
	do j = 1, n2
		if(i >= j .and. (i-j+1) <= np) then
			rp(i-half) = rp(i-half) + r3p(i-j+1)*r2(j)
		endif
	enddo
 enddo

 dvp = stp
 vp0 = vp(1)

end subroutine CONVOLVE