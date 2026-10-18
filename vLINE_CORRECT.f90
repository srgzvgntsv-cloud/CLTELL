subroutine VLINE_CORRECT
 use OBSERVED
 use LINES
 use MODELS
 implicit none
 integer(4)              :: i, j, persentold
 real(8)                 :: persent, FMIN, tolerant = 1.d-5, dumf, dev, d
 real(8), allocatable    :: weights(:)
 external                :: vCORR_FUN

 call INIT_PROF_CORR

 open(6,carriagecontrol = 'fortran')
 
 call QUICK_SORT(npoint,nlines)

 do j = 1, iterations
    persentold = 0
	do i = 1,nlines
		persent = (real(i)/real((nlines)*iterations))*1.d2 - (1.d2/j-1.d2)
		if(int(persent) > persentold)then
			write(6,'(1h+"      Correct vLINE profiles:        ",F8.3,1x,"%")')persent
			persentold = int(persent)
		endif
		line_active = bi(i)
		lrv(line_active) = FMIN(vmin,vmax,vCORR_FUN,tolerant,dumf,0)
	enddo
 enddo
 if(allocated(bi)) deallocate(bi)

 close(13, status='delete')

 call LINE_ABSORBTION('skip','skip')
 allocate(weights(nobs)); weights = 1.d0
 call deviat(robs,model,weights,nobs,d,dev)     
 dev = sqrt(dev); deallocate(weights)
 write(6,"(1h+'      Correct vLINE profiles:          OK          Mean DEV.=',F7.3,'%')") dev*100
 write(*,*)

end subroutine VLINE_CORRECT


real(8) function vCORR_FUN(x,n)
 use LINES
 use OBSERVED
 use CONSTANTS
 use MODELS
 use INSTRUMENT
 implicit none
 integer              :: n, i, j, e1, e2, p, ios
 real(8)              :: vw, ri, x, a
 real(8), allocatable :: tmodel(:)
 real(8) :: cor
 real(8) :: f1, f2
 real(8), allocatable :: work(:)
 integer :: k, m, k1, k2, no

 vCORR_FUN = 0.d0; a = 0.d0 
 lrv(line_active) = x

 e1 = line_left(line_active); e2 = line_right(line_active)

 if(e2 <= e1) then
	vCORR_FUN = 1.d-30
	return
 endif

 allocate(tmodel(e2-e1+1)); tmodel = 0.d0

 rewind(13)

 do i = 1, nlines
	if(line_right(i) < e1 .or. line_left(i) > e2) then
		read(13, iostat=ios)
		if(ios /= 0) stop 'vCORR_FUN: cannot skip line profile'
		cycle
	endif

	allocate (vp(npoint_new(i)), rp(npoint_new(i)))

	read(13, iostat=ios) vp, rp
	if(ios /= 0) stop 'vCORR_FUN: cannot read line profile'

	vp = vp+lrv(i)

	cor = 1.d0
	if(allocated(corr_id) .and. allocated(corr_fac)) then
		p = corr_id(i)
		if(p > 0) cor = corr_fac(p)
	else
		if(allocated(line_corr)) cor = line_corr(i) 
	endif

	do j = e1, e2 ! Include given line profile in model spectrum
		vw = ((wobs(j)-wave_center(i))/wave_center(i))*cvel ! Coordinate of observed spectrum point in velocity space
		if(vw < vp(1) .or. vw > vp(npoint_new(i))) cycle
		call qinter(vp,rp,npoint_new(i),vw,ri,1) ! Interpolate profile to observed spectrum point
		tmodel(j-e1+1) = tmodel(j-e1+1) + ri * cor * line_icorr(i)        ! Model spectrum
	enddo
	deallocate(vp,rp)
 enddo

 tmodel = 1.d0 - exp(-tmodel)

 no = e2-e1+1
 allocate(work(no))
 do i = 1, no
	k1 = i - nip
	if(k1 < 1) k1 = 1
	k2 = i + nip
	if(k2 > no) k2 = no
	f1 = 0.d0; f2 = 0.d0
	do k = k1, k2
		m = abs(k - i) + 1
		if(m > nip) cycle
		f2 = f2 + inst_prof(m)
		f1 = f1 + tmodel(k) * inst_prof(m)
	enddo
	if(f2 <= 0.d0) then
		work(i) = 0.d0
	else
		work(i) = f1 / f2
	endif
 enddo
 tmodel(1:no) = work(1:no)


 do i = e1,e2
	a = a+(robs(i)-tmodel(i-e1+1))**2
 enddo

 vCORR_FUN = sqrt(a/(e2-e1+1))
 deallocate(tmodel)

end function vCORR_FUN

