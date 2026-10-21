subroutine ILINE_CORRECT		
 use OBSERVED
 use LINES
 use MODELS
 implicit none
 integer(4)              :: i, j, persentold
 real(8)                 :: persent, FMIN, tolerant = 1.d-5, dumf, dev, d
 real(8), allocatable    :: weights(:)
 external                :: iCORR_FUN

 call INIT_PROF_CORR

 open(6,carriagecontrol = 'fortran')
 
 call QUICK_SORT(npoint,nlines)

 do j = 1, iterations
    persentold = 0
	do i = 1,nlines
		persent = (real(i)/real((nlines)*iterations))*1.d2 - (1.d2/j-1.d2)
		if(int(persent) > persentold)then
			write(6,'(1h+"      Correct iLINE profiles:        ",F8.3,1x,"%")')persent
			persentold = int(persent)
		endif
		if(corr_target /= line_family(i)) cycle
		line_active = bi(i)
		line_icorr(line_active) = FMIN(pmin,pmax,iCORR_FUN,tolerant,dumf,0)
	enddo
 enddo
 if(allocated(bi)) deallocate(bi)
 close(13, status='delete')
 rewind(7)

 call LINE_ABSORBTION('skip','skip')
 allocate(weights(nobs)); weights = 1.d0
 call deviat(robs,model,weights,nobs,d,dev)     
 dev = sqrt(dev); deallocate(weights)
 write(6,"(1h+'      Correct iLINE profiles:          OK          Mean DEV.=',F7.3,'%')") dev*100
 write(*,*)

end subroutine ILINE_CORRECT


real(8) function iCORR_FUN(x,n)
 use LINES
 use OBSERVED
 use CONSTANTS
 use MODELS
 use INSTRUMENT
! use TELLURICA
! use WATER_VAROP_CONT
 implicit none
 integer              :: n, i, j, e1, e2, p, ios
 real(8)              :: vw, ri, x, a
 real(8), allocatable :: tmodel(:)
 real(8) :: cor
 real(8) :: f1, f2
 real(8), allocatable :: work(:)
 integer :: k, m, k1, k2, no

 iCORR_FUN = 0.d0; a = 0.d0 
 line_icorr(line_active) = x

 e1 = line_left(line_active); e2 = line_right(line_active)
 
 if(e2 <= e1) then
	iCORR_FUN = 1.d-30
	return
 endif

 allocate(tmodel(e2-e1+1)); tmodel = 0.d0

 rewind(13)

 do i = 1, nlines
	if(line_right(i) < e1 .or. line_left(i) > e2) then
		read(13, iostat=ios)
		if(ios /= 0) stop 'iCORR_FUN: cannot skip line profile'
		cycle
	endif

	allocate(vp(npoint_new(i)), rp(npoint_new(i)))

	read(13, iostat=ios) vp, rp
	if(ios /= 0) stop 'iCORR_FUN: cannot read line profile'

	vp = vp + lrv(i)

	cor = 1.d0

	if(allocated(corr_id) .and. allocated(corr_fac)) then
		p = corr_id(i)
		if(p > 0) cor = corr_fac(p)
	else
		if(allocated(line_corr)) cor = line_corr(i)
	endif

	do j = e1, e2
		vw = ((wobs(j)-wave_center(i))/wave_center(i))*cvel

		if(vw < vp(1) .or. vw > vp(npoint_new(i))) cycle

		call qinter(vp, rp, npoint_new(i), vw, ri, 1)

		tmodel(j-e1+1) = tmodel(j-e1+1) + ri * cor * line_icorr(i)
	enddo

	deallocate(vp, rp)
 enddo

 tmodel = 1.d0 - exp(-tmodel)
! do i = e1,e2
!	tmodel(i-e1+1) = tmodel(i-e1+1) + water_cont_tau(i) + rayscat_cont_tau(i)
! enddo

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

 iCORR_FUN = sqrt(a/(e2-e1+1))
 deallocate(tmodel, work)

end function iCORR_FUN



subroutine INIT_PROF_CORR
 use LINES
 use SLD
 use OBSERVED
 use CONSTANTS
 implicit none

 integer(4) :: i, il, ios
 real(8), allocatable :: v_prof(:), r_prof(:)
 real(8) :: stp
 integer :: LOCATE_LEFT, LOCATE_RIGHT
 real(8) :: wmin, wmax, tmp
 logical :: file_exists

 if(allocated(npoint_new)) deallocate(npoint_new)
 allocate(npoint_new(nlines))

 if(.not. allocated(line_streng)) allocate(line_streng(nlines))
 if(.not. allocated(line_left)) allocate(line_left(nlines))
 if(.not. allocated(line_right)) allocate(line_right(nlines))

 rewind(7)

 do i = 1, nlines
	allocate(v_prof(npoint(i)), r_prof(npoint(i)))

	read(7, iostat=ios) v_prof, r_prof
	if(ios /= 0) stop 'INIT_PROF_CORR: cannot read profile'

	if(line_family(i) == 1) stp = step_tell
	if(line_family(i) == 2) stp = step_star

	npoint_new(i) = (SLD_npoint-1) + (int((v_prof(npoint(i))-v_prof(1))/stp)/2)*2 + 1

	deallocate(v_prof, r_prof)
 enddo

 inquire(file='prof_sld.tmp', exist=file_exists)

 if(file_exists) then
	open(13, file='prof_sld.tmp', status='old', iostat=ios)
	if(ios == 0) close(13, status='delete')
 endif

 open(13, file='prof_sld.tmp', form='unformatted', access='sequential', status='new', action='readwrite', iostat=ios)

 if(ios /= 0) stop 'INIT_PROF_CORR: cannot open prof_sld.tmp'

 rewind(7)

 do i = 1, nlines
	il = line_sld(i)

	allocate(v_prof(npoint(i)), r_prof(npoint(i)))

	read(7, iostat=ios) v_prof, r_prof
	if(ios /= 0) stop 'INIT_PROF_CORR: cannot reread profile'

	if(line_family(i) == 1) stp = step_tell
	if(line_family(i) == 2) stp = step_star

	call CONVOLVE(stp, npoint(i), v_prof, r_prof, SLD_npoint, SLD_intensity(il,1:SLD_npoint))

	write(13) vp(1:npoint_new(i)), rp(1:npoint_new(i))

	line_streng(i) = maxval(rp(1:npoint_new(i)))

	wmin = wave_center(i) * (1.d0 + vp(1)/cvel)
	wmax = wave_center(i) * (1.d0 + vp(npoint_new(i))/cvel)

	if(wmin > wmax) then
		tmp = wmin
		wmin = wmax
		wmax = tmp
	endif

	line_left(i) = LOCATE_LEFT(wobs, nobs, wmin)
	line_right(i) = LOCATE_RIGHT(wobs, nobs, wmax)

	if(line_left(i) < 1) line_left(i) = 1
	if(line_right(i) > nobs) line_right(i) = nobs

	deallocate(v_prof, r_prof, vp, rp)
 enddo

 endfile(13)
 rewind(13)

end subroutine INIT_PROF_CORR
