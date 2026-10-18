subroutine FIND_RV
 use MODELS; use CONSTANTS; use OBSERVED; use LINES
 implicit none
 integer :: i
 integer, external :: INDWAV
 real(8) :: ccf_best
 real(8), parameter :: vmn = -200.d0
 real(8), parameter :: vmx = 200.d0
 real(8), parameter :: dv_coarse = 1.0d0
 real(8), parameter :: dv_fine = 0.05d0
 real(8), parameter :: vwin = 3.0d0

 allocate(lrv(nlines)); lrv = 0.d0

 call LINE_ABSORBTION('skip','tell')
 call FIND_RV_ONE(vmn, vmx, dv_coarse, dv_fine, vwin, RV_tell, ccf_best)
! RV_tell = 0.d0
 do i = 1, nlines
	if(line_family(i) == 1) line_rv(i) = RV_tell
 enddo

 call LINE_ABSORBTION('skip','star')
 call FIND_RV_ONE(vmn, vmx, dv_coarse, dv_fine, vwin, RV_star, ccf_best)
 RV_star = 0.d0
 do i = 1, nlines
	if(line_family(i) == 2) line_rv(i) = RV_star
 enddo

 do i = 1, nlines
	index_wavelenght(i) = INDWAV(wave_center(i), line_rv(i))
 enddo

 write(*,'("  RV of telluric lines is defined as: [",f7.2,"] km/s",/)') RV_tell
 write(*,'("  RV of star lines is defined as:     [",f7.2,"] km/s",/)') RV_star

end subroutine FIND_RV


subroutine FIND_RV_ONE(vminl, vmaxl, dv1, dv2, win, rv_out, ccf_out)
 implicit none
 real(8), intent(in) :: vminl, vmaxl, dv1, dv2, win
 real(8), intent(out) :: rv_out, ccf_out
 real(8) :: v, vbest, ccf, ccfmax
 real(8) :: vstart, vend
 real(8) :: y1, y2, y3, denom, off, rv2, ccf2
 real(8) :: CCF_AT
 
 ccfmax = -1.d30
 vbest = 0.d0
 v = vminl
 do while (v <= vmaxl + 1.d-12)
	ccf = CCF_AT(v)
	if(ccf > ccfmax) then
		ccfmax = ccf
		vbest = v
	endif
	v = v + dv1
 enddo
 ccfmax = -1.d30
 rv_out = vbest
 vstart = max(vminl, vbest - win)
 vend   = min(vmaxl, vbest + win)
 v = vstart
 do while (v <= vend + 1.d-12)
	ccf = CCF_AT(v)
	if(ccf > ccfmax) then
		ccfmax = ccf
		rv_out = v
	endif
	v = v + dv2
 enddo
 if(rv_out - dv2 >= vminl .and. rv_out + dv2 <= vmaxl) then
	y1 = CCF_AT(rv_out - dv2)
	y2 = CCF_AT(rv_out)
	y3 = CCF_AT(rv_out + dv2)
	denom = (y1 - 2.d0*y2 + y3)
	if(abs(denom) > 1.d-12) then
		off = 0.5d0*dv2*(y1 - y3)/denom
		if(off >  dv2) off =  dv2
		if(off < -dv2) off = -dv2
		rv2 = rv_out + off
		ccf2 = CCF_AT(rv2)
		if(ccf2 > y2) then
			rv_out = rv2
			ccfmax = ccf2
		else
			ccfmax = y2
		endif
	else
		ccfmax = y2
	endif
 endif
 ccf_out = ccfmax

end subroutine FIND_RV_ONE


real(8) function CCF_AT(v)
 use OBSERVED; use MODELS; use CONSTANTS
 implicit none
 real(8), intent(in) :: v
 integer :: k, m
 real(8) :: mean_t, mean_o, sx, sy, sxy
 real(8) :: xo, yt
 real(8), allocatable :: xint(:), yint(:)

 allocate(xint(nobs), yint(nobs))
 do k = 1, nobs
	xint(k) = wobs(k) * (1.d0 - v/cvel)
 enddo
 call QINTER(wobs, model, nobs, xint, yint, nobs)
 mean_t = 0.d0; mean_o = 0.d0; m = 0
 do k = 1, nobs
	if(xint(k) >= wobs(1) .and. xint(k) <= wobs(nobs)) then
		mean_t = mean_t + yint(k)
		mean_o = mean_o + robs(k)
		m = m + 1
	endif
 enddo
 if(m < 3) then
	CCF_AT = -1.d30
	return
 endif
 mean_t = mean_t / dble(m)
 mean_o = mean_o / dble(m)
 sx = 0.d0; sy = 0.d0; sxy = 0.d0
 do k = 1, nobs
	if(xint(k) >= wobs(1) .and. xint(k) <= wobs(nobs)) then
		yt = yint(k) - mean_t
		xo = robs(k) - mean_o
		sx = sx + xo*xo
		sy = sy + yt*yt
		sxy = sxy + xo*yt
	endif
 enddo
 if(sx <= 0.d0 .or. sy <= 0.d0) then
	CCF_AT = -1.d30
 else
	CCF_AT = sxy / (sqrt(sx)*sqrt(sy))
 endif
 deallocate(xint, yint)

end function CCF_AT


subroutine TELL_MODEL_RV 
 use MODELS; use OBSERVED; use TELLURICA; use LINES
 implicit none
 integer :: i, k, j
 real(8) :: tau, wnf, wn0
 real(8) :: vprof, VOIGT_WN
 real(8), allocatable :: xint(:), yint(:)
 real(8) :: ccf_best
 real(8), parameter :: vmn = -200.d0
 real(8), parameter :: vmx = 200.d0
 real(8), parameter :: dv_coarse = 1.0d0
 real(8), parameter :: dv_fine = 0.05d0
 real(8), parameter :: vwin = 3.0d0

 model = 0.d0
 do i = 1, nobs
	tau = 0.d0
	wnf = 1.d8 / wobs(i)
	do k = 1, nlin_tell
		wn0 = 1.d8 / wave_center(k)
		do j = 1, nlayer
			vprof = VOIGT_WN(wnf, wn0, dopcoef(k,j), lorcoef(k,j))
			tau = tau + abscoef_layer(k,j) * vprof
		enddo
	enddo
	model(i) = 1.d0 - exp(-tau)
 enddo

 allocate(xint(nobs), yint(nobs))

 call FIND_RV_ONE(vmn, vmx, dv_coarse, dv_fine, vwin, RV_tell, ccf_best)

 write(*,'("  RV of telluric lines is defined as: [",f7.2,"] km/s",/)') RV_tell

end subroutine TELL_MODEL_RV