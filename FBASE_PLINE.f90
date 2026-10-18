subroutine FBASE_PLINE
 use LINES
 use CONSTANTS
 use TELLURICA
 use OBSERVED
 implicit none

 integer :: i, il, ios, j, k, n1, n2, maxlen, l
 integer, allocatable :: i1(:), i2(:), lp(:)
 real(8), allocatable :: rbuf(:), vbuf(:)
 real(8) :: cw, invcw, wn0, wn, sumv, maxp
 real(8), external :: VOIGT_WN
 real(8), allocatable :: vpl(:), rpl(:)

 real(8), allocatable :: tmp3(:), tmp4(:), tmp5(:)
 integer, allocatable :: tmp1(:), tmp2(:)
 character(4), allocatable :: tmp7(:)

 integer :: current_star_record, target_star_record
 logical :: file_exists

 allocate(i1(nlines), i2(nlines))
 
 allocate(line_intensity(nlines)); line_intensity = 0.d0

 maxlen = 1
 do il = 1, nlines
	if(line_family(il) == 1) then
		cw = wave_center(il)
		wn0 = 1.d8 / cw
		n1 = 1
		n2 = nobs
		do j = index_wavelenght(il), 1, -1
			call LINE_ABSORBTION_VALUE_WN(j, wn0, il)
			if(if_stop == 1 .or. j == 1) then
				n1 = j+1
				exit
			endif
		enddo
		do j = index_wavelenght(il)+1, nobs
			call LINE_ABSORBTION_VALUE_WN(j, wn0, il)
			if(if_stop == 1 .or. j == nobs) then
				n2 = j-1
				exit
			endif
		enddo
		i1(il) = n1
		i2(il) = n2
		npoint(il) = n2 - n1 + 1
	endif
 enddo

 maxlen = maxval(npoint(1:nlines))
 allocate(vpl(maxlen), rpl(maxlen))

 inquire(file='prof_lines.tmp', exist=file_exists)

 if(file_exists) then
	open(7, file='prof_lines.tmp', status='old', iostat=ios)
	if(ios == 0) close(7, status='delete')
 endif

 open(7, file='prof_lines.tmp', form='unformatted', &
 access='sequential', status='new', action='readwrite', iostat=ios)

 if(ios /= 0) stop 'FBASE_PLINE: cannot open prof_lines.tmp'
 
 allocate(lp(nlines)); lp = 1; l = 0

 do il = 1, nlines
	if(line_family(il) == 1) then
		if(npoint(il) < 1) then
			lp(il) = 0
			cycle
		endif
		l = l + 1
		allocate(rbuf(npoint(il)), vbuf(npoint(il)))
		n1 = i1(il)
		n2 = i2(il)
		rbuf = 0.d0
		vbuf = 0.d0
		cw  = wave_center(il)
		wn0 = 1.d8 / cw
		invcw = cvel / cw
		k = 0
		do i = n1, n2
			k = k + 1
			wn = 1.d8 / wobs(i)
			sumv = 0.d0
			do j = 1, nlayer
				sumv = sumv + abscoef_layer(il,j) * VOIGT_WN(wn, wn0, dopcoef(il,j), lorcoef(il,j))
			enddo
			rbuf(k) = sumv
			vbuf(k) = (wobs(i) - cw) * invcw
		enddo
		line_intensity(il) = exp(-maxval(rbuf))
		write(7) vbuf, rbuf
		deallocate(rbuf, vbuf)
	endif
 enddo
 nlin_tell = l
 deallocate(i1, i2)
 
 current_star_record = 0
 k = 0

 do il = 1, nlines
	if(line_family(il) == 2) then
		k = k + 1
		if(npoint(il) < 1) then
			lp(il) = 0
			cycle
		endif
		target_star_record = nrec_star(k)
		if(target_star_record <= current_star_record) then
			stop 'FBASE_PLINE: invalid star profile record order'
		endif
		do while(current_star_record < target_star_record - 1)
			read(35, iostat=ios)
			if(ios /= 0) stop 'FBASE_PLINE: cannot skip star profile'
			current_star_record = current_star_record + 1
		enddo
		l = l + 1
		read(35, iostat=ios) vpl(1:npoint(il)), rpl(1:npoint(il))
		if(ios /= 0) stop 'FBASE_PLINE: cannot read star profile'
		current_star_record = current_star_record + 1
		write(7) vpl(1:npoint(il)), -log(max(rpl(1:npoint(il)), 1.d-300))
	endif
 enddo

 close(35)
 endfile(7)
 rewind(7)

 deallocate(vpl, rpl)
 
! --------------------------------------------------------------------

 allocate(tmp1(nlines),tmp2(nlines),tmp3(nlines))
 k = 0
 do i = 1, nlines
	if(lp(i) == 1) then
		k = k + 1
		tmp1(k) = line_molid(i)
		tmp2(k) = line_family(i)
		tmp3(k) = line_intensity(i)
	endif
 enddo
 deallocate(line_molid, line_family, line_intensity)
 allocate(line_molid(l), line_family(l), line_intensity(l))
 line_molid = tmp1(1:k); line_family = tmp2(1:k); line_intensity = tmp3(1:k)
 tmp1 = 0; tmp2 = 0; deallocate(tmp3)
 
 k = 0
 do i = 1, nlines
	if(lp(i) == 1) then
		k = k + 1
		tmp1(k) = index_wavelenght(i)
		tmp2(k) = npoint(i)
	endif
 enddo
 deallocate(index_wavelenght, npoint)
 allocate(index_wavelenght(l), npoint(l))
 index_wavelenght = tmp1(1:k); npoint = tmp2(1:k)
 deallocate(tmp1, tmp2)

 allocate(tmp4(nlines),tmp5(nlines))
 k = 0
 do i = 1, nlines
	if(lp(i) == 1) then
		k = k + 1
		tmp4(k) = line_corr(i)
		tmp5(k) = line_icorr(i)
	endif
 enddo
 deallocate(line_corr, line_icorr)
 allocate(line_corr(l), line_icorr(l))
 line_corr = tmp4(1:k); line_icorr = tmp5(1:k)
 tmp4 = 0; tmp5 = 0

 k = 0
 do i = 1, nlines
	if(lp(i) == 1) then
		k = k + 1
		tmp4(k) = line_rv(i)
		tmp5(k) = wave_center(i)
	endif
 enddo
 deallocate(line_rv, wave_center)
 allocate(line_rv(l), wave_center(l))
 line_rv = tmp4(1:k); wave_center = tmp5(1:k)
 deallocate(tmp4, tmp5)

 allocate(tmp7(nlines))
 k = 0
 do i = 1, nlines
	if(lp(i) == 1) then
		k = k + 1
		tmp7(k) = line_species(i)
	endif
 enddo
 deallocate(line_species)
 allocate(line_species(l))
 line_species = tmp7(1:k)
 deallocate(tmp7)

! --------------------------------------------------------------------------------------

 nlines = l
 nlin_star = nlines - nlin_tell

end subroutine FBASE_PLINE


subroutine LINE_ABSORBTION_VALUE_WN(iw, wn0, nl)
 use TELLURICA
 use OBSERVED
 implicit none

 integer, intent(in) :: iw, nl
 real(8), intent(in) :: wn0

 integer :: j
 real(8) :: maxv, tmp, wn
 real(8), parameter :: thr = 1.d-4
 real(8), external :: VOIGT_WN

 if_stop = 0
 maxv = -1.d300

 wn = 1.d8 / wobs(iw)

 do j = 1, nlayer
	tmp = abscoef_layer(nl,j) * VOIGT_WN(wn, wn0, dopcoef(nl,j), lorcoef(nl,j))
	if(tmp > maxv) maxv = tmp
 enddo

 if(maxv < thr) if_stop = 1

end subroutine LINE_ABSORBTION_VALUE_WN