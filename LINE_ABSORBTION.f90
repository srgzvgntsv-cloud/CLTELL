subroutine LINE_ABSORBTION(mode,mode2,Jacobian)
 use LINES
 use CONSTANTS
 use SLD
 use MODELS
 use OBSERVED
 use INSTRUMENT
 use WATER_VAROP_CONT
 use TELLURICA
 implicit none

 integer :: j, i, k, l, il, u, ios
 integer :: i1, i2, npar, p
 integer :: LOCATE_LEFT, LOCATE_RIGHT
 real(8) :: vw, ri, ya, di, arg, stp
 real(8) :: wmin, wmax, tmp, invcw
 real(8) :: cor
 character(4), optional :: mode, mode2
 real(8), optional :: Jacobian(:,:)
 real(8), allocatable :: vpt(:), rpt(:)
 character(4) :: cha
 character(4) :: lmode, lmode2
 character(1) :: csld
 logical :: ok

 lmode = 'func'; lmode2 = 'skip'
 if(present(mode)) lmode = adjustl(mode)
 if(present(mode2)) lmode2 = adjustl(mode2)

 model = 0.d0; model_tmp = 0.d0

 if(lmode == 'grad') then
	if(.not. present(Jacobian)) stop 'LINE_ABSORBTION: Jacobian is not present in grad mode'
	Jacobian = 0.d0; if(lmode2 == 'fcor') Jacobian = 1.d-5
 endif
 
 rewind(7)
 
 do il = 1, nlines

 	if(lmode2 == 'tell' .and. line_family(il) == 2) then
		read(7, iostat=ios)
		if(ios /= 0) stop 'LINE_ABSORBTION: cannot skip profile'
		cycle
	endif

	if(lmode2 == 'star' .and. line_family(il) == 1) then
		read(7, iostat=ios)
		if(ios /= 0) stop 'LINE_ABSORBTION: cannot skip profile'
		cycle
	endif

	allocate(vpt(npoint(il)), rpt(npoint(il)))

	read(7, iostat=ios) vpt, rpt
	if(ios /= 0) stop 'LINE_ABSORBTION: cannot read line profile'

    if(line_family(il) == 1) stp = step_tell
	if(line_family(il) == 2) stp = step_star

	call CONVOLVE(stp, npoint(il), vpt, rpt, SLD_npoint, SLD_intensity(line_sld(il),1:SLD_npoint))


	vp = vp+lrv(il)

	invcw = cvel / wave_center(il)

	wmin = wave_center(il) * (1.d0 + vp(1)/cvel)
	wmax = wave_center(il) * (1.d0 + vp(np)/cvel)
	if(wmin > wmax) then
		tmp = wmin
		wmin = wmax
		wmax = tmp
	endif

	i1 = LOCATE_LEFT(wobs, nobs, wmin)
	i2 = LOCATE_RIGHT(wobs, nobs, wmax)
	if(i1 < 1) i1 = 1
	if(i2 > nobs) i2 = nobs

	cor = 1.d0
	if(allocated(corr_id) .and. allocated(corr_fac)) then
		p = corr_id(il)
		if(p > 0) cor = corr_fac(p)
	else
		if(allocated(line_corr)) cor = line_corr(il) 
	endif

	if(i2 >= i1) then
		do i = i1, i2
			vw = (wobs(i) - wave_center(il)) * invcw
			call QINTER(vp, rp, np, vw, ri, 1)
			model(i) = model(i) + ri * cor * line_icorr(il) 
			if(lmode == 'grad' .and. lmode2 == 'fsld') then
				u = (line_sld(il)-1)*SLD_npoint
				do l = 1, SLD_npoint
					arg = vw - SLD_velocity(line_sld(il),l)
					call LININT_UNIFORM(arg, vp0, dvp, np, r3p, ya, ok)
					if(ok) Jacobian(i,l+u) = Jacobian(i,l+u) + ya * cor * line_icorr(il)
				enddo
			endif
			if(lmode == 'grad' .and. lmode2 == 'fcor') then
				if(allocated(corr_id) .and. allocated(corr_fac)) then
					p = corr_id(il)
					if(p > 0) Jacobian(i,p) = Jacobian(i,p) + ri * line_icorr(il)
				else
					Jacobian(i,il) = Jacobian(i,il) + ri * line_icorr(il)
				endif
			endif
		enddo
	endif
	deallocate(vpt, rpt)

	if(allocated(vp)) deallocate(vp)
	if(allocated(rp)) deallocate(rp)
	if(allocated(r3p)) deallocate(r3p)

 enddo

 if(lmode == 'grad') then
	npar = size(Jacobian,2)
	do i = 1, nobs
		di = exp(-model(i))
		Jacobian(i,1:npar) = Jacobian(i,1:npar) * di
	enddo
 endif

! model = model + water_cont_tau + rayscat_cont_tau
do i = 1,nobs
	model(i) = 1.d0 - exp(-model(i))
	model_tmp(i) = model(i)
enddo

 k = index(observed_names,'.',back=.true.)

 call INSTRUMENTAL('skip')

! if(lmode == 'grad' .and. lmode2 == 'fsld') then
!    npar = size(Jacobian,2)
!    call INSTRUMENTAL_SLD_JACOBIAN(Jacobian,npar)
! endif

 if(lmode2 /= 'tell' .and. lmode2 /= 'star') then
	open(901,file=trim(adjustl(synth_path))//observed_names(1:k)//'result')
	do i = 1, nobs
		write(901,*) wobs(i), 1.d0-model(i)
	enddo
	close(901)
 endif

 if(lmode == 'tell' .and. lmode2 == 'tell') then
	open(901,file=trim(adjustl(synth_path))//observed_names(1:k)//'no_tell')
	do i = 1, nobs
		write(901,*) wobs(i), 1.d0-(robs(i)-model(i))
	enddo
	close(901)
 endif

 if(lmode2 /= 'tell' .and. lmode2 /= 'star') then
	do i = 1, nsld
		do j = 1, nlines
			if(line_sld(j) == i) then
				cha = line_species(j)
				write(csld,'(i1)')line_sld(j)
				exit
			endif
		enddo 
		open(110+i,file=trim(adjustl(synth_path))//observed_names(1:k-1)//'_'//cha//'#'//csld//'.sld')
		do j = 1, SLD_npoint
			write(110+i,*) SLD_velocity(i,j), 1.d0-SLD_intensity(i,j)
		enddo
		close(110+i)
	enddo
 endif



end subroutine LINE_ABSORBTION