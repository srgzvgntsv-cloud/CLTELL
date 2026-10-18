subroutine LINE_CORRECT(cm,ct)
 use OBSERVED; use LINES
 implicit none

 integer(4) :: nX, cm, ct
 real(8), allocatable :: x(:), weights(:)
 external :: CORR_FUN

 corr_mode = cm; corr_target = ct
 if(.not. allocated(corr_id)) stop 'LINE_CORRECT: corr_id not allocated. Call SETUP_LINE_CORRECTION first.'
 if(.not. allocated(corr_fac)) stop 'LINE_CORRECT: corr_fac not allocated. Call SETUP_LINE_CORRECTION first.'

 nX = ncor
 if(nX <= 0) stop 'LINE_CORRECT: ncor <= 0 (nothing to fit)'

 allocate(weights(nobs), x(nX))
 weights = 1.d0

 x(1:nX) = corr_fac(1:nX)

 call MINI(x, robs, weights, nX, nobs, CORR_FUN, 2, 0)

 corr_fac(1:nX) = x(1:nX)

 deallocate(weights, x)

end subroutine LINE_CORRECT



subroutine CORR_FUN(x, f, Jacobian, n, nwl, mode)
 use OBSERVED; use MODELS; use LINES
 implicit none

 integer(4), intent(in) :: n, nwl
 real(8), intent(in) :: x(n)
 real(8), intent(out) :: f(nwl)
 real(8), intent(inout) :: Jacobian(nwl,n)
 character(4), intent(in) :: mode

 if(.not. allocated(corr_fac)) stop 'CORR_FUN: corr_fac not allocated'
 if(n /= ncor) stop 'CORR_FUN: size(x) /= ncor. Check LINE_CORRECT and SETUP_LINE_CORRECTION.'

 corr_fac(1:n) = x(1:n)

 call LINE_ABSORBTION(mode, 'fcor', Jacobian)

 f = model

end subroutine CORR_FUN


subroutine SETUP_LINE_CORRECTION
 use LINES
 implicit none

 integer :: il, p, found
 logical :: sel
 integer, allocatable :: fam_tmp(:)
 character(4), allocatable :: lab_tmp(:)

 if(nlines <= 0) then
	ncor = 0
	if(allocated(corr_id)) deallocate(corr_id)
	if(allocated(corr_fac)) deallocate(corr_fac)
	if(allocated(corr_label)) deallocate(corr_label)
	if(allocated(corr_family)) deallocate(corr_family)
	return
 endif

 if(.not. allocated(line_family)) stop 'SETUP_LINE_CORRECTION: line_family not allocated'
 if(.not. allocated(line_species)) stop 'SETUP_LINE_CORRECTION: line_species not allocated'


 if(.not. allocated(npoint)) stop 'SETUP_LINE_CORRECTION: npoint not allocated'

 if(allocated(corr_id)) deallocate(corr_id)
 allocate(corr_id(nlines))
 corr_id = 0

 if(allocated(corr_fac)) deallocate(corr_fac)
 if(allocated(corr_label)) deallocate(corr_label)
 if(allocated(corr_family)) deallocate(corr_family)

 allocate(fam_tmp(nlines), lab_tmp(nlines))
 fam_tmp = 0
 lab_tmp = ' '

 ncor = 0

 do il = 1, nlines
    if(npoint(il) <= 1) cycle
	sel = .true.
	select case(corr_target)
		case(0)
			sel = .true.                     ! star + tell
		case(1)
			if(line_family(il) /= 1) sel = .false.   ! tell
		case(2)
			if(line_family(il) /= 2) sel = .false.   ! star
		case default
			stop 'SETUP_LINE_CORRECTION: unknown corr_target'
	end select
	if(.not. sel) cycle
	if(corr_mode == 0) then
		ncor = ncor + 1
		corr_id(il) = ncor
		fam_tmp(ncor) = line_family(il)
		lab_tmp(ncor) = line_species(il)
	else if(corr_mode == 1) then
		found = 0
		do p = 1, ncor
			if(fam_tmp(p) == line_family(il) .and. lab_tmp(p) == line_species(il)) then
				found = p
				exit
			endif
		enddo
		if(found == 0) then
			ncor = ncor + 1
			fam_tmp(ncor) = line_family(il)
			lab_tmp(ncor) = line_species(il)
			corr_id(il) = ncor
		else
			corr_id(il) = found
		endif
	else if(corr_mode > 2) then
		stop 'SETUP_LINE_CORRECTION: unknown corr_mode'
	endif
 enddo

 if(ncor > 0) then
	allocate(corr_fac(ncor), corr_label(ncor), corr_family(ncor))
	corr_fac = 1.d0
	corr_label(1:ncor) = lab_tmp(1:ncor)
	corr_family(1:ncor) = fam_tmp(1:ncor)
 else
	allocate(corr_fac(1), corr_label(1), corr_family(1))
	corr_fac(1) = 1.d0
	corr_label(1) = 'NONE'
	corr_family(1) = 0
 endif

 deallocate(fam_tmp, lab_tmp)

end subroutine SETUP_LINE_CORRECTION
