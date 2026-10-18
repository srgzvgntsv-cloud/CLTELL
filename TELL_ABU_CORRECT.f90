module TELL_ABU_FIT_MOL
 use OBSERVED; use LINES; use TELLURICA; use CONSTANTS
 implicit none
 integer :: ntell = 0
 integer, allocatable :: tell_il(:) 
 integer, allocatable :: tell_molid(:) 
 real(8), allocatable :: tell_wn0(:)
 real(8), allocatable :: abs_base(:,:) ! (ntell,nlayer)
 real(8), allocatable :: dop_base(:,:) ! (ntell,nlayer)
 real(8), allocatable :: lor_base(:,:) ! (ntell,nlayer)
 integer :: nfit = 0
 real(8), allocatable :: wfit(:), rfit(:), wnfit(:), wgt(:)
 integer :: npar_mol = 0
 integer :: mol_parid(9) = 0 ! mol -> p (1..npar_mol)
 integer, allocatable :: par_molid(:) ! p -> mol (1..9)
 real(8), external :: VOIGT_WN, HUMLIC

 contains

 subroutine PREP_TELL_ABU_FIT_MOL(vwin_kms)
  implicit none
  real(8), intent(in) :: vwin_kms
  integer, external :: LOCATE_LEFT, LOCATE_RIGHT
  integer :: il, k, j, i1, i2, i, kk, m, p
  real(8) :: cw, cw_shift, wmin, wmax, tmp, rvfac
  logical, allocatable :: mask(:)
  logical :: mol_ok(9)

  if(nlines <= 0) stop 'PREP_TELL_ABU_FIT_MOL: nlines <= 0'
  if(nlayer <= 0) stop 'PREP_TELL_ABU_FIT_MOL: nlayer <= 0'
  if(nobs <= 0) stop 'PREP_TELL_ABU_FIT_MOL: nobs <= 0'
  if(.not. allocated(line_family)) stop 'PREP_TELL_ABU_FIT_MOL: line_family not allocated'
  if(.not. allocated(line_molid))  stop 'PREP_TELL_ABU_FIT_MOL: line_molid not allocated'
  if(.not. allocated(wave_center)) stop 'PREP_TELL_ABU_FIT_MOL: wave_center not allocated'
  if(.not. allocated(abscoef_layer)) stop 'PREP_TELL_ABU_FIT_MOL: abscoef_layer not allocated'
  if(.not. allocated(dopcoef)) stop 'PREP_TELL_ABU_FIT_MOL: dopcoef not allocated'
  if(.not. allocated(lorcoef)) stop 'PREP_TELL_ABU_FIT_MOL: lorcoef not allocated'

  ntell = 0
  do il = 1, nlines
	if(line_family(il) /= 1) cycle
	if(line_molid(il) < 1 .or. line_molid(il) > 9) cycle
	ntell = ntell + 1
  enddo
  if(ntell <= 0) stop 'PREP_TELL_ABU_FIT_MOL: no telluric lines found'
  if(allocated(tell_il)) deallocate(tell_il, tell_molid, tell_wn0)
  allocate(tell_il(ntell), tell_molid(ntell), tell_wn0(ntell))
  k = 0
  do il = 1, nlines
	if(line_family(il) /= 1) cycle
	if(line_molid(il) < 1 .or. line_molid(il) > 9) cycle
	k = k + 1
	tell_il(k) = il
	tell_molid(k) = line_molid(il)
	tell_wn0(k) = 1.d8 / wave_center(il)
  enddo

  if(allocated(abs_base)) deallocate(abs_base, dop_base, lor_base)
  allocate(abs_base(ntell,nlayer), dop_base(ntell,nlayer), lor_base(ntell,nlayer))
  do k = 1, ntell
	il = tell_il(k)
	do j = 1, nlayer
		abs_base(k,j) = abscoef_layer(il,j)
		dop_base(k,j) = dopcoef(il,j)
		lor_base(k,j) = lorcoef(il,j)
    enddo
  enddo

  allocate(mask(nobs))
  mask = .false.
  rvfac = 1.d0 + RV_tell/cvel
  do k = 1, ntell
	il = tell_il(k)
	cw = wave_center(il)
	cw_shift = cw * rvfac
	wmin = cw_shift * (1.d0 - vwin_kms/cvel)
	wmax = cw_shift * (1.d0 + vwin_kms/cvel)
	if(wmin > wmax) then
		tmp = wmin; wmin = wmax; wmax = tmp
	endif
	i1 = LOCATE_LEFT(wobs, nobs, wmin)
	i2 = LOCATE_RIGHT(wobs, nobs, wmax)
	if(i1 < 1) i1 = 1
	if(i2 > nobs) i2 = nobs
	if(i2 >= i1) mask(i1:i2) = .true.
  enddo
  nfit = 0
  do i = 1, nobs
	if(mask(i)) nfit = nfit + 1
  enddo
  if(nfit <= 20) stop 'PREP_TELL_ABU_FIT_MOL: nfit too small'
  if(allocated(wfit)) deallocate(wfit, rfit, wnfit, wgt)
  allocate(wfit(nfit), rfit(nfit), wnfit(nfit), wgt(nfit))
  kk = 0
  do i = 1, nobs
	if(.not. mask(i)) cycle
	kk = kk + 1
	wfit(kk)  = wobs(i)
	rfit(kk)  = robs(i)
	wnfit(kk) = 1.d8 / wobs(i)
	wgt(kk)   = 1.d0
  enddo
  deallocate(mask)

  mol_ok = .false.
  do k = 1, ntell
	mol_ok(tell_molid(k)) = .true.
  enddo
  mol_parid = 0
  npar_mol = 0
  do m = 1, 9
	if(mol_ok(m)) then
		npar_mol = npar_mol + 1
		mol_parid(m) = npar_mol
	endif
  enddo
  if(allocated(par_molid)) deallocate(par_molid)
  allocate(par_molid(npar_mol))
  do m = 1, 9
	p = mol_parid(m)
	if(p > 0) par_molid(p) = m
  enddo
  write(6,'("  PREP_TELL_ABU_FIT_MOL: ntell =",i4,", nfit =",i5,", RV_tell =",f8.2,", npar_mol =",i1,/)') ntell, nfit, RV_tell, npar_mol

 end subroutine PREP_TELL_ABU_FIT_MOL


 subroutine TELL_ABU_MOL_CORRECT(vwin_kms)
  implicit none
  real(8), intent(in) :: vwin_kms
  integer(4) :: nX, p, i, ios, nmol
  real(8), allocatable :: x(:), weights(:)
  real(8) :: ab

  if(abu_fit == .true.) then
	call PREP_TELL_ABU_FIT_MOL(vwin_kms)
	if(npar_mol <= 0) then
		write(6,'("  TELL_ABU_MOL_CORRECT: nothing to fit (npar_mol = 0)")')
		return
	endif
	nX = npar_mol
	allocate(x(nX), weights(nfit))
	weights = wgt
	do p = 1, nX
		if(scale_mol(par_molid(p)) <= 0.d0) scale_mol(par_molid(p)) = 1.d0
		x(p) = log(scale_mol(par_molid(p)))
	enddo
	call MINI(x, rfit, weights, nX, nfit, MOL_FUN, 3, 1)
	do p = 1, nX
		scale_mol(par_molid(p)) = exp(x(p))
	enddo
	deallocate(x, weights)
	write(6,'("  TELL_ABU_MOL_CORRECT: scale_mol(1..9) = ")')
	write(6,'(9(1x,1pe10.3))') scale_mol(1:9)
  endif

  call APPLY_MOL_SCALE_TO_ABUDANCE
 
 end subroutine TELL_ABU_MOL_CORRECT

 subroutine MOL_FUN(x, f, Jacobian, n, nwl, mode)
  implicit none
  integer(4), intent(in) :: n, nwl
  real(8), intent(in) :: x(n)
  real(8), intent(out) :: f(nwl)
  real(8), intent(inout) :: Jacobian(nwl,n)
  character(4), intent(in) :: mode
  integer :: p, m
  real(8) :: loc_scale(9)

  if(n /= npar_mol) stop 'MOL_FUN: n /= npar_mol'
  if(nwl /= nfit) stop 'MOL_FUN: nwl /= nfit'
  loc_scale = scale_mol
  do p = 1, npar_mol
	m = par_molid(p)
	loc_scale(m) = exp(x(p))
  enddo
  Jacobian = 0.d0
  call TELL_MODEL_AND_JAC_MOL(loc_scale, f, Jacobian)

 end subroutine MOL_FUN


 subroutine TELL_MODEL_AND_JAC_MOL(loc_scale, f, Jac)
  implicit none
  real(8), intent(in) :: loc_scale(9)
  real(8), intent(out) :: f(nfit)
  real(8), intent(inout) :: Jac(nfit,npar_mol)
  integer :: i, k, j, p, m
  real(8) :: tau, facI
  real(8) :: wn0s, vprof
  real(8) :: rvfac_inv
  real(8) :: dtau_p

  rvfac_inv = 1.d0 / (1.d0 + RV_tell/cvel)
  do i = 1, nfit
	tau = 0.d0
	do k = 1, ntell
		m = tell_molid(k)
		wn0s = tell_wn0(k) * rvfac_inv
		do j = 1, nlayer
			vprof = VOIGT_WN(wnfit(i), wn0s, dop_base(k,j), lor_base(k,j))
			tau = tau + abs_base(k,j) * loc_scale(m) * vprof
		enddo
	enddo
	f(i) = 1.d0 - exp(-tau)
	facI = exp(-tau)
	do p = 1, npar_mol
		m = par_molid(p)
		dtau_p = 0.d0
		do k = 1, ntell
			if(tell_molid(k) /= m) cycle
			wn0s = tell_wn0(k) * rvfac_inv
			do j = 1, nlayer
				vprof = VOIGT_WN(wnfit(i), wn0s, dop_base(k,j), lor_base(k,j))
				dtau_p = dtau_p + abs_base(k,j) * loc_scale(m) * vprof
			enddo
		enddo
		Jac(i,p) = facI * dtau_p
	enddo
  enddo

 end subroutine TELL_MODEL_AND_JAC_MOL

 subroutine APPLY_MOL_SCALE_TO_ABUDANCE
  implicit none
  integer :: m, j

  if(size(abudance_mol,1) < 9) stop 'APPLY_MOL_SCALE_TO_ABUDANCE: abudance_mol first dim < 9'
  open(55,file = trim(adjustl(synth_path))//'scale_abudance.dat')
  do m = 1, 9
	do j = 1, nlayer
		abudance_mol(m,j) = abudance_mol(m,j) * scale_mol(m)
	enddo
	write(55,*)m, scale_mol(m)
  enddo
  close(55)
  write(6,'("  Abudance_mol updated by scale_mol",/)')
 end subroutine APPLY_MOL_SCALE_TO_ABUDANCE

end module TELL_ABU_FIT_MOL


