subroutine READ_HITRAN
 use OBSERVED
 use TELLURICA
 use MOLECULE_ISOTOP
 use LINES
 use CONSTANTS
 implicit none

 type(isotop) :: molecule
 integer :: ios
 integer :: im, i, k, n, kmax
 integer, allocatable :: index_line(:)
 integer, external :: INDWAV

 real(8) :: fr, ab, ga, gs, en, ns, da
 real(8) :: wav
 logical :: keep
 real(8) :: depth
 real(8), parameter :: thr = 1.d-30

 open(1, file='Hitran_database.par', status='old', iostat=ios)
 if(ios /= 0) stop 'file hitran_database.par not exist'

 kmax = 0
 do
	read(1,"(i3,f12.6,e10.3,10x,f5.4,f5.3,f10.4,f4.2,f8.6,92x)", iostat=ios) im, fr, ab, ga, gs, en, ns, da
	if(ios /= 0) exit
	kmax = kmax + 1
 enddo
 rewind(1)

 allocate(index_line(kmax))
 index_line = 0
 nlin_tell = 0

 do k = 1, kmax
	read(1,"(i3,f12.6,e10.3,10x,f5.4,f5.3,f10.4,f4.2,f8.6,92x)", iostat=ios) im, fr, ab, ga, gs, en, ns, da
	if(ios /= 0) exit
	molecule = get_id(im)
	if(molecule.id_molecule == 0) cycle
	if(fr <= 0.d0) cycle
	wav = 1.d8/fr
	if(wav < wfirst) cycle
	if(wav > wlast) cycle
	if(INDWAV(wav, 0.d0) == 0) cycle
	call HITRAN_DEPTH_FILTER(thr, molecule.id_molecule, molecule.id_isotop, fr, ab, ga, gs, en, ns, da, keep, depth)
	if(.not. keep) cycle
	nlin_tell = nlin_tell + 1
	index_line(k) = 1
 enddo

 rewind(1)

 if(allocated(imol)) deallocate(imol)
 if(allocated(wave_center_tell)) deallocate(wave_center_tell)
 if(allocated(abscoef)) deallocate(abscoef)
 if(allocated(gamma_air)) deallocate(gamma_air)
 if(allocated(gamma_self)) deallocate(gamma_self)
 if(allocated(energyi)) deallocate(energyi)
 if(allocated(n_air)) deallocate(n_air)
 if(allocated(d_air)) deallocate(d_air)
 if(allocated(line_species_tell)) deallocate(line_species_tell)
 if(allocated(line_rv_tell)) deallocate(line_rv_tell)
 if(allocated(index_wavelenght_tell)) deallocate(index_wavelenght_tell)

 allocate(imol(nlin_tell), wave_center_tell(nlin_tell), abscoef(nlin_tell), gamma_air(nlin_tell), gamma_self(nlin_tell), energyi(nlin_tell), n_air(nlin_tell), d_air(nlin_tell))
 allocate(line_species_tell(nlin_tell), line_rv_tell(nlin_tell))
 allocate(index_wavelenght_tell(nlin_tell))

 n = 0
 do k = 1, kmax
	read(1,"(i3,f12.6,e10.3,10x,f5.4,f5.3,f10.4,f4.2,f8.6,92x)", iostat=ios) im, fr, ab, ga, gs, en, ns, da
	if(ios /= 0) exit
	if(index_line(k) == 1) then
		n = n + 1
		wav = 1.d8/fr
		i = nlin_tell + 1 - n
		imol(i) = im
		wave_center_tell(i) = wav
		abscoef(i) = ab
		gamma_air(i) = ga
		gamma_self(i) = gs
		energyi(i) = en
		n_air(i) = ns
		d_air(i) = da
		index_wavelenght_tell(i) = INDWAV(wav, 0.d0)
	endif
 enddo

 deallocate(index_line)
 close(1)

 do i = 1, nlin_tell
	molecule = get_id(imol(i))
	line_species_tell(i) = molecule.name_molecule
	line_rv_tell(i) = 0.d0
 enddo

 contains

 subroutine HITRAN_DEPTH_FILTER(thr, id_molecule, id_isotop, nu0, Sref, gamma_air, gamma_self, elow, n_air, d_air, keep, depth)
  use TELLURICA
  use CONSTANTS
  implicit none

  real(8), intent(in) :: thr
  integer, intent(in) :: id_molecule, id_isotop
  real(8), intent(in) :: nu0, Sref, gamma_air, gamma_self, elow, n_air, d_air
  logical, intent(out) :: keep
  real(8), intent(out) :: depth

  integer :: j
  real(8) :: tau_peak
  real(8) :: nu_shift, Q_Q, E_E, T_T, ax
  real(8) :: dopw, lorw, vo0
  real(8) :: abloc, abloc2
  real(8) :: PARTF, VOIGT_WN

  if(id_molecule < 1 .or. id_molecule > 9) then
	keep = .false.
	depth = 0.d0
	return
  endif

  if(Sref <= 0.d0 .or. nu0 <= 0.d0) then
	keep = .false.
	depth = 0.d0
	return
  endif

  if(nlayer <= 0) stop 'HITRAN_DEPTH_FILTER: nlayer <= 0'
  if(.not. allocated(press)) stop 'HITRAN_DEPTH_FILTER: press not allocated'
  if(.not. allocated(teff)) stop 'HITRAN_DEPTH_FILTER: teff not allocated'
  if(.not. allocated(abudance_mol)) stop 'HITRAN_DEPTH_FILTER: abudance_mol not allocated'

  tau_peak = 0.d0

  do j = 1, nlayer

	abloc = abudance_mol(id_molecule, j)
	abloc2 = abudance(id_molecule, j)
	if(abloc <= 0.d0) cycle
	if(abloc2 <= 0.d0) cycle

	nu_shift = nu0 + d_air*(press(j)/Pref)

	Q_Q = PARTF(id_molecule, id_isotop, Tref) / PARTF(id_molecule, id_isotop, teff(j))
	E_E = exp(-cl2*elow/teff(j)) / exp(-cl2*elow/Tref)
	T_T = (1.d0 - exp(-cl2*nu_shift/teff(j))) / (1.d0 - exp(-cl2*nu_shift/Tref))
	ax = Q_Q * E_E * T_T

	dopw = (nu0/vellgt) * const1 * sqrt(teff(j)/meanmass(id_molecule,id_isotop))

	lorw = ((1.d0-abloc2)*gamma_air + abloc2*gamma_self) * (press(j)/Pref) * (Tref/teff(j))**n_air

	vo0 = VOIGT_WN(nu_shift, nu_shift, dopw, lorw)

	tau_peak = tau_peak + vo0 * (ax * Sref * abloc)

  enddo

  if(tau_peak <= 0.d0) then
	depth = 0.d0
	keep = .false.
	return
  endif

  depth = 1.d0 - exp(-tau_peak)
  keep = (depth >= thr)

 end subroutine HITRAN_DEPTH_FILTER

end subroutine READ_HITRAN



