subroutine GET_ABSCOEF(wr)
 use TELLURICA
 use PART_FUNCTION
 use CONSTANTS
 use MOLECULE_ISOTOP
 use LINES
 implicit none

 type(isotop) :: mol
 integer :: i, j, wr
 real(8) :: wave_cm, Qref, Qj, PARTF
 real(8) :: invTref, invTj, dInvT
 real(8) :: denomT, numT
 real(8) :: expE, expTerm
 real(8) :: a0, E0
 real(8), parameter :: tiny = 1.d-300

 if(allocated(abscoef_layer)) deallocate(abscoef_layer)
 allocate(abscoef_layer(nlin_tell,nlayer)); abscoef_layer = 0.d0

 invTref = 1.d0 / Tref

 do i = 1, nlin_tell
	wave_cm = 1.d8 / wave_center_tell(i) 
	mol = get_id(imol(i))
	Qref = PARTF(mol.id_molecule, mol.id_isotop, Tref)
	denomT = 1.d0 - exp(-cl2*wave_cm*invTref)
	if(abs(denomT) < tiny) denomT = tiny
	a0 = abscoef(i)
	E0 = energyi(i)
	do j = 1, nlayer
		invTj = 1.d0 / teff(j)
		dInvT = invTj - invTref
		Qj = PARTF(mol.id_molecule, mol.id_isotop, teff(j))
		if(Qj < tiny) Qj = tiny
		expE = exp(-cl2*E0*dInvT)
		numT = 1.d0 - exp(-cl2*wave_cm*invTj)
		expTerm = (Qref/Qj) * expE * (numT/denomT)
		abscoef_layer(i,j) = a0 * expTerm * abudance_mol(mol.id_molecule, j)
	enddo	
 enddo

 if(wr == 1) write(6,'(1h+" Preparation of Absorption coefficient by layers completed successfully                    ")')

end subroutine GET_ABSCOEF


subroutine GET_DOPCOEF(wr)
 use TELLURICA
 use CONSTANTS
 use MOLECULE_ISOTOP
 use LINES
 implicit none

 type(isotop) :: mol
 integer :: i, j, wr
 real(8) :: wave_cm, factor
 real(8), allocatable :: sqrtT(:)

 if(allocated(dopcoef)) deallocate(dopcoef)
 allocate(dopcoef(nlin_tell,nlayer)); dopcoef = 0.d0

 allocate(sqrtT(nlayer))
 do j = 1, nlayer
	sqrtT(j) = sqrt(teff(j))
 enddo

 do i = 1, nlin_tell
	wave_cm = 1.d8 / wave_center_tell(i)
	mol = get_id(imol(i))
	factor = (wave_cm/vellgt) * const1 / sqrt(meanmass(mol.id_molecule, mol.id_isotop))
	do j = 1, nlayer
		dopcoef(i,j) = factor * sqrtT(j)
	enddo
 enddo
 deallocate(sqrtT)

 if(wr == 1) write(6,'(1h+,a)')' Preparation of Doppler    coefficient by layers completed successfully                     '

end subroutine GET_DOPCOEF


subroutine GET_LORCOEF(wr)
 use TELLURICA
 use CONSTANTS
 use MOLECULE_ISOTOP
 use LINES
 implicit none

 type(isotop) :: mol
 integer :: i, j, wr
 real(8) :: pRatio, logTr, powTr
 real(8) :: ab, gmix
 real(8) :: nai
 real(8), allocatable :: pr(:), logTrefOverT(:)

 if(allocated(lorcoef)) deallocate(lorcoef)
 allocate(lorcoef(nlin_tell,nlayer)); lorcoef = 0.d0

 allocate(pr(nlayer), logTrefOverT(nlayer))
 do j = 1, nlayer
	pr(j) = press(j) / Pref
	logTrefOverT(j) = log(Tref / teff(j))
 enddo

 do i = 1, nlin_tell
	mol = get_id(imol(i))
	nai = n_air(i)
	do j = 1, nlayer
		ab = abudance(mol.id_molecule, j)
		gmix = (1.d0 - ab)*gamma_air(i) + ab*gamma_self(i)
		powTr = exp(nai * logTrefOverT(j))
		lorcoef(i,j) = gmix * pr(j) * powTr
	enddo
 enddo

 deallocate(pr, logTrefOverT)

 if(wr == 1) write(6,'(1h+,a)')' Preparation of Lorentz    coefficient by layers completed successfully                       '

end subroutine GET_LORCOEF

