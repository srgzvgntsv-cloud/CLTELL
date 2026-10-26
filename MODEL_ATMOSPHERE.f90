subroutine READ_MODEL_EARLTH_ATMOSPHERE
 use TELLURICA
 use CONSTANTS
 use EARTHLY_MODEL_ATM
 use OBSERVED
 use TELL_ABU_FIT_MOL
 implicit none

 integer :: i, j, nold, ios, nmol
 real(8), allocatable :: col(:)
 real(8) :: ab           

 nold = size(height)

 if(allocated(teff)) deallocate(teff)
 if(allocated(press)) deallocate(press)
 if(allocated(abudance)) deallocate(abudance)
 if(allocated(abudance_mol)) deallocate(abudance_mol)

 allocate(teff(nlayer), press(nlayer), abudance(9,nlayer), abudance_mol(9,nlayer))
 allocate(col(nlayer))

 call QINTER(height, temperature(1:nold), nold, altitude, teff, nlayer)
 call QINTER(height, pressure(1:nold), nold, altitude, press, nlayer)

 do i = 1, 9
	call QINTER(height, vmr(i,1:nold), nold, altitude, abudance(i,1:nlayer), nlayer)
 enddo

 do i = 1, nlayer
	col(i) = N0 * (press(i)/Pref) * (273.15d0/teff(i)) * astep * 1.d2
 enddo

 do j = 1, 9
	do i = 1, nlayer
		abudance_mol(j,i) = col(i) * abudance(j,i)
	enddo
 enddo

 deallocate(col)

 abu_fit = .false.
 open(55,file = trim(adjustl(synth_path))//'scale_abudance.dat', status = 'old', iostat = ios)
 if(ios /= 0) then
	abu_fit = .true.
 else
	do i = 1, 9
		abu_fit = .true.
		read(55,*,iostat = ios) nmol, ab
		if(ios /= 0) exit
		if(ab < 0.d0) exit
		scale_mol(i) = ab 
		abu_fit = .false.
	enddo
	if(abu_fit == .true.) then
		scale_mol = 1.d0
		write(*,*)'Missing lines in scale_abudance.dat or negative scale_mol values'
	endif
	close(55)
	call APPLY_MOL_SCALE_TO_ABUDANCE	
 endif

! write(6,"(1h+' Prepared all parameters from the atmosphere model ')")

end subroutine READ_MODEL_EARLTH_ATMOSPHERE


subroutine GET_LAYERS
 use TELLURICA
 implicit none

 integer :: k
 real(8) :: z0, ztop
 character(80)  :: string1, string2, string3

 ztop = 120000.d0
 z0 = height0
 if(z0 < 0.d0) z0 = 0.d0
 if(z0 > ztop) z0 = ztop

 if(allocated(altitude)) deallocate(altitude)

 nlayer = int((ztop - z0)/astep) + 1
 if(nlayer < 1) nlayer = 1

 allocate(altitude(nlayer))

 do k = 1, nlayer
	altitude(k) = z0 + dble(k-1)*astep
 enddo
 altitude(nlayer) = ztop

 write(string1,'(i5)') nlayer
 write(string2,'(f6.2)') altitude(1)/1.d3
 write(string3,'(f6.2)') altitude(nlayer)/1.d3
 write(*,'("  Atmosphere is divided into [",a,"] layers from [",a,"] km to [",a,"] km",/)') &
 trim(adjustl(string1)), trim(adjustl(string2)), trim(adjustl(string3))

end subroutine GET_LAYERS