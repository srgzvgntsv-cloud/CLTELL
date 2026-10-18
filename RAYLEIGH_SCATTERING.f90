subroutine READ_RELATIVE_HUMIDITY
 use TELLURICA
 use CONSTANTS
 use OBSERVED, only : file_relative_humidity
 implicit none
 integer :: n, ios, i
 real(8) :: h, r
 real(8) :: Es_hPa, Es_Pa, Pw_Pa, RHpct, RHfrac
 real(8) :: dz_cm, kB, n_m3, n_cm3, ln10, expo
 real(8), allocatable :: hr(:), rh(:)

 ln10 = 2.302585092994046d0
 dz_cm = astep * 1.d2
 kB = (rgas/1000.d0) / Na
 if(allocated(relative_humidity)) deallocate(relative_humidity)
 allocate(relative_humidity(nlayer))
 relative_humidity = 0.d0
 open(9, file=file_relative_humidity, status='old', iostat=ios)
 if(ios == 0) then
	n = 0
	do
		read(9,*,iostat=ios) h, r
		if(ios > 0) cycle
		if(ios == -1) exit
		n = n + 1
	enddo
	if(n >= 2) then
		rewind(9)
		allocate(hr(n), rh(n))
		n = 0
		ios = 0
		do
			read(9,*,iostat=ios) h, r
			if(ios > 0) cycle
			if(ios == -1) exit
			n = n + 1
			hr(n) = h
			rh(n) = r
		enddo
		rh(1) = rh(2)
		call LINTER(hr, rh, n, altitude, relative_humidity, nlayer)
		do i = 1, nlayer
			RHpct = relative_humidity(i)
			if(RHpct < 0.d0) RHpct = 0.d0
			if(RHpct > 100.d0) RHpct = 100.d0
			RHfrac = RHpct / 100.d0
			expo = (7.5d0*teff(i) - 2048.625d0) / (teff(i) - 35.85d0)
			Es_hPa = 6.1078d0 * exp(ln10*expo)
			Es_Pa  = Es_hPa * 100.d0
			Pw_Pa = RHfrac * Es_Pa
			if(teff(i) > 0.d0) then
				n_m3  = Pw_Pa / (kB * teff(i))
				n_cm3 = n_m3 / 1.d6
				abudance_mol(1,i) = n_cm3 * dz_cm
			endif
		enddo
		deallocate(hr, rh)
	endif
	close(9)
 else
	close(9)
 endif
 do i = 1, nlayer
	expo = (7.5d0*teff(i) - 2048.625d0) / (teff(i) - 35.85d0)
	Es_hPa = 6.1078d0 * exp(ln10*expo)
	Es_Pa = Es_hPa * 100.d0
	if(Es_Pa > 0.d0 .and. teff(i) > 0.d0 .and. dz_cm > 0.d0) then
		n_cm3 = abudance_mol(1,i) / dz_cm
		n_m3  = n_cm3 * 1.d6
		Pw_Pa = n_m3 * kB * teff(i)
		relative_humidity(i) = 100.d0 * Pw_Pa / Es_Pa
	else
		relative_humidity(i) = 0.d0
	endif
	if(relative_humidity(i) < 0.d0) relative_humidity(i) = 0.d0
	if(relative_humidity(i) > 100.d0) relative_humidity(i) = 100.d0
 enddo

end subroutine READ_RELATIVE_HUMIDITY



 subroutine GET_RAYSCAT(wr)
 use TELLURICA
 use OBSERVED
 implicit none

 integer :: i, j, wr
 real(8), allocatable :: abud_cm(:), invwl4(:)
 real(8) :: factor, wl4

 if(allocated(rayscat)) deallocate(rayscat)
 if(allocated(rayscat_cont_tau)) deallocate(rayscat_cont_tau)
 allocate(rayscat(nobs,nlayer),rayscat_cont_tau(nobs))

 allocate(abud_cm(nlayer))
 do i = 1, nlayer
	abud_cm(i) = sum(abudance_mol(1:9,i))/1.d15
	if(abud_cm(i) <= 0.d0) abud_cm(i) = 1.d-300
 enddo
 allocate(invwl4(nobs))
 do j = 1, nobs
	invwl4(j) = 1.d16 / (wobs(j)*wobs(j)*wobs(j)*wobs(j))
 enddo
 do i = 1, nlayer
	do j = 1, nobs
		factor = (82.68340447d0 * (refind(j,i)*refind(j,i) - 1.d0)**2) / (abud_cm(i)) * invwl4(j)
		rayscat(j,i) = 1.061d0 * factor * 1.d8
	enddo
 enddo
 rayscat_cont_tau = sum(rayscat, dim=2)
 deallocate(invwl4, abud_cm)

 if(wr == 1) write(6,"(1h+' Preparation of Rayleigh scattering coefficient by wavelengths has been successfully completed')")

 end subroutine GET_RAYSCAT



 subroutine REFRACTIVE_INDEX
! The equation for the index of refraction n is taken from: Edlen,K.(1966) The Refractive Index of Air, Metrologia 2:12
 use TELLURICA
 use OBSERVED
 implicit none
 real(8) :: Pw_hPa, Es_hPa, RHfrac, ln10, expo
 real(8), allocatable :: n2(:), n4(:), n1arr(:), n3arr(:)
 real(8) :: wave_number
 integer :: i, j

 ln10 = 2.302585092994046d0

 if(allocated(n2)) deallocate(n2)
 if(allocated(n4)) deallocate(n4)
 allocate(n2(nlayer), n4(nlayer))

 if(allocated(refind)) deallocate(refind)
 allocate(refind(nobs,nlayer))

 allocate(n1arr(nobs), n3arr(nobs))

 do i = 1, nlayer
	expo = (7.5d0*teff(i)-2048.625d0)/(teff(i)-35.85d0)
	Es_hPa = 6.1078d0 * exp(ln10*expo)
	RHfrac = relative_humidity(i)/100.d0
	if(RHfrac < 0.d0) RHfrac = 0.d0
	if(RHfrac > 1.d0) RHfrac = 1.d0
	Pw_hPa = RHfrac * Es_hPa
	n2(i) = ((press(i)/100.d0 - Pw_hPa)/1013.25d0) * 296.15d0/teff(i)
	n4(i) = Pw_hPa/1013.25d0
 enddo
 do j = 1, nobs
	wave_number = 1.d8/wobs(j)
	n1arr(j) = 83.43d0 + 185.08d0/(1.d0-(wave_number/1.14d5)**2) + 4.11d0/(1.d0-(wave_number/6.24d4)**2)
	n3arr(j) = 4349.d0 - (wave_number/1.7d4)**2
 enddo
 do i = 1, nlayer
	do j = 1, nobs
		refind(j,i) = (n1arr(j)*n2(i) + n3arr(j)*n4(i))/1.d6 + 1.d0
	enddo
 enddo
 deallocate(n1arr, n3arr, n2, n4)
 if(allocated(relative_humidity)) deallocate(relative_humidity)
! write(6,"(1h+' Preparation of Rayleigh scattering coefficient by layers has been successfully completed')")
 end subroutine REFRACTIVE_INDEX