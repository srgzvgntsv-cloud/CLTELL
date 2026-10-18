subroutine LOAD_OBSERVED
 use OBSERVED
 use MODELS
 implicit none

 integer :: bitpix, naxis1, nb, if_FITS
 real(8) :: dwave, w0, w, r, first_pix
 integer :: ios, i, k, n, length_Header
 byte, allocatable :: buf(:)
 real(4), allocatable :: fl(:)
 character(8) :: keyword
 character(160) :: string
 character(80)  :: string1, string2, string3, string4, string5

 if_FITS = 0
 k = index(observed_names,'.fits')
 if(k > 0) if_FITS = 1

 if(allocated(wobs)) deallocate(wobs)
 if(allocated(robs)) deallocate(robs)
 if(allocated(model)) deallocate(model)
 if(allocated(model_tmp)) deallocate(model_tmp)

 select case(if_FITS)

	case(0)
		open(10,file=trim(adjustl(synth_path))//observed_names,status='old',iostat=ios)
		if(ios /= 0) then
			write(*,'(" This file not exist: ",a)') trim(adjustl(synth_path))//observed_names
			stop
		endif
		nobs = 0
		do
			read(10,*,iostat=ios) w, r
			if(ios /= 0) exit
			if(w < wfirst) cycle
			if(w > wlast) exit
			nobs = nobs + 1
		enddo
		if(nobs <= 0) then
			close(10)
			stop 'LOAD_OBSERVED: no points in selected regions'
		endif
		allocate(wobs(nobs), robs(nobs), model(nobs), model_tmp(nobs))
		rewind(10)
		nobs = 0
		do
			read(10,*,iostat=ios) w, r
			if(ios /= 0) exit
			if(w < wfirst) cycle
			if(w > wlast) exit
			nobs = nobs + 1
			wobs(nobs) = w
			robs(nobs) = r
			model(nobs) = 0.d0
			model_tmp(nobs) = 0.d0
		enddo
		close(10)

	case(1)
		open(10,file=trim(adjustl(synth_path))//observed_names,status='old',iostat=ios,access='direct',form='binary',recl=80)
		if(ios /= 0) then
			write(*,'(" This file not exist: ",a)') trim(adjustl(synth_path))//observed_names
			stop
		endif
		string(1:80) = ''
		n = 0
		do
			n = n + 1
			read(10,rec=n) string
			keyword = string(1:8)
			k = index(string,'/')
			if(k == 0) k = 81
			select case(keyword)
				case('BITPIX')
					read(string(10:k-1),*) bitpix
				case('CRVAL1')
					read(string(10:k-1),*) w0
				case('CRPIX1')
					read(string(10:k-1),*) first_pix
				case('NAXIS1')
					read(string(10:k-1),*) naxis1
				case('CDELT1')
					read(string(10:k-1),*) dwave
				case('END')
					exit
			end select
		enddo
		length_Header = n*80
		if(mod(length_Header,2880) == 0) then
			nb = length_Header
		else
			nb = (length_Header/2880 + 1)*2880
		endif
		close(10)
		if(first_pix == 0.d0) first_pix = 1.d0
		w0 = w0 - first_pix*dwave
		allocate(buf(nb), fl(naxis1))
		open(10,file=trim(adjustl(synth_path))//observed_names,form='binary',convert='BIG_ENDIAN',status='old',iostat=ios)
		if(ios /= 0) stop 'LOAD_OBSERVED: cannot reopen FITS’'
		read(10) buf, fl
		close(10)
		deallocate(buf)
		nobs = 0
		do i = 1, naxis1
			w = w0 + dble(i)*dwave
			if(w < wfirst) cycle
			if(w > wlast) exit
			nobs = nobs + 1
		enddo
		if(nobs <= 0) then
			deallocate(fl)
			stop 'LOAD_OBSERVED: no FITS points in selected regions'
		endif
		allocate(wobs(nobs), robs(nobs), model(nobs), model_tmp(nobs))
		nobs = 0
		do i = 1, naxis1
			w = w0 + dble(i)*dwave
			if(w < wfirst) cycle
			if(w > wlast) exit
			nobs = nobs + 1
			wobs(nobs) = w
			robs(nobs) = dble(fl(i))
			model(nobs) = 0.d0
			model_tmp(nobs) = 0.d0
		enddo
		deallocate(fl)

 end select

 robs = 1.d0 - robs

 open(10,file='Vega/vega_cont.dat',status='old',iostat=ios)
 n = 0
 do
	read(10,*,iostat=ios)w,r
	if(ios /= 0) exit
	if(w < wfirst) cycle
	if(w > wlast) exit
	n = n + 1
 enddo
 allocate(continuum(n))
 rewind(10)
 n = 0
 do
	read(10,*,iostat=ios) w, r
	if(ios /= 0) exit
	if(w < wfirst) cycle
	if(w > wlast) exit
	n = n + 1
	continuum(n) = r
 enddo
 continuum = 1.-continuum
 close(10)


 call GAUS_FOR_RESOLUTION

 write(string1,'(i8)') nobs
 write(string2,'(f10.3)') wobs(1)
 write(string3,'(f10.3)') wobs(nobs)
 write(string4,'(f7.3)') minval(1.d0-robs)
 write(string5,'(f6.3)') maxval(1.d0-robs)
 write(*,'("  Spectrum is loaded. Npoint: [",a,"], Wave (A): [",a," <> ",a,"], Intensity: [",a," <> ",a,"]",/)') &
 trim(adjustl(string1)), trim(adjustl(string2)), trim(adjustl(string3)), trim(adjustl(string4)), trim(adjustl(string5))

end subroutine LOAD_OBSERVED

   

