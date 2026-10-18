subroutine Hprofile(wavelength_c,nHyd,i1,i2) ! Calculate profile of given Hydrogen line
 use synth_data
 use HlineTab
 use radiation_field
 use atmosphere_model
 implicit real(8) (a-h,o-z)
 character(4) id_l 
 character :: id_ser(4) = (/'L','H','P','B'/)
 real(8), parameter :: c = 2.997925d5
 real(8) :: dwH = 300.d0 ! maximum distyance from line center
 integer(4), parameter :: nHydr = 60001 ! maximum number of points in Hydrogen line profile buffers 
 real(8) wbuf(nHydr),profbuf(nHydr) ! Place for save calculated profile
 real(8), allocatable :: wprofile(:),profile(:)

 nHyd=0
 icentre=nHydr/2+1
 wave=wavelength_c-dwH-wavelength_step
 do i=1,nHydr
   wave=wave+wavelength_step  
   call update_radiation_frequency(wave)
   call get_continuum_opacity(wave)
   call  get_continuum_flux(wave)         ! Interpolate early calculated continuum flux and intensities to the given wavelength
   Hnu_c=Hnu(1)
   if(wave > 3600.and.wave < 4200)then
      call kapp
      Aline=0.d0
      call Josh
      Hnu_c=Hnu(1)
   endif
   call Hline_profile(wave,i1,i2)
   Aline=AHline
   call Josh

   residual=Hnu(1)/Hnu_c
   if(i == icentre)residual_c=residual
   if(nHyd == 0)wleft=wave
   if(wleft > wavelength_last.or.(2.d0*wavelength_c-wleft) < wavelength_first)exit ! Eliminate lines outside wavelengths region
   nHyd=nHyd+1
   wbuf(nHyd)=wave
   profbuf(nHyd)=residual
 enddo

 if(nHyd == 0)return
 wright=wave

 k = 0
 do i = 1,nHyd
	if(wbuf(i) >= wavelength_first .and. wbuf(i) <= wavelength_last) k = k+1
 enddo
 allocate (wprofile(k),profile(k))
 k = 0
 do i = 1,nHyd
	if(wbuf(i) >= wavelength_first .and. wbuf(i) <= wavelength_last) then
		k = k+1
		wprofile(k) = wbuf(i)
		profile(k) = profbuf(i)
    endif
 enddo
 nHyd = k

 write(id_l,'(a4)')'H  1'
 write(132)wavelength_c,id_l,residual_c,nHyd,wprofile(1),wprofile(nHyd)
 write(132)wprofile,profile
 deallocate (wprofile,profile)   

end


subroutine Hline_profile(wave,nser,line) ! Calculate individual Hydrogen line profile
 use radiation_field
 use atmosphere_model
 use temporary_arrays
 use HlineTab
 implicit real(8) (a-h,o-z)
 real(8) Prof,nmax
 real(8), allocatable, save :: wwij(:,:)
 integer :: itemp1 = 0

  if (itemp1 == 0 ) then ! Independent upon wavelengths
    itemp1 = 1
    do j = 1, nrhox
      Mlast(j) = 1100.d0 / xne(j)**(2.d0/15.d0) ! last visible line
    enddo
    do n = 1, nHydSeries
      do j = 1, nrhox
        ep      = 13.595d0 - 13.595d0 / (n*n)
        g       = 2.d0 * (n*n)
        xHh(n,j) = exp( - ep / Tkev(j) ) * g * xDens(j,1,1) / rho(j)
      enddo
    enddo

  nlevels=nHydSeries+nHydLines
  if(.not.allocated(wwij))allocate (wwij(nlevels,Nrhox))

  do i = 1, nlevels
   do j = 1, Nrhox
    wwij(i,j) = occupationProbability(1d0,0d0,dble(i),0d0,t(j),xNe(j),xDens2(J,1,1)*0.d0)
   end do
  end do 
 endif

   Ahline = 0.d0
 
  do j = 1, nrhox
    wi = wwij(nser,j)
       nup = nser + line
       wj  =  wwij(nup,j)
        del = abs( wave - wAir(Line,nser) )
        vnm = 2.99792458d18 / wAir(Line,nser)
        call VCS(nser,del,Prof,Line,j)    ! Vidal, Cuper, Smith
        if ( nser == 2 .and. del > 0.2d0 .and. Line <= 4 ) then
         reson   = ( wAir(Line,nser)**2 ) * 1.71d-27 * xDens(j,1,1) * 2.d0
         a       = radamp(Line) + reson
         pr2     = a / 3.1416d0 / ( (del**2) + a*a )
         Prof = ( Prof + pr2 ) * fHyd(Line,nser)
        else
         Prof = Prof * fHyd(Line,nser)
        endif
       stiml     = 1.d0 - exp( ( - 4.7993d-11 * vnm ) / T(j) )
       Prof   = Prof * stiml * xHh(nser,j) / (vnm**2) * 7.957d16
       Ahline(j) =Prof*wj
       Sline(j) = Bnu(j)
  enddo

end


subroutine line_profile(ww1,ww2,wc,id_l,rc,Stark_switch) ! Calculate non Hydrogen line profile from ww1 to ww2
 use synth_data
 use radiation_field
 use atmosphere_model
 use lines_data
 use Hetab
 implicit real(8) (a-h,o-z)
 real(8) kappa
 character(4) id_l
 real(8), allocatable :: wprofile(:),profile(:)
 real(8) wbuf(30000),pbuf(30000)
 real(8) :: c=2.997925d5
 integer(2) :: Stark_switch

 if(ww2 < wavelength_first .or. ww1 > wavelength_last)return ! skip lines outside wavelengths region

 rlimit=0.9998d0 ! profile intensity limit
 m=0
 w=wc+wavelength_step
 do  ! Calculate left part of line profile
   w=w-wavelength_step
   if(w <= wavelength_first) exit
   do j=1,Nrhox ! Calculate line absorbtion vs. depth
     x=w-wc            ! Distance from line centre
     if(Stark_switch == 0)then
       xDop=abs(x)/10.d0*Dopdump(j)           ! used in nm
	   kappa=al_int(j)*Voigt(xDop,Adump(j))
     else
	   call He_prof(Stark_switch,j,x,kappa)
     endif
	 Aline(j)=kappa                              
   enddo
   call update_radiation_frequency(w)
   call get_continuum_opacity(w)
   call  get_continuum_flux(w)         ! Interpolate early calculated continuum flux and intensities to the given wavelength
   Hnu_c=Hnu(1)
   call Josh
   r=Hnu(1)/Hnu_c
   if(r >= rlimit)exit
   m=m+1
   pbuf(m)=r
   wbuf(m)=w ! velocities as abscices
 enddo
 ww1=w+wavelength_step; m1=m

 w=wc
 do  ! Calculate right part of line profile
   w=w+wavelength_step
   if(w >= wavelength_last) exit
   do j=1,Nrhox ! Calculate line absorbtion vs. depth
     x=w-wc            ! Distance from line centre
     if(Stark_switch == 0)then
       xDop=abs(x)/10.d0*Dopdump(j)           ! used in nm
	   kappa=al_int(j)*Voigt(xDop,Adump(j))
     else
	   call He_prof(Stark_switch,j,x,kappa)
     endif
	 Aline(j)=kappa                              
   enddo
   call update_radiation_frequency(w)
   call get_continuum_opacity(w)
   call  get_continuum_flux(w)         ! Interpolate early calculated continuum flux and intensities to the given wavelength
   Hnu_c=Hnu(1)
   call Josh
   r=Hnu(1)/Hnu_c
   if(r >= rlimit)exit
   m=m+1
   pbuf(m)=r
   wbuf(m)=w ! velocities as abscices
 enddo
 ww2 = w-wavelength_step

 if(m <= 3)return ! eliminate too slow lines

 allocate (wprofile(m),profile(m))
 do i=1,m1
   wprofile(i)=wbuf(m1-i+1)
   profile(i)=pbuf(m1-i+1)
 enddo
 wprofile(m1+1:m)=wbuf(m1+1:m); profile(m1+1:m)=pbuf(m1+1:m)
 nlines_db=nlines_db+1
 write(131)wc,id_l,rc,m,ww1,ww2
 write(131)wprofile,profile

 deallocate (wprofile,profile)

end

subroutine update_radiation_frequency(wave)
 use radiation_field
 use atmosphere_model
 implicit none

 real(8), intent(in) :: wave
 integer :: k

 freq = 2.997925d18 / wave
 freqlg = log(freq)

 do k = 1, Nrhox
   EhvkT(k) = exp(-freq * hkT(k))
   stim(k) = 1.d0 - EhvkT(k)
   Bnu(k) = 1.47439d-47 * freq**3 * EhvkT(k) / stim(k)
 enddo
end subroutine update_radiation_frequency


subroutine form_database
 use synth_data
 use HlineTab
 implicit none

 character(4) :: id_l
 integer(4) :: i, j
 integer(4) :: nHydr, n_Hydrogen
 integer(4) :: buffer_size
 integer(4) :: number_to_copy
 real(8) :: wave_c, rc, w1, w2
 real(8), allocatable :: wprofile(:), profile(:)
 real(8), parameter :: dwH = 600.d0

 n_Hydrogen = 0

 do i = 1, nHydSeries
	do j = 1, nHydLines
		if(wAir(j,i) + dwH < wavelength_first) cycle
		if(wAir(j,i) - dwH > wavelength_last) cycle
		call Hprofile(wAir(j,i), nHydr, i, j)
		if(nHydr <= 0) cycle
		n_Hydrogen = n_Hydrogen + 1
	enddo
 enddo

 rewind(131)
 rewind(132)

 write(134) nlines_db, n_Hydrogen

 buffer_size = 0

 do i = 1, nlines_db
	read(131) wave_c, id_l, rc, number_to_copy, w1, w2
	if(number_to_copy <= 0) then
		write(*,*) 'form_database: invalid profile size: ', number_to_copy
		stop
	endif
	if(number_to_copy > buffer_size) then
		if(allocated(wprofile)) deallocate(wprofile)
		if(allocated(profile)) deallocate(profile)
		buffer_size = number_to_copy
		allocate(wprofile(buffer_size))
		allocate(profile(buffer_size))
	endif

	read(131) wprofile(1:number_to_copy), profile(1:number_to_copy)

	write(134) wave_c, id_l, rc, number_to_copy, w1, w2

	write(133) wprofile(1:number_to_copy), profile(1:number_to_copy)
 enddo

 do i = 1, n_Hydrogen
	read(132) wave_c, id_l, rc, number_to_copy, w1, w2

	if(number_to_copy <= 0) then
		write(*,*) 'form_database: invalid hydrogen profile size: ', number_to_copy
		stop
	endif

	if(number_to_copy > buffer_size) then
		if(allocated(wprofile)) deallocate(wprofile)
		if(allocated(profile)) deallocate(profile)
		buffer_size = number_to_copy
		allocate(wprofile(buffer_size))
		allocate(profile(buffer_size))
	endif

	read(132) wprofile(1:number_to_copy), profile(1:number_to_copy)

	write(134) wave_c, id_l, rc, number_to_copy, w1, w2

	write(133) wprofile(1:number_to_copy), profile(1:number_to_copy)
 enddo

 if(allocated(wprofile)) deallocate(wprofile)
 if(allocated(profile)) deallocate(profile)

 close(131, status = 'delete')
 close(132, status = 'delete')

 rewind(133)
 rewind(134)

 call get_alphas

end subroutine form_database


   
subroutine get_alphas
 use synth_data
 implicit none

 real(8), allocatable :: vprof(:), rprof(:)
 real(8), allocatable :: alphas(:), wprof(:)
 real(8), allocatable :: alphas_sp(:)
 real(8), allocatable :: wlines(:), simple_sum(:)
 real(8), allocatable :: cRs(:), w1s(:), w2s(:)

 integer(4), allocatable :: npointss(:), nbl(:)
 character(4), allocatable :: temp_idents(:)

 character(8), parameter :: database_magic = 'SPROFDB2'
 integer(4), parameter :: database_version = 2

 integer(4) :: n_lines, n_hyd, nlin
 integer(4) :: max_points
 integer(4) :: np
 integer(4) :: i, j
 integer(4) :: ilo, ihi
 integer(4) :: id, map1
 integer(4) :: ios

 logical :: file_exists

 real(8) :: profint, vp
 real(8), parameter :: vc = 2.997925d5
 real(8), parameter :: min_profile_sum = 1.d-30

 synth_R = 1.d0 - synth_R

 read(134, iostat = ios) n_lines, n_hyd

 if(ios /= 0) then
	write(*,*) 'get_alphas: cannot read database header'
	stop
 endif

 nlin = n_lines + n_hyd

 if(n_lines < 0 .or. n_hyd < 0 .or. nlin < 0) then
	write(*,*) 'get_alphas: invalid number of lines'
	stop
 endif

 allocate(cRs(nlin))
 allocate(w1s(nlin))
 allocate(w2s(nlin))
 allocate(npointss(nlin))
 allocate(temp_idents(nlin))
 allocate(wlines(nlin))

 do j = 1, nlin
	read(134, iostat = ios) wlines(j), temp_idents(j), cRs(j), npointss(j), w1s(j), w2s(j)

	if(ios /= 0) then
		write(*,*) 'get_alphas: cannot read line information: ', j
		stop
	endif

	if(npointss(j) <= 0) then
		write(*,*) 'get_alphas: invalid profile size: ', j, npointss(j)
		stop
	endif
 enddo

 close(134, status = 'delete')

 if(nlin > 0) then
	max_points = maxval(npointss)
 else
	max_points = 1
 endif

 allocate(wprof(max_points))
 allocate(rprof(max_points))
 allocate(vprof(max_points))
 allocate(alphas(max_points))

 allocate(simple_sum(number_of_wavelengths))
 allocate(alphas_sp(number_of_wavelengths))
 allocate(nbl(number_of_wavelengths))

 alphas_sp = 1.d0
 simple_sum = 0.d0
 nbl = 0

 do j = 1, nlin
	np = npointss(j)

	read(133, iostat = ios) wprof(1:np), rprof(1:np)

	if(ios /= 0) then
		write(*,*) 'get_alphas: cannot read temporary profile: ', j
		stop
	endif

	rprof(1:np) = 1.d0 - rprof(1:np)

	call find_wavelength_range(synth_wl, number_of_wavelengths, w1s(j), w2s(j), ilo, ihi)

	if(ilo <= ihi) then
		nbl(ilo:ihi) = nbl(ilo:ihi) + 1

		do i = ilo, ihi
			id = map1(wprof, rprof, np, synth_wl(i), profint, 1)

			if(profint > 0.d0) then
				simple_sum(i) = simple_sum(i) + profint
			endif
		enddo
	endif
 enddo

 rewind(133)

 do i = 1, number_of_wavelengths
	if(nbl(i) > 1 .and. simple_sum(i) > min_profile_sum) then
		alphas_sp(i) = synth_R(i) / simple_sum(i)
	endif
 enddo

 deallocate(nbl)
 deallocate(simple_sum)

 inquire(file = trim(adjustl(proffilename_inf)), exist = file_exists)

 if(file_exists) then
	open(646, file = trim(adjustl(proffilename_inf)), form = 'unformatted', status = 'old', iostat = ios)

	if(ios == 0) close(646, status = 'delete')
 endif

 inquire(file = trim(adjustl(proffilename_bin)), exist = file_exists)

 if(file_exists) then
	open(656, file = trim(adjustl(proffilename_bin)), form = 'unformatted', status = 'old', iostat = ios)

	if(ios == 0) close(656, status = 'delete')
 endif

 open(646, file = trim(adjustl(proffilename_inf)), form = 'unformatted', access = 'sequential', status = 'new', action = 'write', iostat = ios)

 if(ios /= 0) then
	write(*,*) 'get_alphas: cannot create information file:'
	write(*,*) trim(adjustl(proffilename_inf))
	stop
 endif

 open(656, file = trim(adjustl(proffilename_bin)), form = 'unformatted', access = 'sequential', status = 'new', action = 'write', iostat = ios)

 if(ios /= 0) then
	write(*,*) 'get_alphas: cannot create profile file:'
	write(*,*) trim(adjustl(proffilename_bin))
	close(646, status = 'delete')
	stop
 endif

 write(646) database_magic, database_version, n_lines, n_hyd

 do i = 1, nlin
	write(646) wlines(i), temp_idents(i), cRs(i), npointss(i), w1s(i), w2s(i)
 enddo

 do i = 1, nlin
	np = npointss(i)

	read(133, iostat = ios) wprof(1:np), rprof(1:np)

	if(ios /= 0) then
		write(*,*) 'get_alphas: cannot reread profile: ', i
		stop
	endif

	vprof(1:np) = ((wprof(1:np) - wlines(i)) / wlines(i)) * vc

	id = map1(synth_wl, alphas_sp, number_of_wavelengths, wprof, alphas, np)

	rprof(1:np) = 1.d0 - (1.d0 - rprof(1:np)) * alphas(1:np)

	write(656, iostat = ios) vprof(1:np), rprof(1:np)

	if(ios /= 0) then
		write(*,*) 'get_alphas: cannot write profile: ', i
		stop
	endif
 enddo

 close(133, status = 'delete')
 close(646)
 close(656)

 if(itest == 1) then
	allocate(simple_sum(number_of_wavelengths))
	simple_sum = 0.d0

	open(656, file = trim(adjustl(proffilename_bin)), form = 'unformatted', access = 'sequential', status = 'old', action = 'read', iostat = ios)

	if(ios /= 0) then
		write(*,*) 'get_alphas: cannot open profile file for test'
		stop
	endif

	do j = 1, nlin
		np = npointss(j)

		read(656, iostat = ios) vprof(1:np), rprof(1:np)

		if(ios /= 0) then
			write(*,*) 'get_alphas: error reading profile: ', j
			stop
		endif

		rprof(1:np) = 1.d0 - rprof(1:np)

		where(rprof(1:np) < 0.d0)
			rprof(1:np) = 0.d0
		endwhere

		call find_wavelength_range(synth_wl, number_of_wavelengths, w1s(j), w2s(j), ilo, ihi)

		if(ilo <= ihi) then
			do i = ilo, ihi
				vp = ((synth_wl(i) - wlines(j)) / wlines(j)) * vc

				id = map1(vprof, rprof, np, vp, profint, 1)

				if(profint > 0.d0) then
					simple_sum(i) = simple_sum(i) + profint
				endif
			enddo
		endif
	enddo

	close(656)

	inquire(file = 'testtest1.dat', exist = file_exists)

	if(file_exists) then
		open(665, file = 'testtest1.dat', &
		     status = 'old', iostat = ios)

		if(ios == 0) close(665, status = 'delete')
	endif

	open(665, file = 'testtest1.dat', &
	     status = 'new', action = 'write', iostat = ios)

	if(ios /= 0) then
		write(*,*) 'get_alphas: cannot create testtest1.dat'
		stop
	endif

	do i = 1, number_of_wavelengths
		write(665,*) synth_wl(i), 1.d0 - simple_sum(i)
	enddo

	close(665)
	deallocate(simple_sum)
 endif

 deallocate(wprof)
 deallocate(rprof)
 deallocate(vprof)
 deallocate(alphas)
 deallocate(alphas_sp)

 deallocate(wlines)
 deallocate(temp_idents)
 deallocate(cRs)
 deallocate(npointss)
 deallocate(w1s)
 deallocate(w2s)

contains

 subroutine find_wavelength_range(x, nx, xmin, xmax, first, last)
  implicit none

  integer(4), intent(in) :: nx
  real(8), intent(in) :: x(nx)
  real(8), intent(in) :: xmin, xmax

  integer(4), intent(out) :: first, last
  integer(4) :: left, right, middle

  first = 1
  last = 0

  if(nx <= 0) return
  if(xmax < x(1)) return
  if(xmin > x(nx)) return

  left = 1
  right = nx

  do while(left < right)
	middle = left + (right - left) / 2

	if(x(middle) < xmin) then
		left = middle + 1
	else
		right = middle
	endif
  enddo

  first = left

  left = first
  right = nx

  do while(left < right)
	middle = left + (right - left + 1) / 2

	if(x(middle) > xmax) then
		right = middle - 1
	else
		left = middle
	endif
  enddo

  last = left

  if(first > last) then
	first = 1
	last = 0
  endif

 end subroutine find_wavelength_range

end subroutine get_alphas
