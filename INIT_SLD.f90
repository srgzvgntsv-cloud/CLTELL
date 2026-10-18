subroutine INIT_SLD
 use OBSERVED
 use SLD
 use LINES

 implicit none

 integer :: i, j, k, m, isld
 integer :: limit_index, intensity_group
 integer :: nsld_tell, star_sld
 integer :: number_of_lines

 real(8) :: gaussx(3)

 character(4), allocatable :: species(:)
 character(4) :: cha

 integer, allocatable :: profile_group(:)

 allocate(species(nlines))
 allocate(profile_group(nlines))

 if(allocated(line_sld)) deallocate(line_sld)
 allocate(line_sld(nlines))

 nsld = 0
 nsld_tell = 0
 star_sld = 0

 line_sld = 0
 species = '    '
 profile_group = 0

 do i = 1, nlines

	if(line_family(i) /= 1) cycle

	intensity_group = 1

	do m = 1, n_sld_group_species
		if(line_species(i) /= sld_group_species(m)) cycle
		intensity_group = sld_group_count(m)+1
		do j = 1, sld_group_count(m)
			limit_index = sld_group_first(m)+j-1
			if(line_intensity(i) <= sld_group_limits(limit_index)) then
				intensity_group = j
				exit
			endif
		enddo
		exit
	enddo

	k = 0

	do j = 1, nsld_tell
		if(line_species(i) == species(j) .and. intensity_group == profile_group(j)) then
			k = j
			exit
		endif
	enddo

	if(k == 0) then
		nsld = nsld+1
		nsld_tell = nsld_tell+1
		species(nsld_tell) = line_species(i)
		profile_group(nsld_tell) = intensity_group
		k = nsld_tell
	endif

	line_sld(i) = k
 enddo

 if(allocated(sld_group_species)) deallocate(sld_group_species)
 if(allocated(sld_group_first)) deallocate(sld_group_first)
 if(allocated(sld_group_count)) deallocate(sld_group_count)
 if(allocated(sld_group_limits)) deallocate(sld_group_limits)

 n_sld_group_species = 0
 n_sld_group_limits = 0

 do i = 1, nlines
	if(line_family(i) == 2) then
		nsld = nsld+1
		star_sld = nsld
		exit
	endif
 enddo

 if(star_sld > 0) then
	do i = 1, nlines
		if(line_family(i) == 2) line_sld(i) = star_sld
	enddo
 endif

 write(*,'("  Reading data from HITRAN file was successful, extracted lines:  [",i5,"]")') nlin_tell

 do i = 1, nsld_tell
	number_of_lines = 0; cha = '    '; isld = 0
	do j = 1, nlines
		if(line_family(j) == 1 .and. line_sld(j) == i) then
			number_of_lines = number_of_lines+1
			cha = line_species(j)
			isld = line_sld(j)
		endif
	enddo
	if(i == 1) then
		write(*,'(55x,"Of these:  [",i5,"] - ",a4,", #SLD:",i3)') number_of_lines, cha, isld
	else
		write(*,'(66x,"[",i5,"] - ",a4,", #SLD:",i3)') number_of_lines, cha, isld
	endif
 enddo
 
 do i = 1, nlines
	if(line_family(i) == 2) then
		isld = line_sld(i)
		exit
	endif
 enddo
 write(*,'("  Extracted from the star lines database:                         [",i5,"] -       #SLD:",i3)') nlin_star, isld
 write(*,*)

 deallocate(species)
 deallocate(profile_group)

 if(nsld <= 0) then
	write(*,'("  Error: no SLD profiles were created")')
	stop
 endif

 if(allocated(SLD_velocity)) deallocate(SLD_velocity)
 if(allocated(SLD_intensity)) deallocate(SLD_intensity)

 allocate(SLD_velocity(nsld,SLD_npoint))
 allocate(SLD_intensity(nsld,SLD_npoint))

 SLD_velocity = 0.d0
 SLD_intensity = 0.d0

 gaussx(1) = 5.d-1
 gaussx(2) = 0.d0
 gaussx(3) = SLD_npoint/32.d0

 do i = 1, nsld_tell
	SLD_velocity(i,SLD_npoint/2+1) = 0.d0
	do j = 1, SLD_npoint/2
		SLD_velocity(i,SLD_npoint/2+1-j) = -dble(j)*step_tell
		SLD_velocity(i,SLD_npoint/2+1+j) =  dble(j)*step_tell
	enddo
	call FUNCTION_GAUSS(SLD_npoint,3,gaussx,SLD_velocity(i,1:SLD_npoint),SLD_intensity(i,1:SLD_npoint))
 enddo

 if(star_sld > 0) then
	SLD_velocity(star_sld,SLD_npoint/2+1) = 0.d0
	do i = 1, SLD_npoint/2
		SLD_velocity(star_sld,SLD_npoint/2+1-i) = -dble(i)*step_star
		SLD_velocity(star_sld,SLD_npoint/2+1+i) =  dble(i)*step_star
	enddo
	call INIT_SLD_STAR_ROT(star_sld,vsini_user,0.6d0)
 endif

 call SETUP_LINE_CORRECTION

end subroutine INIT_SLD


subroutine INIT_SLD_STAR_ROT(isld, vsini, eps)
 use SLD
 implicit none
 integer, intent(in) :: isld
 real(8), intent(in) :: vsini
 real(8), intent(in), optional :: eps

 integer :: j, mid
 real(8) :: e, v, x, y, dv, pi, norm

 if(isld < 1 .or. isld > nsld) stop 'INIT_SLD_STAR_ROT: isld out of range'

 pi = 3.1415926535897932384626433832795d0
 e = 0.6d0
 if(present(eps)) e = eps
 if(e < 0.d0) e = 0.d0
 if(e > 1.d0) e = 1.d0

 if(allocated(SLD_intensity) .eqv. .false.) stop 'INIT_SLD_STAR_ROT: SLD_intensity not allocated'
 if(allocated(SLD_velocity) .eqv. .false.) stop 'INIT_SLD_STAR_ROT: SLD_velocity not allocated'
 if(SLD_npoint < 3) stop 'INIT_SLD_STAR_ROT: SLD_npoint too small'

 dv = abs(SLD_velocity(isld,2) - SLD_velocity(isld,1))
 if(dv <= 0.d0) stop 'INIT_SLD_STAR_ROT: invalid SLD_velocity step'

 SLD_intensity(isld,1:SLD_npoint) = 0.d0

 if(vsini <= dv) then
	mid = SLD_npoint/2 + 1
	SLD_intensity(isld,mid) = 1.d0/dv
	return
 endif

 do j = 1, SLD_npoint
	v = SLD_velocity(isld,j)
	x = v / vsini
	if(abs(x) < 1.d0) then
		y = 1.d0 - x*x
		SLD_intensity(isld,j) = ( 2.d0*(1.d0-e)*sqrt(y) + 0.5d0*pi*e*y ) / ( pi*vsini*(1.d0 - e/3.d0) )
	else
		SLD_intensity(isld,j) = 0.d0
	endif
 enddo

 norm = sum(SLD_intensity(isld,1:SLD_npoint)) * dv
 if(norm > 0.d0) then
	SLD_intensity(isld,1:SLD_npoint) = SLD_intensity(isld,1:SLD_npoint) / norm
 else
	mid = SLD_npoint/2 + 1
	SLD_intensity(isld,mid) = 1.d0/dv
 endif

end subroutine INIT_SLD_STAR_ROT