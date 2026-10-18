subroutine READ_SLD_GROUP_LINE

 use SLD; use OBSERVED

 implicit none
 character(80)        :: line, limits_string, token
 character(4)         :: species
 real(8), allocatable :: temporary_limits(:)
 real(8), allocatable :: parsed_limits(:)
 logical              :: logic
 integer              :: i, j, ios, n
 character(4)         :: molsp(9)
 data (molsp(n),n=1,9) / 'H2O ', 'CO2 ', 'O3  ', 'N2O ', 'CO  ', 'CH4 ', 'O2  ', 'OH  ', 'N2  ' /
 integer              :: q1, q2, q2_relative
 integer              :: colon_position, comment_position, comma_position
 integer              :: number_of_limits, old_number_limits, species_length

 if(allocated(sld_group_species)) deallocate(sld_group_species)
 if(allocated(sld_group_first)) deallocate(sld_group_first)
 if(allocated(sld_group_count)) deallocate(sld_group_count)
 if(allocated(sld_group_limits)) deallocate(sld_group_limits)

 n_sld_group_species = 0
 n_sld_group_limits = 0

 if(ngroup <= 0) return

 allocate(sld_group_species(ngroup), sld_group_first(ngroup), sld_group_count(ngroup))
 sld_group_species = '    '; sld_group_first = 0; sld_group_count = 0

 do i = 1, ngroup

	line = config_line(i)
	comment_position = index(line,'!')
	if(comment_position == 1) cycle
	if(comment_position > 1) then
		line = line(1:comment_position-1)
	endif

	line = adjustl(line)
	if(len_trim(line) == 0) cycle

	q1 = index(line,"'")

	if(q1 == 0) then
		write(*,'("  Invalid SLD grouping line: ",a)') trim(line)
		write(*,'("  Expected format: ''H2O '': 0.5, 1.0")')
		stop
	endif

	q2_relative = index(line(q1+1:),"'")

	if(q2_relative == 0) then
		write(*,'("  Missing closing quote in line: ",a)') trim(line)
		stop
	endif

	q2 = q1+q2_relative

	colon_position = index(line(q2+1:),':')

	if(colon_position == 0) then
		write(*,'("  Missing colon in line: ",a)') trim(line)
		stop
	endif

	colon_position = q2+colon_position
	species_length = q2-q1-1

	if(species_length /= 4) then
		write(*,'("  Invalid molecular name in line: ",a)') trim(line)
		write(*,'("  Molecular name must contain exactly 4 characters")')
		write(*,'("  Examples: ''H2O '', ''O2  '', ''CO2 ''")')
		stop
	endif

	species = line(q1+1:q2-1)

	logic = .false.

	do j = 1, 9
		if(species == molsp(j)) then
			logic = .true.
			exit
		endif
	enddo

	if(.not. logic) then
		write(*,'("  Invalid molecular name in line: ",a)') trim(line)
		write(*,'("  Supported molecules: H2O, CO2, O3, N2O, CO, CH4, O2, OH, N2")')
		stop
	endif

	do j = 1, n_sld_group_species
		if(species == sld_group_species(j)) then
			write(*,'("  Duplicate SLD grouping rule for molecule: ",a4)') species
			stop
		endif
	enddo

	limits_string = adjustl(line(colon_position+1:))

	if(len_trim(limits_string) == 0) then
		write(*,'("  No intensity limits specified for molecule: ",a4)') species
		stop
	endif

	number_of_limits = 1

	do j = 1, len_trim(limits_string)
		if(limits_string(j:j) == ',') then
			number_of_limits = number_of_limits+1
		endif
	enddo

	allocate(parsed_limits(number_of_limits),stat=ios)

	if(ios /= 0) then
		write(*,'("  Memory allocation error for molecule: ",a4)') species
		stop
	endif

	do j = 1, number_of_limits

		comma_position = index(limits_string,',')

		if(comma_position > 0) then
			if(comma_position == 1) then
				write(*,'("  Empty intensity limit for molecule: ",a4)') species
				stop
			endif
			token = adjustl(limits_string(1:comma_position-1))
			limits_string = adjustl(limits_string(comma_position+1:))
		else
			token = adjustl(limits_string)
			limits_string = ' '
		endif

		if(len_trim(token) == 0) then
			write(*,'("  Empty intensity limit for molecule: ",a4)') species
			stop
		endif

		read(token,*,iostat=ios) parsed_limits(j)

		if(ios /= 0) then
			write(*,'("  Invalid intensity limit for molecule ", a4,": ",a)') species, trim(token)
			stop
		endif

		if(parsed_limits(j) <= 0.d0) then
			write(*,'("  Intensity limit for molecule ",a4," must be greater than zero: ",g15.7)') species, parsed_limits(j)
			stop
		endif

	enddo

	do j = 2, number_of_limits
		if(parsed_limits(j) <= parsed_limits(j-1)) then
			write(*,'("  Intensity limits for molecule ",a4," must be strictly increasing")') species
			stop
		endif
	enddo

	n_sld_group_species = n_sld_group_species+1

	sld_group_species(n_sld_group_species) = species
	sld_group_first(n_sld_group_species) = n_sld_group_limits+1
	sld_group_count(n_sld_group_species) = number_of_limits

	old_number_limits = n_sld_group_limits

	allocate(temporary_limits(old_number_limits+number_of_limits))

	if(old_number_limits > 0) then
		temporary_limits(1:old_number_limits) = sld_group_limits(1:old_number_limits)
	endif

	temporary_limits(old_number_limits+1:old_number_limits+number_of_limits) = parsed_limits(1:number_of_limits)

	if(allocated(sld_group_limits)) then
		deallocate(sld_group_limits)
	endif

	allocate(sld_group_limits(old_number_limits+number_of_limits))

	sld_group_limits(1:old_number_limits+number_of_limits) = temporary_limits(1:old_number_limits+number_of_limits)

	n_sld_group_limits = old_number_limits+number_of_limits

	deallocate(temporary_limits)
	deallocate(parsed_limits)

 enddo

 if(n_sld_group_species > 0) then
	do i = 1, n_sld_group_species
		if(i == 1) write(*,'("  SLD intensity groups:   | ",a4,": ")',	advance='no') sld_group_species(i)
		if(i > 1) write(*,'("                          | ",a4,": ")',	advance='no') sld_group_species(i)
		do j = 1, sld_group_count(i)
			if(j > 1) then
				write(*,'(", ")',advance='no')
			endif
			write(*,'(f4.2)',advance='no') sld_group_limits(sld_group_first(i)+j-1)
		enddo
		write(*,*)
	enddo
	write(*,*)
 endif

end subroutine READ_SLD_GROUP_LINE