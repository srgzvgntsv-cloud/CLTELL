subroutine LOAD_STAR_LINE
 use LINES
 use OBSERVED
 use SLD
 implicit none

 integer :: ios, n, i, reg, nprof_i
 integer :: n1, n2, database_version
 integer, external :: INDWAV
 real(8) :: w, r, d1, d2
 character(4) :: e
 character(8) :: database_magic
 integer, allocatable :: ifline(:)

 open(25, file=trim(adjustl(synth_path))//trim(adjustl(file_info_database_star)), &
 form='unformatted', access='sequential', status='old', action='read', iostat=ios)
 if(ios /= 0) stop 'file with database information not exist'

 read(25, iostat=ios) database_magic, database_version, n1, n2
 if(ios /= 0) stop 'LOAD_STAR_LINE: cannot read database header'
 if(database_magic /= 'SPROFDB2') stop 'LOAD_STAR_LINE: invalid database format'
 if(database_version /= 2) stop 'LOAD_STAR_LINE: unsupported database version'

 open(35, file=trim(adjustl(synth_path))//trim(adjustl(file_database_star)), &
 form='unformatted', access='sequential', status='old', action='read', iostat=ios)
 if(ios /= 0) stop 'file with database not exist'

 n = n1 + n2
 if(n < 0) stop 'LOAD_STAR_LINE: invalid number of records'

 allocate(ifline(n))
 ifline = 0

 nlin_star = 0

 do i = 1, n
	read(25, iostat=ios) w, e, r, nprof_i, d1, d2
	if(ios /= 0) stop 'LOAD_STAR_LINE: cannot read line information'

	if(w < wfirst) cycle
	if(w > wlast) cycle

	if(INDWAV(w, 0.d0) == 0) cycle

	nlin_star = nlin_star + 1
	ifline(i) = 1
 enddo

 if(allocated(wave_center_star)) deallocate(wave_center_star)
 if(allocated(line_species_star)) deallocate(line_species_star)
 if(allocated(npoint_star)) deallocate(npoint_star)
 if(allocated(nrec_star)) deallocate(nrec_star)
 if(allocated(index_wavelenght_star)) deallocate(index_wavelenght_star)
 if(allocated(line_rv_star)) deallocate(line_rv_star)

 allocate(wave_center_star(nlin_star))
 allocate(line_species_star(nlin_star))
 allocate(npoint_star(nlin_star))
 allocate(nrec_star(nlin_star))
 allocate(index_wavelenght_star(nlin_star))
 allocate(line_rv_star(nlin_star))

 rewind(25)

 read(25, iostat=ios) database_magic, database_version, n1, n2
 if(ios /= 0) stop 'LOAD_STAR_LINE: cannot reread database header'

 nlin_star = 0

 do i = 1, n
	read(25, iostat=ios) w, e, r, nprof_i, d1, d2
	if(ios /= 0) stop 'LOAD_STAR_LINE: cannot reread line information'

	if(ifline(i) == 1) then
		nlin_star = nlin_star + 1
		wave_center_star(nlin_star) = w
		line_species_star(nlin_star) = e
		npoint_star(nlin_star) = nprof_i
		nrec_star(nlin_star) = i
		index_wavelenght_star(nlin_star) = INDWAV(w, 0.d0)
	endif
 enddo

 if(size(line_rv_star) > 0) line_rv_star = 0.d0

 deallocate(ifline)
 close(25)

end subroutine LOAD_STAR_LINE