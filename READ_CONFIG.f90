subroutine READ_CONFIG ! Initialize input parameters
 use SLD; use OBSERVED; use TELLURICA; use LINES
 implicit none
 integer                   :: i, k, ios
 character(20)             :: string
 character(80)             :: string1, string2

 open(1,file='SLD.config',status='old',iostat=ios) 
 if(ios /= 0) stop 'Configuration file not exist'

 read(1,'((a))') synth_path ! path to database
 k = index(synth_path,' '); synth_path = trim(adjustl(synth_path(1:k-1)))

 read(1,'((a))') observed_names
 k = index(observed_names,' '); observed_names = trim(adjustl(observed_names(1:k-1)))
 write(*,'("  File name spectra:      | ",a)')trim(adjustl(observed_names))

 read(1,'((a))') file_relative_humidity
 k = index(file_relative_humidity,' '); file_relative_humidity = trim(adjustl(file_relative_humidity(1:k-1)))
 write(*,'("  File relative humidity: | ",a)')trim(adjustl(file_relative_humidity))
 write(*,'("                          | ")')

 read(1,*) corr_mode, corr_target 
 if(corr_mode /= 0 .and. corr_mode /= 1 .and. corr_mode /= 2) stop 'Line correction parameters are set incorrectly'
 if(corr_target /= 0 .and. corr_target /= 1 .and. corr_target /= 2) stop 'Line correction parameters are set incorrectly'
 write(*,'("  Correct lines:          | ",i1,2x,i1)') corr_mode, corr_target
 write(*,'("                          | ")')

 read(1,*,iostat = ios) wfirst, wlast
 if(ios /= 0) stop 'Error reading wavelength region limits'
 write(string1,"(f9.2)")wfirst; 	write(string2,"(f9.2)")wlast
 string1 = trim(adjustl(string1))//" A -- "//trim(adjustl(string2))//" A"
 write(*,'("  Selected region:        | ",a)')trim(adjustl(string1))
 write(*,'("                          | ")')
 string1 = ''; string2 = ''
 
 read(1,*) vsini_user
 if(vsini_user < 0 .or. vsini_user > 500) stop 'Error: The Vsini must be greater than zero'
 write(string1,"(f6.2)") vsini_user
 write(*,'("  Vsini:                  | ",a)') trim(adjustl(string1))
 string1 = ''

 read(1,*) Resolution
 if(Resolution <= 0 .or. Resolution > 1000000) stop 'Error: The Resolution must be greater than zero'
 write(string1,"(i6)") Resolution/1000
 write(*,'("  Resolution:             | ",a,"k")') trim(adjustl(string1))
 write(*,'("                          | ")')
 string1 = ''

 if(height0 < 0.e0 .or. height0 > 1.2E+05) stop 'Error: Altitude must be between 0 km and 120 km'
 write(string1,"(f5.2)") height0/1.d3
 write(*,'("  Observatory height:     | ",a," km")') trim(adjustl(string1)) 
 string1 = ''

 if(nmodel < 1 .or. nmodel > 3) stop 'Error: Atmosphere Model Number must be between 1 and 3'
 if(nmodel == 1)string = 'Equatorial diurnal'; if(nmodel == 2)string = 'Mean latitude night'; if(nmodel == 3)string = 'Polar winter'
 write(*,'("  Atmosphere Model:       | ",a)') trim(adjustl(string))

 if(astep <= 0 .or. astep > 120000) stop 'Error: The altitude step in the atmosphere must be greater than zero'
 write(string1,"(f5.2)") real(astep)/1.d3
 write(*,'("  Altitude step:          | ",a," km")') trim(adjustl(string1)) 
 write(*,'("                          | ")')
 string1 = ''

 read(1,*) SLD_npoint
 if(SLD_npoint <= 10) stop 'Not SLD_npoint <= 10'
 if(mod(SLD_npoint,2) == 0) SLD_npoint = SLD_npoint+1
 write(string1,"(i3)") SLD_npoint
 write(*,'("  Points SLD-profile:     | ",a)') trim(adjustl(string1))
 write(*,'("                          | ")')
 string1 = ''

 read(1,'((a))') file_database_star
 k = index(file_database_star,' ') 
 file_info_database_star = trim(adjustl(file_database_star(1:k-1)))//'.dir'
 file_database_star = trim(adjustl(file_database_star(1:k-1)))//'.dbs'
 write(*,'("  File database star line:| ",a,/)')trim(adjustl(file_database_star))
 
 read(1,*,iostat = ios) ngroup
 if(ios /= 0) then
	ngroup = 0
 else
	allocate(config_line(ngroup))
	config_line = ''
	do i = 1,ngroup
		read(1,'(a)',iostat = ios) config_line(i)
		if(ios /= 0) then
			stop '   Error reading SLD grouping configuration'
		endif
	enddo
 endif

 close(1)

end subroutine READ_CONFIG
