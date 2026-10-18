program CLTELL
 use LINES
 use TELL_ABU_FIT_MOL
 implicit none
 character(8) :: time_start

 write(*,'(a,/)')'CLTELL: version 3.0, 02.05.2026'

 call START_TIMER(time_start)         ! Start of the program running time countdown                 

 call INIT_CONFIG

 call READ_CONFIG                     ! Reading a configuration file

 call READ_SLD_GROUP_LINE

 call LOAD_OBSERVED                   ! Loading the observed spectrum into memory

 call GET_LAYERS                      ! Stratification of the atmosphere into a given number of layers

 call READ_MODEL_EARLTH_ATMOSPHERE    ! Reading data from the earth's atmosphere model

 call READ_RELATIVE_HUMIDITY          ! Reading data from a relative humidity file

! call REFRACTIVE_INDEX                ! Calculation of refractive index

! --------------------------------------------------------------------------------------------------------

 call READ_HITRAN

 call LOAD_STAR_LINE

 call MERGING_LIST_LIN

! --------------------------------------------------------------------------------------------------------

 if(abu_fit == .true.) then

	call GET_ABSCOEF(0)                     ! Calculation of the absorption coefficient in the line for a given temperature and pressure in the layer

	call GET_DOPCOEF(0)                     ! Calculating HWHM for Doppler Profile

	call GET_LORCOEF(0)                     ! Calculating HWHM for Lorenz Profile

	call TELL_MODEL_RV 

	call TELL_ABU_MOL_CORRECT(40.0d0)
 
 endif

 ! --------------------------------------------------------------------------------------------------------
 
 call GET_ABSCOEF(0)                     ! Calculation of the absorption coefficient in the line for a given temperature and pressure in the layer

 call GET_DOPCOEF(0)                     ! Calculating HWHM for Doppler Profile

 call GET_LORCOEF(0)                     ! Calculating HWHM for Lorenz Profile

! call GET_WATER_CONT_COEF(0)

! call GET_RAYSCAT(0)
 
 call FBASE_PLINE                     ! Formation of a database with profiles of telluric lines

 call INIT_SLD               ! Removing weak lines from the telluric line profile database

 call FIND_RV                         ! Search for radial velocity for telluric lines

 ! ------------------------------------------------------------------------------------------------------------------------------

 call FORM_TELLURIC                   ! Correction and fitting of the terrestrial absorption spectrum to the observed one. Removal of telluric absorption from the spectrum of a star

 if(corr_vel == 1) call VLINE_CORRECT

 if(corr_mode == 2) then
	call ILINE_CORRECT
 endif

 if(corr_mode == 0 .or. corr_mode == 1) then
	call LINE_CORRECT(corr_mode, corr_target)
 endif

 call FORM_TELLURIC

 call FIND_INSTRUMENTAL               ! Search and refinement of the instrumental profile

 ! --------------------------------------------------------------------------------------------------------

 call LINE_ABSORBTION('tell','tell')

 close(7, status='delete')            ! Deleting the file with the telluric line profile database

 call DEALLOCATE_ARRAY

 call SHOW_DURATION(time_start)       ! Finish of the program running time countdown

 stop

end program CLTELL



