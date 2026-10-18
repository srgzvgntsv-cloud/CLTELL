program SynthV ! Calculate synthetic spectrum for given stellar atmosphere model with (optional) database for SLD
use synth_data
use radiation_field
use dfport
character(8)time_start,time_finish

 write(*,*)'SynthVdb 04.05.19'
 time_start=clock() 

 call SynthV_input  ! Read input data and model atmosphere. Open databases and output files

 call number_density ! Calculate number densities of all species

 call calculate_tau5000 ! Calculate opacity depths at 5000 A

 call continuum_opacity ! Calculate opacity in continuum for given wavelengths interval

 call lines_opacity     ! Preselection of spectral lines and calculation lines opacity

 call calculate_spectrum
 if(if_database /= 0) then
	write(*,*)'Start form DataBase'
	call form_database
 endif

 time_finish=clock()
 write(*,*)'Start at ',time_start
 write(*,*)'Finish at ',time_finish
 read(time_start,'(i2,1x,i2,1x,i2)')ih1,im1,is1
 read(time_finish,'(i2,1x,i2,1x,i2)')ih2,im2,is2
 t1=ih1*3600+im1*60+is1; t2=ih2*3600+im2*60+is2
 write(*,*)'Duration: ',t2-t1,' sec'

 close(200,status='delete')
 call arrays_deallocation
 call deallocate_continuum_arrays
end