module OBSERVED
 implicit none
! Paths / file names
 character(160) :: observed_names
 character(160) :: synth_path
 character(160) :: file_relative_humidity
 character(160) :: file_database_star
 character(160) :: file_info_database_star
 character(80), allocatable :: config_line(:)
! Spectral regions / sampling
 integer :: Resolution
 integer :: nobs
 real(8)              :: wfirst, wlast
 real(8), allocatable :: wobs(:), robs(:), continuum(:)
! Global RV / star init params
 real(8) :: RV_tell
 real(8) :: RV_star
 real(8) :: vsini_user
! Global work flag
 integer :: if_stop
 integer, allocatable  :: bi(:)
end module OBSERVED


module LINES
 implicit none
! Convolution workspace
 integer :: np
 real(8) :: vp0, dvp
 real(8), allocatable :: vp(:), rp(:), r3p(:)
! HITRAN / generic line lists 

 integer :: nlin_tell
 integer, allocatable :: index_wavelenght_tell(:)
 integer, allocatable :: npoint_tell(:)
 real(8), allocatable :: wave_center_tell(:)
 real(8), allocatable :: line_rv_tell(:)
 character(4), allocatable :: line_species_tell(:)

! Star line database index
 integer :: nlin_star
 real(8), allocatable :: wave_center_star(:)
 real(8), allocatable :: line_rv_star(:)
 character(4), allocatable :: line_species_star(:)
 integer, allocatable :: npoint_star(:)
 integer, allocatable :: nrec_star(:)
 integer, allocatable :: index_wavelenght_star(:)

 integer :: irecl, nlines, line_active
 integer, allocatable :: line_molid(:), line_family(:), line_sld(:)
 real(8), allocatable :: line_corr(:), line_icorr(:)
 integer, allocatable :: index_wavelenght(:)
 integer, allocatable :: npoint(:), npoint_new(:), line_left(:), line_right(:)
 real(8), allocatable :: wave_center(:)
 real(8), allocatable :: line_rv(:)
 character(4), allocatable :: line_species(:)
 real(8), allocatable :: line_streng(:), line_intensity(:)
 real(8), allocatable :: lrv(:)

 integer :: corr_mode = 0
 integer :: corr_target = 0
 integer :: corr_vel = 0
 integer :: ncor = 0
 integer :: iterations
 real(8) :: vmin, vmax, pmin, pmax

 integer, allocatable :: corr_id(:)
 real(8), allocatable :: corr_fac(:)
 character(4), allocatable :: corr_label(:)
 integer, allocatable :: corr_family(:)
 interface
	subroutine LINE_ABSORBTION(mode,mode2,Jacobian)
		real(8), optional :: Jacobian(:,:)
		character(4), optional :: mode, mode2
	end subroutine LINE_ABSORBTION
 end interface

end module LINES


module SLD
 implicit none
 integer              :: SLD_npoint, nsld
 real(8)              :: step_star, step_tell
 real(8), allocatable :: SLD_velocity(:,:)
 real(8), allocatable :: SLD_intensity(:,:)
 real(8)              :: regpar(2)
 integer              :: ngroup = 0
 integer              :: n_sld_group_species = 0
 integer              :: n_sld_group_limits = 0
 character(4), allocatable :: sld_group_species(:)
 integer, allocatable :: sld_group_first(:)
 integer, allocatable :: sld_group_count(:)
 real(8), allocatable :: sld_group_limits(:)
end module SLD


module MODELS
 implicit none
 real(8), allocatable :: model(:), model_tmp(:)
end module MODELS

module TELLURICA
 implicit none
! Atmospheric layering
 integer :: nlayer
 integer :: nmodel
 real(8) :: astep ! meters
 real(8) :: height0 ! meters
 real(8), allocatable :: altitude(:), teff(:), press(:)
! Abundances
 real(8), allocatable :: abudance(:,:) ! (species, layer) VMR or similar
 real(8), allocatable :: abudance_mol(:,:) ! (species, layer) column or number-related quantity
 real(8), allocatable :: relative_humidity(:)
! Line parameters (per line)
 integer, allocatable :: imol(:)
 real(8), allocatable :: abscoef(:)
 real(8), allocatable :: gamma_air(:), gamma_self(:)
 real(8), allocatable :: energyi(:), n_air(:), d_air(:)
! Line broadening / opacities per (line, layer)
 real(8), allocatable :: dopcoef(:,:), lorcoef(:,:), abscoef_layer(:,:)
! Refraction / scattering
 real(8), allocatable :: refind(:,:) ! (wavelength index, layer)
 real(8), allocatable :: rayscat(:,:) ! (wavelength index, layer)
 real(8), allocatable :: rayscat_cont_tau(:)
 logical              :: abu_fit
 real(8)              :: scale_mol(9) = 1.d0
end module TELLURICA

module INSTRUMENT
 implicit none
 integer :: nip
 real(8), allocatable :: inst_prof(:)
 interface
	subroutine INSTRUMENTAL(mode,Jacobian2)
		real(8), optional :: Jacobian2(:,:)
		character(4), optional :: mode
	end subroutine INSTRUMENTAL
 end interface
end module INSTRUMENT


subroutine INIT_CONFIG
 use OBSERVED; use LINES; use SLD; use TELLURICA
 implicit none 
 integer(4)      :: ios

 open(1,file='config.ini',status='old',iostat=ios)
   if(ios /= 0) then
		corr_vel = 1
		vmin = -3.d0; vmax = 3.d0
		pmin = 1.d-3; pmax = 1.d3
		iterations = 2
		step_star = 1.d0
		step_tell = 1.d-1
		regpar(1) = 0.01; regpar(2) = 0.0
		height0 = 0
		nmodel = 2
		astep = 2000
   else
		read(1,*) iterations
		read(1,*) corr_vel
		read(1,*) vmin
		read(1,*) vmax
		read(1,*) pmin
		read(1,*) pmax
		read(1,*) step_star
		read(1,*) step_tell
		read(1,*) regpar
		read(1,*) height0
		read(1,*) nmodel
		read(1,*) astep
   endif
   close(1)

end subroutine INIT_CONFIG