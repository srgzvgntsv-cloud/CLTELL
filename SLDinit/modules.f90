! All common data for subroutines are placed in different modules

module synth_data ! Information about synthetic spectrum
real(8) wavelength_first,wavelength_last,wavelength_step
integer number_of_wavelengths
integer atomic_lines, molecular_lines ! Switch: will be given type of lines absorbtion included in calculations
integer if_ev  ! In VALD list: 1 -ev, o - cm-1
integer first_record_atomic,first_record_molecular ! record numbers in lines lists for first line
real(8) V_micro ! microturbulent velocity
real(8),allocatable :: Vturb(:),Vrad(:)  ! Microturbulent velocity and radial velocity vs. depths
integer Vrad_switch                      ! If this value =1, Vrad vs. depths implemented
real(8),allocatable :: tau5000(:),kappa5000(:),logtau5000(:)        ! Optical depths at 5000 A, corresponded opacity and log10 from optical depths
integer number_stratified ! Number of chemical elements with stratified abundace
real(8),allocatable :: xab_Strat(:,:) ! Array with stratified abundencies
integer atomic_numbers_stratified(99) ! Atomic numbers of chemical elements with stratified abundance
integer number_of_lines
character(3) output_format
data atomic_lines/0/, molecular_lines/0/ ! Absence as default
data number_of_lines/0/                  ! Absence as default
integer if_database,nlines_db ! switch for calculating lines profiles database: 0 or /= 0, numbers of non Hydrogen lines in it
data if_database/0/

character(100) proffilename,  proffilename_bin,  proffilename_inf
real(8), allocatable :: synth_wl(:), synth_R(:)
integer itest

end module synth_data


module atmosphere_model
real(8) Teff,grav,glog,scale_log !Effective temperature, surface gravity and it's log, metallicity
character(74) Title ! Model atmosphere title
integer Nrhox ! Nuber of depths
real(8),allocatable :: T(:),Tkev(:),Tk(:),hkT(:),Tlog(:) ! Temperature arrays (Tk really kT)
real(8),allocatable :: P(:),Rhox(:),xNe(:),xNatom(:),rho(:) ! Hressure,column mass,Ne, Ntotal,density
real(8),allocatable :: xDens(:,:,:),xDens2(:,:,:)   ! Number densities of atoms and ions: xDens are N/partition function
real(8),allocatable :: xnmol(:,:),x_mol(:,:) ! same for molecules
real(8),allocatable :: abRoss(:),tauros(:) ! Rosseland means and optical depths
real(8) xAbund(99),wtmole,xScale ! Accepted abundancies,molecular weight and metallicity scale
real(8) model_abundance(99) 
end module atmosphere_model

module radiation_field
real(8) Freq,Freqlg  ! Frequency and it log 
real(8),allocatable :: EhvkT(:),Stim(:),Bnu(:) !Ehv/kT, 1-Ehv/kT, Plank function
real(8),allocatable :: Taunu(:),Snu(:),Hnu(:),Jnu(:),Jmins(:),Knu(:) ! Optical depth,source fuction and integral characteristics of radiation field 
real(8),allocatable :: Acont(:),Aline(:),AHline(:),Abtot(:),sigmaC(:) ! Absorbtion and scattering coefficients
real(8),allocatable :: Aline_mu(:,:) ! A special array for Vrad /= 0 (line absorbtion at fixed cos(theta))
real(8),allocatable :: Scont(:),Sline(:),Alpha(:) !Source functions and albedo
real(8) Inu(7)       ! Specific intensities at 7 fixed cos(theta)
real(8) :: mu_fix(7) = (/0.9636d0,0.8864d0,0.8018d0,0.7071d0,0.5976d0,0.4629d0,0.2673d0/) ! Fixed cos(theta) values for specific intensities
real(8),allocatable :: H_cont(:),S_cont(:),J_cont(:) ! Continuum radiation field
real(8) effective_depth_tau(1),effective_depth_number(1) ! effective depth data 
end module radiation_field

module lines_data ! Data and arrays for spectral lines absorbtion: data for GIVEN line
real(8),allocatable :: al_int(:),Adump(:),Dopdump(:) ! Integral absorbtion and lines dumping parameters
real(8),allocatable :: aline_buffer(:)     ! Given line absorbtion vs. depths 
real(8) central_wavelength,gf,Elow,code,gamma_Rad,gamma_Stark,gamma_WdW ! wavelength of line centre,gf,excitation energy,ion code,dumping constants
integer atomic_number, ionization
character(200) VALD_string  ! VALD format of given line
character(4) line_identifier
end module lines_data

module temporary_arrays
! The set of arrays, used mainly in KAPP and other Kurucz's subroutines
real(8),allocatable :: frac(:,:),PFPLUS(:),PFMIN(:),Told(:)
real(8),allocatable :: Aa(:),Bb(:),Cc(:),xh(:),xjs(:),CTWO(:),B2CT(:),B2CT1(:)
real(8),allocatable :: THETA(:),FFTHETA(:)
real(8),allocatable :: BOLT(:,:),EXLIM(:),BOLTN(:,:),BOLTEX(:),FREET(:),w(:),dt(:)
real(8),allocatable :: C1218(:),C1420(:),xHmin(:),C1240(:),C1444(:),bolt2(:,:)
real(8),allocatable :: C1130(:),C1020(:),C1169(:),ahot(:),xnz(:,:),answer(:,:)
real(8),allocatable :: E1(:,:),E2(:,:),Tx(:),vt2(:)
real(8),allocatable :: xHh(:,:),Mlast(:),ProfIn(:,:,:,:)
integer(4),allocatable :: nt(:)
real(8),allocatable :: xnebuf(:)
 save
end module temporary_arrays
