subroutine DEALLOCATE_ARRAY
 use OBSERVED; use LINES; use SLD; use MODELS; use TELLURICA; use INSTRUMENT
 implicit none
! --------------------------------------------------------------------------------------------------
 if(allocated(wobs))                    deallocate(wobs)					! module OBSERVED
 if(allocated(robs))                    deallocate(robs)
! --------------------------------------------------------------------------------------------------
 np = 0; vp0 = 0; dvp = 0
 if(allocated(vp))                      deallocate(vp)		                 ! module LINES
 if(allocated(rp))                      deallocate(rp)
 if(allocated(r3p))                     deallocate(r3p)
 if(allocated(index_wavelenght_tell))   deallocate(index_wavelenght_tell)
 if(allocated(npoint_tell))             deallocate(npoint_tell)
 if(allocated(wave_center_tell))        deallocate(wave_center_tell)
 if(allocated(line_rv_tell))            deallocate(line_rv_tell)
 if(allocated(line_species_tell))       deallocate(line_species_tell)
 nlin_tell = 0; nlin_star = 0
 if(allocated(wave_center_star))        deallocate(wave_center_star)
 if(allocated(line_rv_star))            deallocate(line_rv_star)
 if(allocated(line_species_star))       deallocate(line_species_star)
 if(allocated(npoint_star))             deallocate(npoint_star)
 if(allocated(nrec_star))               deallocate(nrec_star)
 if(allocated(index_wavelenght_star))   deallocate(index_wavelenght_star)
 irecl = 0; nlines = 0
 if(allocated(line_molid))              deallocate(line_molid)
 if(allocated(line_family))             deallocate(line_family)
 if(allocated(line_sld))                deallocate(line_sld)
 if(allocated(line_corr))               deallocate(line_corr)
 if(allocated(line_icorr))              deallocate(line_icorr)
 if(allocated(lrv))                     deallocate(lrv)
 if(allocated(index_wavelenght))        deallocate(index_wavelenght)
 if(allocated(npoint))                  deallocate(npoint)
 if(allocated(npoint_new))              deallocate(npoint_new)
 if(allocated(wave_center))             deallocate(wave_center)
 if(allocated(line_rv))                 deallocate(line_rv)
 if(allocated(line_species))            deallocate(line_species)
 corr_mode = 0; corr_target = 0; ncor = 0
 if(allocated(corr_id))                 deallocate(corr_id)
 if(allocated(corr_fac))                deallocate(corr_fac)
 if(allocated(corr_label))              deallocate(corr_label)
 if(allocated(corr_family))             deallocate(corr_family)
! --------------------------------------------------------------------------------------------------
 SLD_npoint = 0; nsld = 0
 if(allocated(SLD_intensity))           deallocate(SLD_intensity)            ! module SLD
 if(allocated(SLD_velocity))            deallocate(SLD_velocity)
! --------------------------------------------------------------------------------------------------
 if(allocated(model))                   deallocate(model)                    ! module MODELS
 if(allocated(model_tmp))               deallocate(model_tmp)
! --------------------------------------------------------------------------------------------------
 nlayer = 0; astep = 0.d0; height0 = 0.d0
 if(allocated(altitude))                deallocate(altitude)                 ! module TELLURICA
 if(allocated(teff))                    deallocate(teff)
 if(allocated(press))                   deallocate(press)
 if(allocated(abudance_mol))            deallocate(abudance_mol)
 if(allocated(abudance))                deallocate(abudance)
 if(allocated(relative_humidity))       deallocate(relative_humidity)
 if(allocated(imol))                    deallocate(imol)                
 if(allocated(abscoef))                 deallocate(abscoef)
 if(allocated(gamma_air))               deallocate(gamma_air)
 if(allocated(gamma_self))              deallocate(gamma_self)
 if(allocated(energyi))                 deallocate(energyi)
 if(allocated(n_air))                   deallocate(n_air)
 if(allocated(d_air))                   deallocate(d_air)
 if(allocated(dopcoef))                 deallocate(dopcoef)
 if(allocated(lorcoef))                 deallocate(lorcoef)
 if(allocated(abscoef_layer))           deallocate(abscoef_layer)
 abu_fit = .false.
! --------------------------------------------------------------------------------------------------
 nip = 0
 if(allocated(inst_prof))               deallocate(inst_prof)                 ! module INSTRUMENT
! --------------------------------------------------------------------------------------------------
 close(11,status='delete')

end subroutine DEALLOCATE_ARRAY




