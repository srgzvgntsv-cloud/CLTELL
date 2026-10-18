subroutine MERGING_LIST_LIN
 use LINES; use SLD; use MOLECULE_ISOTOP; use TELLURICA
 implicit none
 integer  :: il, i, j, k
 type(isotop) :: mol

 nlines = nlin_tell+nlin_star

 if(allocated(npoint)) deallocate(npoint)
 allocate(npoint(nlines))

 if(allocated(line_species)) deallocate(line_species)
 allocate(line_species(nlines))

 if(allocated(index_wavelenght)) deallocate(index_wavelenght)
 allocate(index_wavelenght(nlines))

 if(allocated(wave_center)) deallocate(wave_center)
 allocate(wave_center(nlines))

 if(allocated(line_rv)) deallocate(line_rv)
 allocate(line_rv(nlines))

 if(allocated(line_family)) deallocate(line_family)
 allocate(line_family(nlines))

 do il = 1, nlin_tell
	line_species(il) = line_species_tell(il)
	index_wavelenght(il) = index_wavelenght_tell(il)
	line_family(il) = 1
	wave_center(il) = wave_center_tell(il)
	line_rv(il) = line_rv_tell(il)
 enddo

 do il = 1, nlin_star
	line_species(il+nlin_tell) = line_species_star(il)
	index_wavelenght(il+nlin_tell) = index_wavelenght_star(il)
	line_family(il+nlin_tell) = 2
	wave_center(il+nlin_tell) = wave_center_star(il)
	line_rv(il+nlin_tell) = line_rv_star(il)
	npoint(il+nlin_tell) = npoint_star(il)
 enddo

! deallocate(line_species_tell,index_wavelenght_tell,wave_center_tell,line_rv_tell)
! deallocate(line_species_star,index_wavelenght_star,wave_center_star,line_rv_star)


 if(allocated(line_corr)) deallocate(line_corr)
 allocate(line_corr(nlines))
 line_corr = 1.d0

 if(allocated(line_icorr)) deallocate(line_icorr)
 allocate(line_icorr(nlines))
 line_icorr = 1.d0

 if(allocated(line_molid)) deallocate(line_molid)
 allocate(line_molid(nlines))
 line_molid = 0

 do i = 1,nlines
	if(line_family(i) == 1) then
		mol = get_id(imol(i))
		line_molid(i) = mol.id_molecule
	endif	
 enddo

end subroutine MERGING_LIST_LIN