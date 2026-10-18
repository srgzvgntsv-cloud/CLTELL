!*******************************************************************************
!*************    Hydrogen Lines Absorption    *********************************
!*******************************************************************************
subroutine Hlinop( wave ) ! Calculate  Hydrogen lines opacity
 use radiation_field
 use atmosphere_model
 use temporary_arrays
 use HlineTab
 implicit real(8) (a-h,o-z)
 real(8) Prof,nmax
 real(8), allocatable, save :: wij(:,:)
 
 integer :: itemp1 = 0
  if (itemp1 == 0 ) then ! Independent upon wavelengths
    itemp1 = 1
    do j = 1, nrhox
      Mlast(j) = 1100.d0 / xne(j)**(2.d0/15.d0) ! last visible line
    enddo
    do n = 1, nHydSeries
      do j = 1, nrhox
        ep      = 13.595d0 - 13.595d0 / (n*n)
        g       = 2.d0 * (n*n)
        xHh(n,j) = exp( - ep / Tkev(j) ) * g * xDens(j,1,1) / rho(j)
      enddo
    enddo

  nlevels=nHydSeries+nHydLines
  if(.not.allocated(wij))allocate (wij(nlevels,Nrhox))

  do i = 1, nlevels
   do j = 1, Nrhox
    wij(i,j) = occupationProbability(1d0,0d0,dble(i),0d0,t(j),xNe(j),xDens2(J,1,1)*0.d0)
   end do
  end do 
 endif

   Ahline = 0.d0
  nser = int( hydrogenN(freq) )
  if ( nser == 0 .or. nser > nHydSeries ) return
  if ( nser == 1 .and. freq < 2.d15 ) return
 ! if ( nser == 2 .and. freq < 4.44d14 ) return

  NHlines = Nhl( nser )
  term  = (109677.58375d0*2.997924580d10/dble(nser*nser)-freq)
 if ( term > 0d0 ) then
  mfreq = nint( sqrt(109677.58375d0*2.997924580d10/term) )
 else
! patch to avoid negative mfreq at the hydrogen series limits
  mfreq = 1000000
 end if 
  m1    = mfreq-nser
  m2    = m1 + 1
  m1    = max(m1,nser + 1)
  
if(m1 > maxval(mlast))return
  do j = 1, nrhox
!   if ( m1 > Mlast(j)) then
!    Ahline(j) =COULX(nser,3.28806d15/dble(nser*nser),1.d0)*(1.d0-ehvkT(j))*xHh(nser,j)
!    Sline(j) = Bnu(j)
!	cycle
!  endif
   line1=min(nser,4); line2 = line1
   if ( nser > 2 .and. nser <= 4 )line1 = line1 - 1
   do nlo=line1,line2 
      wi = wij(nlo,j)
     do line = 1, nHl(nlo)
       nup = nlo + line
       wj  =  wij(nup,j)
        del = abs( wave - wAir(Line,nlo) )
        vnm = 2.99792458d18 / wAir(Line,nlo)
        call VCS(nlo,del,Prof,Line,j)    ! Vidal, Cuper, Smith
        if ( nser == 2 .and. del > 0.2d0 .and. Line <= 4 ) then
         reson   = ( wAir(Line,nlo)**2 ) * 1.71d-27 * xDens(j,1,1) * 2.d0
         a       = radamp(Line) + reson
         pr2     = a / 3.1416d0 / ( (del**2) + a*a )
         Prof = ( Prof + pr2 ) * fHyd(Line,nlo)
        else
         Prof = Prof * fHyd(Line,nlo)
        endif
       stiml     = 1.d0 - exp( ( - 4.7993d-11 * vnm ) / T(j) )
       Prof   = Prof * stiml * xHh(nlo,j) / (vnm**2) * 7.957d16
       Ahline(j) = Ahline(j) + Prof*wj
	 enddo  
   enddo
    Sline(j) = Bnu(j)
  enddo

end subroutine Hlinop
!*******************************************************************************
subroutine VCS( ns, dw, Prof, n, idepth )
! Get Stark profile points for Hydrogen spectrum
! dw(A); n=1-Alpha,2-Beta etc.;ns=1 -Lyman...
!-----------------------------------------------
! useful only in following ranges:
! 2500 < T < 160,000; 10^10 < Ne < 10^18; 10^-3(-5) < delta(alpha) < 10^3
 use atmosphere_model
 use temporary_arrays
 use HlineTab
 use Lemkesdat
 implicit real(8) (a-h,o-z)
 real(8), parameter :: T0   = 3.39794d0
 real(8), parameter :: aNe0 = 10.d0
 real(8), parameter :: dTt   = 0.30103d0 ! delta[log(T)]
 real(8), parameter :: dNe  = 0.5d0     ! delta[log(Ne)]
 real(8), parameter :: dA   = 0.2d0     ! delta[log{dalta(alpha)}]
 real(8) Fo,Prof,dw
 integer mne(nHydSeries,nHydLines)
 save    Fo,mne
 integer ::  ns, n

 Fo    = 1.25d-9 * (Xne(idepth)**0.666666666667d0)
 if ( Told(idepth) /= T(idepth) ) then   ! Temperature and electron density interpolation
   Told(idepth) = T(idepth)
     aT    = log10(T(idepth))      ! Interpolation coefficients:
     bT    = (aT-T0)/dTt+1.d0
     iT    = bT
     iT    = max(min(iT,6),1)
     wT    = bT - iT
     aNe   = log10(Xne(idepth))
     aNe   = max(10.d0,aNe)
     bNe   = (aNe-aNe0)/dNe+1.d0
     iNe   = bNe
     iNe   = max(min(iNe,16),1)
     wNe   = bNe - iNe
     do ls = 1, nHydSeries
       nline = Nhl( ls )
       do line = 1, nline
         Npr = mp(ls,line)
         do k = 1, Npr	 	
	       ProfIn(ls,k,line,idepth) = (1.d0-wNe)*(1.d0-wT)*sVCS(ls,iT,iNe,k,line)+ &
                                 (1.d0-wNe)*wT*sVCS(ls,iT+1,iNe,k,line)+      &
                                 wNe*(1.d0-wT)*sVCS(ls,iT,iNe+1,k,line)+      &
                                 wNe*wT*sVCS(ls,iT+1,iNe+1,k,line)
	     end do    ! Alpha profile
       end do 	   ! Lines
     end do
 endif
! Alpha interpolation:
   if ( dw == 0.d0 ) then
     Prof = 10**(ProfIn(ns,1,n,idepth)) / Fo
   else
     aAlp = log10(abs(dw) / Fo)
     bAlp =(aAlp-a0(ns,n)) / dA + 1.d0
     iAlp = bAlp
     if ( iAlp <= 1 ) then
       Prof = 10**(ProfIn(ns,1,n,idepth)) / Fo
     else if ( iAlp > (mp(ns,n)-1) ) then
       Prof = 10**(ProfIn(ns,mp(ns,n),n,idepth) +  &
                 2.5d0 * (a0(ns,n) + (mp(ns,n) - 1) * dA - aAlp)) / Fo
     else
       wAlp = bAlp - iAlp
       Prof = 10**((1.d0-wAlp)*ProfIn(ns,iAlp,n,idepth) +   &
                 wAlp*ProfIn(ns,iAlp+1,n,idepth))/Fo
     end if
   end if
end subroutine VCS
!*******************************************************************************
