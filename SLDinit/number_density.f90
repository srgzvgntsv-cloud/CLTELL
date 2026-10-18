subroutine number_density   ! Calculate the number densities
 use synth_data
 use atmosphere_model
 use chemical_elements
 use temporary_arrays
 implicit real(8) (a-h,o-z)
 integer(4) LOCZ(29)/1,3,6,10,14,18,22,27,33,39,45,51,57,63,69,75,81,86,91, &
                96,101,106,111,116,121,126,131,136,141/

 if(abund(1).lt.0.)abund(1)=10**(abund(1))  ! Abundancies:
 if(abund(2).lt.0.)abund(2)=10**(abund(2))
 xabund(1)=abund(1)
 xabund(2)=abund(2)
 do j=3,99
  if(abund(j).gt.0.)abund(j)=log10(abund(j))
  xabund(j)=10**(abund(j))
  abund(j)=log10(xabund(j))
 enddo
 wtmole=0.d0   ! Molecular weight:
  do j=1,99
    wtmole=wtmole+xabund(j)*atmass(j)
  enddo
 do j=1,nrhox             ! Temperature and pressure depended quantities
  Tk(J)=1.38054E-16*T(J)  ! kT
  hkT(J)=6.6256E-27/Tk(J) ! hkT
  Tkev(J)=8.6171E-5*T(J)  ! In ev
  Tlog(J)=LOG(T(J))
  xNatom(J)=P(J)/Tk(J)-xNe(J)  ! Total number of atoms
  rho(J)=xNatom(J)*wtmole*1.660E-24 ! Density
 enddo

 ! call INTEG(rhox,abross,tauros,nrhox,abross(1)*rhox(1)) ! Tau Ross  Change to tau5000!!!!

! Calculate  number densities of plasma's components:
 xdens=0.d0; xDens2=0.d0
  do j=1,Nrhox
   do k=1,99
    if(k.gt.28)Nions=3 ! For heavy elements max third spectra
    if(k.le.28)Nions=LOCZ(K+1)-LOCZ(K)
    if(k.eq.6.or.k.eq.7)Nions=6
    call PFSAHA(J,K,NIONS,11)
    E1=ANSWER
    call PFSAHA(j,k,nions,12)
    E2=ANSWER
     xx=xabund(k)
      if(number_stratified > 0)then ! If stratification are implemented
       do i=1,number_stratified
        if(atomic_numbers_stratified(i) == k)then
          xx=xab_Strat(i,j)
         exit
        endif
       enddo
      endif
      do  m=1,Nions
        C=1.d0
        IF(k.ne.1)C=1.49725d-2/rho(j)
         xDens(j,k,m)=E1(j,m)*xNatom(j)*xx*C
	     xDens2(j,k,m)=E2(j,m)*xNatom(j)*xx
      enddo
    enddo
  enddo
 close (1)

 if(molecular_lines > 0)then
  xnebuf=xne    ! Save the Ne from atmosphere model
  call nmolec  ! Molecular number densities
  xne=xnebuf    ! Restore the Ne from atmosphere model
 endif

end

SUBROUTINE PFSAHA(J,IZ,NION,MODE) ! Calculate number densities, number densities/partition functions for atmoic species
 use atmosphere_model
 use temporary_arrays
 use PfSahaTab

 implicit real(8) (a-h,o-z)
  DIMENSION F(6),PART(6),POTLO(6),LOCZ(29)
  REAL(8) IP(6)
  data LOCZ/1,3,6,10,14,18,22,27,33,39,45,51,57,63,69,75,81,86,91, &
                96,101,106,111,116,121,126,131,136,141/
!     MODE 1 RETURNS IONIZATION FRACTION /PARTITION FUNCTION
!     MODE 2 RETURNS IONIZATION FRACTION
!     MODE 3 RETURNS PARTITION FUNCTION
!     MODE 4 RETURNS NUMBER OF ELECTRONS PRODUCED
!     MODE 5 RETURNS ANSWER(ION)=PF   ANSWER(ION+7)=IP
!     MODE + 10 RETURN ALL IONS TO NION.   MODE ALONE RETURN NION ONLY.
  MODE1=MODE
   IF(MODE1 > 10)MODE1=MODE1-10
!     LOWERING OF THE IONIZATION POTENTIAL IN VOLTS FOR UNIT ZEFF
  CHARGE=XNE(J)*2.d0
  EXCESS=2.d0*XNE(J)-P(J)/TK(J)
!     ALLOWANCE FOR DOUBLY IONIZED HELIUM
   IF(EXCESS.GT.0.)CHARGE=CHARGE+2.d0*EXCESS
  DEBYE=SQRT(TK(J)/2.8965d-18/CHARGE)
  POTLOW=MIN(1.D0,1.44d-7/DEBYE)
  TV=TKEV(J)
   if(iz == 6)then
    n=354
	nions=6
   else
    if(iz == 7)then
	  n=360
	  nions=6
	else
	 if(iz <= 28)then     	
      n=LOCZ(iz)
	  nions=locz(iz+1)-n
     else
	  n=3*iz+54
      nions=3
     endif
	endif
   endif	
   NION2=MIN(NION + 2,NIONS)
   N=N-1
  DO ION=1,NION2
   Z=ION
   POTLO(ION)=POTLOW*Z
   N=N+1
   NNN100=NNN(6,N)/100
   IP(ION)=NNN100/1000.d0
   G=NNN(6,N)-NNN100*100
    IF(N == 1)then  ! Hydrogen
     PART(1)=2.d0
     PART(1)=PART(1)+8.d0*EXP(-10.196d0/TV)+18.d0*&
      EXP(-12.084d0/TV)+32.d0*EXP(-12.745d0/TV)+50.d0*&
      EXP(-13.051d0/TV)+72.d0*EXP(-13.217d0/TV)
     D1=13.595d0/6.5d0/6.5d0/TV
     D2=POTLO(1)/TV
    else
      T2000=IP(ION)*2000.d0/11.d0
      IT=MAX(1,MIN(9,INT(T(J)/T2000-0.5d0)))
      DTt=T(J)/T2000-float(IT)-0.5d0
      PMIN=1.d0
      I=(IT+1)/2
      K1=NNN(I,N)/100000
      K2=NNN(I,N)-K1*100000
      K3=K2/10
      KSCALE=K2-K3*10
      IF(MOD(IT,2) == 0)then
       P1=K3*SCALE(KSCALE)
       K1=NNN(I+1,N)/100000
       KSCALE=MOD(NNN(I+1,N),10)
       P2=K1*SCALE(KSCALE)
      else
       P1=K1*SCALE(KSCALE)
       P2=K3*SCALE(KSCALE)
        IF(DTt < 0..and.kscale < 1)then
         KP1=P1
         IF(KP1 == INT(P2+.5))PMIN=KP1
         endif
       endif
      PART(ION)=MAX(PMIN,P1+(P2-P1)*DTt)
      IF(G == 0.d0.OR.POTLO(ION) < 0.1 .OR.T(J) < T2000*4.)cycle
      IF(T(J) > (T2000*11.))TV=(T2000*11.d0)*8.6171d-5
      D1=0.1d0/TV
     endif
    D2=POTLO(ION)/TV
    PART(ION)=PART(ION)+G*EXP(-IP(ION)/TV)*(SQRT(13.595d0*Z*Z/TV/D2)**3*&
      (1.d0/3.d0+(1.d0-(.5d0+(1.d0/18.d0+D2/120.d0)*D2)*D2)*D2)&
      -SQRT(13.595d0*Z*Z/TV/D1)**3*(1.d0/3.d0+(1.d0-(0.5d0+&
     (1.d0/18.d0+D1/120.d0)*D1)*D1)*D1))
      TV=TKEV(J)
  enddo
  if(mode1 == 5)then
   ANSWER(7,1)=0.d0
    DO ION=1,NION
     ANSWER(ION,1)=PART(ION)
     ANSWER(ION+7,1)=IP(ION)+ANSWER(ION+6,1)
    enddo
   return
  else
   if(mode1 /=3)then
    N=N-NION2
    CF=2.d0*2.4148d15*T(J)*SQRT(T(J))/XNE(J)
     do ION=2,NION2
      N=N+1
! THE AMIN IS FOR ANY UNFORTUNATE WHO HAS A 360
      F(ION)=CF*PART(ION)/PART(ION-1)*&
       EXP(-MIN((IP(ION-1)-POTLO(ION-1))/TV,100.D0))
     enddo
    F(1)=1.d0
    L=NION2+1
     do ION=2,NION2
      L=L-1
      F(1)=1.+F(L)*F(1)
     enddo
    F(1)=1.d0/F(1)
      do ION=2,NION2
       F(ION)=F(ION-1)*F(ION)
	  enddo
     endif
  endif
  IF(MODE >10)then
   if(mode1 == 1)then
     do ION=1,NION
      ANSWER(J,ION)=F(ION)/PART(ION)
     enddo
   endif
   if(mode1 == 2)then
     do ION=1,NION
      ANSWER(J,ION)=F(ION)
     enddo
   endif
   if(mode1 == 3)then
     do ION=1,NION
      ANSWER(J,ION)=PART(ION)
     enddo
   endif
   if(mode1 == 4)then
    ANSWER(J,1)=0.
     do ION=2,NION2
      ANSWER(J,1)=ANSWER(J,1)+F(ION)*DFLOAT(ION-1)
     enddo
   endif
  else
    if(mode1 == 1)ANSWER(J,1)=F(NION)/PART(NION)
    if(mode1 == 2)ANSWER(J,1)=F(NION)
    if(mode1 == 3)ANSWER(J,1)=PART(NION)
  endif
END

subroutine nmolec ! Solution of equition of state - direct copy of Kurucz's code
 use mol_dat
 use atmosphere_model
 use chemical_elements
 use temporary_arrays

 implicit real(8) (a-h,o-z)
 DIMENSION PFP(13),PFM(13),EION(7)
 EQUIVALENCE (PFP(7),EION(1))
 DIMENSION EQUILJ(160)
 DIMENSION EQ(25),XN(25),XAB(25),IDTERM(25),DEQ(625)
!EQUIVALENCE (FRAC(1,1),DEQ(1))
 DIMENSION EQOLD(25)
! This part placed instead readmol subroutine:
 ifequa=0
 kloc=1
 locj(1)=1
 DO JMOL=1,nummol
   DO ii=1,8
    IF(code(jmol) >=XCODE(II))exit
   enddo
 X=code(jmol)
  DO I=II,8
   ID=X/XCODE(I)+.5
   X=X-ID*XCODE(I)
   IF(ID == 0)ID=100
   IFEQUA(ID)=1
   KCOMPS(KLOC)=ID
   KLOC=KLOC+1
  enddo
  ION=X*100.+.5
   IF(ION >= 1)then
    IFEQUA(100)=1
    IFEQUA(101)=1
     DO I=1,ION
      KCOMPS(KLOC)=101
      KLOC=KLOC+1
     enddo
   endif
   LOCJ(JMOL+1)=KLOC
 enddo
 NLOC=KLOC-1
 IEQUA=1
  DO I=1,100
  IF(IFEQUA(I) == 0)cycle
   IEQUA=IEQUA+1
   IFEQUA(I)=IEQUA
   IDEQUA(IEQUA)=I
  enddo
   NEQUA=IEQUA
   NEQUA1=NEQUA+1
   IFEQUA(101)=NEQUA1
   NEQNEQ=NEQUA**2
   DO KLOC=1,NLOC
    ID=KCOMPS(KLOC)
    KCOMPS(KLOC)=IFEQUA(ID)
   enddo
! This part is copy of nmolec subroutine
  NEQUA1=NEQUA+1
  NEQNEQ=NEQUA**2
  DO K=2,NEQUA
    ID=IDEQUA(K)
    IF(ID < 100)XAB(K)=MAX(XABUND(ID),1.D-20)
  enddo
    IF(ID == 100)XAB(NEQUA)=0.d0
    XNTOT=P(1)/TK(1)
    XN(1)=XNTOT/2.d0
    IF(T(1) < 4000.)XN(1)=XNTOT
    X=XN(1)/10.d0
   DO K=2,NEQUA
    XN(K)=X*XAB(K)
   enddo
   IF(ID == 100)XN(NEQUA)=X
   XNE(1)=X
    DO J=1,NRHOX
     XNTOT=P(J)/TK(J)
     IF(J /= 1)then
      RATIO=P(J)/P(J-1)
      XNE(J)=XNE(J-1)*RATIO
      DO K=1,NEQUA
       XN(K)=XN(K)*RATIO
      enddo
     endif
   DO JMOL=1,NUMMOL
    NCOMP=LOCJ(JMOL+1)-LOCJ(JMOL)
    IF(EQUIL(1,JMOL) /= 0.d0)then
      ION=(CODE(JMOL)-DFLOAT( INT(CODE(JMOL))))*100.+.5
      EQUILJ(JMOL)=0.d0
     IF(T(J) < 10000.)then
      EQUILJ(JMOL)=EXP(EQUIL(1,JMOL)/TKEV(J)-EQUIL(2,JMOL)+&
      (EQUIL(3,JMOL)+(-EQUIL(4,JMOL)+(EQUIL(5,JMOL)-EQUIL(6,JMOL)*&
      T(J))*T(J))*T(J))*T(J)-1.5*(DFLOAT(NCOMP-ION-ION-1))*TLOG(J))
      endif
     cycle
    endif
   IF(NCOMP <= 1)then
      EQUILJ(JMOL)=1.d0
      cycle
   endif
   ID=CODE(JMOL)
   ION=NCOMP-1
    CALL PFSAHA(J,ID,NCOMP,12)
    FRAC=ANSWER
    EQUILJ(JMOL)=FRAC(J,NCOMP)/FRAC(J,1)*XNE(J)**ION
   enddo
   DO K=1,NEQUA
    EQOLD(K)=0.d0
   enddo
!     SET UP 1ST ORDER EQUATIONS FOR THE CHANGE IN NUMBER DENSITY OF EACH ELEMENT.
 do
  deq=0.d0
    EQ(1)=-XNTOT
    K1=1
    KK=1
   DO K=2,NEQUA
    EQ(1)=EQ(1)+XN(K)
    K1=K1+NEQUA
!     K1 IS ACTUALLY 1K
    DEQ(K1)=1.
    EQ(K)=XN(K)-XAB(K)*XN(1)
    KK=KK+NEQUA1
    DEQ(KK)=1.
    DEQ(K)=-XAB(K)
   enddo
   IF(IDEQUA(NEQUA) >= 100)then
    EQ(NEQUA)=-XN(NEQUA)
    DEQ(NEQNEQ)=-1.
   endif
   DO JMOL=1,NUMMOL
    NCOMP=LOCJ(JMOL+1)-LOCJ(JMOL)
    IF(NCOMP == 1)cycle
    TERM=EQUILJ(JMOL)
    LOCJ1=LOCJ(JMOL)
    LOCJ2=LOCJ(JMOL+1)-1
     DO  LOCK=LOCJ1,LOCJ2
      K=KCOMPS(LOCK)
      IF(K /= NEQUA1)then
       TERM=TERM*XN(K)
      else
       TERM=TERM/XN(NEQUA)
      endif
     enddo
  EQ(1)=EQ(1)+TERM
   DO LOCK=LOCJ1,LOCJ2
    K=KCOMPS(LOCK)
    IF(K >= NEQUA1)then
     K=NEQUA
     D=-TERM/XN(K)
    else
     D=TERM/XN(K)
    endif
    EQ(K)=EQ(K)+TERM
  NEQUAK=NEQUA*K-NEQUA
  K1=NEQUAK+1
  DEQ(K1)=DEQ(K1)+D
   DO LOCM=LOCJ1,LOCJ2
    M=KCOMPS(LOCM)
    IF(M == NEQUA1)M=NEQUA
    MK=M+NEQUAK
    DEQ(MK)=DEQ(MK)+D
   enddo
  enddo
!     CORRECTION TO CHARGE EQUATION FOR NEGATIVE IONS
   K=KCOMPS(LOCJ2)
   IF(IDEQUA(K) /= 100)cycle
    DO LOCK=LOCJ1,LOCJ2
     K=KCOMPS(LOCK)
     D=TERM/XN(K)
     IF(K == NEQUA)EQ(K)=EQ(K)-TERM-TERM
     NEQUAK=NEQUA*K-NEQUA
      DO  LOCM=LOCJ1,LOCJ2
       M=KCOMPS(LOCM)
       IF(M /= NEQUA)cycle
       MK=M+NEQUAK
       DEQ(MK)=DEQ(MK)-D-D
      enddo
    enddo
   enddo
    CALL SOLVIT(DEQ,NEQUA,EQ,IDTERM)
  IFERR=0
  SCALE=100.d0
   DO K=1,NEQUA
      RATIO=ABS(EQ(K)/XN(K))
      IF(RATIO > 0.001d0)IFERR=1
      IF(EQOLD(K)*EQ(K) < 0.)EQ(K)=EQ(K)*0.69d0
      XNEQ=XN(K)-EQ(K)
      XN100=XN(K)/100.d0
      IF(XNEQ >= XN100)then
       XN100=XN(K)*100.d0
       XN(K)=XNEQ
      else
       XN(K)=XN(K)/SCALE
       IF(EQOLD(K)*EQ(K) < 0.)SCALE=SQRT(SCALE)
      endif
     EQOLD(K)=EQ(K)
    enddo
      IF(IFERR /= 1)exit
 enddo ! Finish iteration cycle
 do k=1,nequa
  XNZ(J,K)=XN(K)
 enddo
  XNATOM(J)=XN(1)
  RHO(J)=XNATOM(J)*WTMOLE*1.660d-24
  IF(IDEQUA(NEQUA) == 100)XNE(J)=XN(NEQUA)
  DO JMOL=1,NUMMOL
    NCOMP=LOCJ(JMOL+1)-LOCJ(JMOL)
    XNMOL(J,JMOL)=EQUILJ(JMOL)
    LOCJ1=LOCJ(JMOL)
    LOCJ2=LOCJ(JMOL+1)-1
   DO LOCK=LOCJ1,LOCJ2
    K=KCOMPS(LOCK)
    IF(K /= NEQUA1)then
     XNMOL(J,JMOL)=XNMOL(J,JMOL)*XN(K)
    else
     XNMOL(J,JMOL)=XNMOL(J,JMOL)/XN(NEQUA)
    endif
   enddo
  enddo
 enddo
 DO K=2,NEQUA
   ID=IDEQUA(K)
   IF(ID /= 100)then
    DO J=1,NRHOX
!     CALCULATE PARTITION FUNCTIONS
     CALL PFSAHA(J,ID,1,3)
     FRAC=ANSWER
     XNZ(J,K)=XNZ(J,K)/FRAC(J,1)/1.8786E20/SQRT((ATMASS(ID)*T(J))**3)
    enddo
   else
     DO J=1,NRHOX
      XNZ(J,K)=XNZ(J,K)/2./2.4148E15/T(J)/SQRT(T(J))
     enddo
   endif
 enddo
  DO JMOL=1,NUMMOL
   NCOMP=LOCJ(JMOL+1)-LOCJ(JMOL)
    IF(EQUIL(1,JMOL) /= 0.d0)then
     DO  J=1,NRHOX
      XNMOL(J,JMOL)=EXP(EQUIL(1,JMOL)/TKEV(J))
     enddo
    AMASS=0.d0
    LOCJ1=LOCJ(JMOL)
    LOCJ2=LOCJ(JMOL+1)-1
    DO LOCK=LOCJ1,LOCJ2
      K=KCOMPS(LOCK)
     IF(K /= NEQUA1)then
      ID=IDEQUA(K)
      IF(ID < 100)AMASS=AMASS+ATMASS(ID)
       DO J=1,NRHOX
        XNMOL(J,JMOL)=XNMOL(J,JMOL)*XNZ(J,K)
       enddo
     else
       DO  J=1,NRHOX
        XNMOL(J,JMOL)=XNMOL(J,JMOL)/XNZ(J,NEQUA)
       enddo
     endif
    enddo
      DO  J=1,NRHOX
       XNMOL(J,JMOL)=XNMOL(J,JMOL)*1.8786E20*SQRT((AMASS*T(J))**3)
      enddo
     Mol_mass(jmol-natomc)=amass
    else
     ID=CODE(JMOL)
      DO  J=1,NRHOX
       CALL PFSAHA(J,ID,NCOMP,3)
        FRAC=ANSWER
       XNMOL(J,JMOL)=XNMOL(J,JMOL)/FRAC(J,1)
      enddo
    endif
  enddo
!write(101,'(3x,152f12.2)')(code(k),k=1,nummol)
 do j=1,nrhox
!write(101,'(i3,1p152e12.1)')j,(xnmol(j,k),k=1,nummol)
  c=1.49725d-2/rho(j)
 do i=1,nmol
  x_mol(j,i)=xnmol(j,i+natomc)*c
 enddo
!write(120,*)j,x_mol(j,10),xnatom(j)
!write(120,'(i3,1p68e10.2)')j,(x_mol(j,i),i=1,68)
 enddo
END 