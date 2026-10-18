 SUBROUTINE KAPP
 use radiation_field
 use atmosphere_model, only : Nrhox
 implicit real(8) (a-h,o-z)

      do j = 1, Nrhox
       Acont(j)=0.d0
        Scont(j)=Bnu(j)
         sigmaC(j)=0.d0
       enddo
      CALL HOP
      CALL H2PLOP
      CALL HMINOP
      CALL HRAYOP
      CALL HE1OP
      CALL HE2OP
      CALL HEMIOP
      CALL HERAOP
      CALL COOLOP
      CALL LUKEOP
      CALL HOTOP
      CALL ELECOP
 END

real(8) function occupationProbability(z,charge,n,l,t,ne,nh1)
! subroutine from Shulyak's synthetic code
!  Hummer & Milhalas (1988) http://esoads.eso.org/abs/1988ApJ...331..794H, equation 4.71
!  For ions use formula from Hubeny, Hummer & Lanz (1994) http://esoads.eso.org/abs/1994A%26A...282..151H, Appendix A 
implicit none
real(8), parameter :: const = 4d0/3d0*3.1415925d0
real(8), intent(in) :: z, charge, n, l, t, ne, nh1
real(8) :: z_radiator, z_ionic, zeff, rh, k, betac, a, x, f
real(8) :: wijk_neutrals, wijk_ions
real(8) ionRadii

z_radiator = dble(charge)
zeff       = z_radiator + 1d0
rh         = ionRadii(z,zeff,n,l) + ionRadii(z,zeff,1d0,0d0)

!wijk_neutrals = exp( -const*nh1*rh**3 )
wijk_neutrals = 1.d0

if ( n <= 3d0 ) then
  k = 1d0
else
  k = 16d0/3d0*n/(n+1)**2
end if
z_ionic = 1d0
betac   = 8.3d14 * ne**(-2d0/3d0) * z_ionic**3 * k * n**(-4d0)
a       = 0.09d0*ne**(1d0/6d0)/sqrt(t)
x       = (1d0 + a)**3.15d0
f       = 0.1402d0*(x + 4d0*z_radiator*a**3)*betac**3 / (1d0 + 0.1285d0*x*betac**(3d0/2d0))

wijk_ions = f / (1d0 + f)

occupationProbability = wijk_neutrals * wijk_ions
end function occupationProbability
 
real(8) function ionRadii(z,zeff,n,l)
implicit none
real(8), intent(in) :: z, zeff, n, l
real(8), parameter :: a0 =5.291772081145377E-009 ! Bohr radius of 1-st orbit  [cm]

if ( z == 1 .and. zeff == 1d0 ) then
  ionRadii = a0/2d0/zeff * (3d0*n**2 - l*(l+1))
else
  ionRadii = 0d0
end if 

end

real(8) function hydrogenN(freq)
implicit none
real(8), intent(in) :: freq
real(8) :: lambda
!D!
! Due to the numerical uncertainities and input physical constans used the
! H-series jumps and the number of current serie determined by Balmer formulae
! may not coincide. That means that, say, getHn=1 after the limit of Balmer
! serie has already been reached. And vise versa. This causes some artificial 
! picks or depressions in the energy distribution exacetly at the H-series limits...
!
! The following is to avoid such numerical artifacts and to provide correct 
! flux representation in the region between the limit of the current serie and 
! the last visible H-line. Should be acceptable at least for 0.001 A wavelength
! step.
!D!

hydrogenN = sqrt( 109677.58375d0*2.997924580d10/freq )
lambda    = 2.997924580d18 / freq
if ( (lambda >   911.762450d0) .and. (hydrogenN < 2.d0) ) then
  hydrogenN = 1.d0
elseif ( (lambda >=  3647.049001d0) .and. (hydrogenN < 3.d0) ) then
  hydrogenN = 2.d0
elseif ( (lambda >=  8205.871001d0) .and. (hydrogenN < 4.d0) ) then
  hydrogenN = 3.d0
elseif ( (lambda > 14588.199001d0) .and. (lambda < 22794.085d0) .and. (hydrogenN < 5.d0) ) then
  hydrogenN = 4.d0
end if

end function hydrogenN 

  SUBROUTINE HOP ! modified to Hydrogen jumps by Shulyak's synthetic code
  use atmosphere_model
  use radiation_field
  use temporary_arrays
  implicit real(8) (a-h,o-z)
  DIMENSION CONT(8)
  integer(4), parameter :: maxLev = 8
  real(8), parameter :: freq_73852 = 4.05933d13 ! 73852.69 A
  real(8) nmax 

      DO J=1,NRHOX
       DO  N=1,8
        BOLT(J,N)=EXP(-(13.595d0-13.595d0/(N*N))/TKEV(J))*2.d0*(N*N)*&
        xDens(J,1,1)/RHO(J)
       enddo
       FREET(J)=XNE(J)*xDens2(J,1,2)/RHO(J)/SQRT(T(J))
       XR=xDens(J,1,1)*(2.d0/2.d0/13.595d0)*TKEV(J)/RHO(J)
       BOLTEX(J)=EXP(-13.427d0/TKEV(J))*XR
       EXLIM(J)=EXP(-13.595d0/TKEV(J))*XR
      enddo

   DO N=1,8
    CONT(N)=XKARSAS(FREQ,1.D0,N,N)
   enddo
     FREQ3=FREQ**3
     CFREE=3.6919d8/FREQ3
     C=2.815d29/FREQ3
      DO J=1,NRHOX
       if ( freq < freq_73852 ) then
        ex = exlim(j) / ehvkt(j)
       else
        ex = boltex(j)
       end if
    do n = 1, maxLev
      nlo = nint(hydrogenN(freq))
      flim = 109677.58375d0*2.997924580d10/(dble(nlo)**2.d0)
      nmax = 1d0/dble(n)**2 - freq/(109677.58375d0*2.997924580d10)
      if ( (freq < flim) .and. (nmax > 0d0) ) then
        nmax  = 1d0/sqrt(nmax)
        wnmax = occupationProbability(1d0,0d0,nmax,0d0,t(j),xNe(j),xDens2(J,1,1))
        wi    = 1d0
        wij = (wi - wnmax) / wi
        cont(n)  = wij*coulx(nlo,flim,1.d0)
      end if
    end do

       H=(CONT(7)*BOLT(J,7)+CONT(8)*BOLT(J,8)+(EX-EXLIM(J))*C&
         +COULFF(J,1)*FREET(J)*CFREE)*STIM(J)
      DO N=1,6
       H=H+CONT(N)*BOLT(J,N)*(1.d0-EHVKT(J))
      enddo
      aCont(J)=aCont(j)+H
     enddo

      RETURN
      END
!___
 real(8) FUNCTION XKARSAS(FREQ,ZEFF2,N,L)
   use KarsTab
   implicit real(8) (a-h,o-z)
  FREQLG=LOG10(FREQ/ZEFF2)
  XKARSAS=0.d0
   IF(L >= N.OR.N > 6)THEN
    IF(N > 15)then
     FREQN15(29)=LOG10(109677.576d0*2.99792458d10/N**2)
      IF(FREQLG.LT.FREQN15(29))RETURN
       DO I=2,28
        FREQN15(I)=LOG10((EKARSAS(I)+1.d0/N**2)*109677.576d0*2.99792458d10)
        IF(FREQLG > FREQN15(I))exit
       enddo
     I=29
     X=(FREQLG-FREQN15(I))/(FREQN15(I-1)-FREQN15(I))*&
      (XN(I-1,15)-XN(I,15))+XN(I,15)
     XKARSAS=EXP(X*2.30258509299405d0)/ZEFF2
     RETURN
    endif
   IF(FREQLG < FREQN(29,N))RETURN
    DO I=2,29
     IF(FREQLG > FREQN(I,N))exit
    enddo
   X=(FREQLG-FREQN(I,N))/(FREQN(I-1,N)-FREQN(I,N))*(XN(I-1,N)-XN(I,N))+XN(I,N)
   XKARSAS=EXP(X*2.30258509299405d0)/ZEFF2
   RETURN
  ENDIF
   IF(FREQLG < FREQN(29,N))RETURN
    DO I=2,29
      IF(FREQLG > FREQN(I,N))exit
    enddo
   X=(FREQLG-FREQN(I,N))/(FREQN(I-1,N)-FREQN(I,N))*&
     (XL(I-1,L+1,N)-XL(I,L+1,N))+XL(I,L+1,N)
   XKARSAS=EXP(X*2.30258509299405d0)/ZEFF2
  RETURN
 END
!___
 real(8) FUNCTION COULX(N,FREQ,Z)
  IMPLICIT REAL(8) (A-H,O-Z)
  DIMENSION A(6),B(6),C(6)
   DATA A/.9916d0,1.105d0,1.101d0,1.101d0,1.102d0,1.0986d0/
   DATA B/2.719d13,-2.375d14,-9.863d13,-5.765d13,-3.909d13,-2.704d13/
   DATA C/-2.268d30,4.077d28,1.035d28,4.593d27,2.371d27,1.229d27/
   real(8) ni
   ni=dble(n)
    IF(FREQ < Z*Z*3.28805d15/(Ni*Ni))then
     COULX=0.d0
    else
     COULX=2.815d29/FREQ/FREQ/FREQ/(Ni**5.d0)*Z**4.d0
      IF(N > 6)RETURN
      IF(N == 1)THEN
       COULX=COULX*COULBF1S(FREQ,Z)
      else
       COULX=COULX*(A(N)+(B(N)+C(N)*(Z*Z/FREQ))*(Z*Z/FREQ))
       coulx=max(coulx,0.d0)
      endif
    endif
      END
!___
 real(8) FUNCTION COULBF1S(FREQ,Z)
  IMPLICIT REAL(8) (A-H,O-Z)
  DIMENSION GAUNT1S(151)
    DATA GAUNT1S/&
     0.7973,0.8094,0.8212,0.8328,0.8439,0.8548,0.8653,0.8754,0.8852,&
     0.8946,0.9035,0.9120,0.9201,0.9278,0.9351,0.9420,0.9484,0.9544,&
     0.9601,0.9653,0.9702,0.9745,0.9785,0.9820,0.9852,0.9879,0.9903,&
     0.9922,0.9938,0.9949,0.9957,0.9960,0.9960,0.9957,0.9949,0.9938,&
     0.9923,0.9905,0.9884,0.9859,0.9832,0.9801,0.9767,0.9730,0.9688,&
     0.9645,0.9598,0.9550,0.9499,0.9445,0.9389,0.9330,0.9269,0.9206,&
     0.9140,0.9071,0.9001,0.8930,0.8856,0.8781,0.8705,0.8627,0.8546,&
     0.8464,0.8381,0.8298,0.8213,0.8128,0.8042,0.7954,0.7866,0.7777,&
     0.7685,0.7593,0.7502,0.7410,0.7318,0.7226,0.7134,0.7042,0.6951,&
     0.6859,0.6767,0.6675,0.6584,0.6492,0.6401,0.6310,0.6219,0.6129,&
     0.6039,0.5948,0.5859,0.5769,0.5680,0.5590,0.5502,0.5413,0.5324,&
     0.5236,0.5148,0.5063,0.4979,0.4896,0.4814,0.4733,0.4652,0.4572,&
     0.4493,0.4415,0.4337,0.4261,0.4185,0.4110,0.4035,0.3962,0.3889,&
     0.3818,0.3749,0.3680,0.3611,0.3544,0.3478,0.3413,0.3348,0.3285,&
     0.3222,0.3160,0.3099,0.3039,0.2980,0.2923,0.2866,0.2810,0.2755,&
     0.2701,0.2648,0.2595,0.2544,0.2493,0.2443,0.2394,0.2345,0.2298,&
     0.2251,0.2205,0.2160,0.2115,0.2072,0.2029,0.1987/
      COULBF1S=0.d0
      IF(FREQ/Z**2 < 3.28805d15)RETURN
      ELOG=LOG10(FREQ/Z**2/3.28805d15)
      I=ELOG/0.02d0
      I=MAX(MIN(I+1,150),1)
      COULBF1S=GAUNT1S(I)+(GAUNT1S(I+1)-GAUNT1S(I))/.02d0*(ELOG-(I-1)*.02d0)
      RETURN
      END
!__
 real(8)  FUNCTION COULFF(J,NZ)
  use atmosphere_model
  use radiation_field
   implicit real(8) (a-h,o-z)
   DIMENSION Z4LOG(6),A(11,12)
   DATA Z4LOG/0.d0,1.20412d0,1.90849d0,2.40824d0,2.79588d0,3.11261d0/
   DATA A/                                                 &
    5.53,5.49,5.46,5.43,5.40,5.25,5.00,4.69,4.48,4.16,3.85,&
    4.91,4.87,4.84,4.80,4.77,4.63,4.40,4.13,3.87,3.52,3.27,&
    4.29,4.25,4.22,4.18,4.15,4.02,3.80,3.57,3.27,2.98,2.70,&
    3.64,3.61,3.59,3.56,3.54,3.41,3.22,2.97,2.70,2.45,2.20,&
    3.00,2.98,2.97,2.95,2.94,2.81,2.65,2.44,2.21,2.01,1.81,&
    2.41,2.41,2.41,2.41,2.41,2.32,2.19,2.02,1.84,1.67,1.50,&
    1.87,1.89,1.91,1.93,1.95,1.90,1.80,1.68,1.52,1.41,1.30,&
    1.33,1.39,1.44,1.49,1.55,1.56,1.51,1.42,1.33,1.25,1.17,&
    0.90,0.95,1.00,1.08,1.17,1.30,1.32,1.30,1.20,1.15,1.11,&
    0.55,0.58,0.62,0.70,0.85,1.01,1.15,1.18,1.15,1.11,1.08,&
    0.33,0.36,0.39,0.46,0.59,0.76,0.97,1.09,1.13,1.10,1.08,&
    0.19,0.21,0.24,0.28,0.38,0.53,0.76,0.96,1.08,1.09,1.09/
!     ERROR CORRECTED 13APR88
!     A0.45,0.48,0.52,0.60,0.75,0.91,1.15,1.18,1.15,1.11,1.08,
      GAMLOG=10.39638d0-TLOG(J)/1.15129d0+Z4LOG(NZ)
      IGAM=MAX(MIN( INT(GAMLOG+7.d0),10),1)
      HVKTLG=(FREQLG-TLOG(J))/1.15129d0-20.63764d0
      IHVKT=MAX(MIN( INT(HVKTLG+9.d0),11),1)
      Pv=GAMLOG-(IGAM-7)
      Qv=HVKTLG-(IHVKT-9)
      COULFF = (1.d0-Pv)*((1.d0-Qv)*A(IGAM,IHVKT)+Qv*A(IGAM,IHVKT+1))+&
      Pv*((1.d0-Qv)*A(IGAM+1,IHVKT)+Qv*A(IGAM+1,IHVKT+1))
      RETURN
      END
!___
  SUBROUTINE H2PLOP
  use atmosphere_model
  use radiation_field
  implicit real(8) (a-h,o-z)
   IF(FREQ > 3.28805d15)RETURN
   FR=-3.0233d3+(3.7797d2+(-1.82496d1+(3.9207d-1-3.1672d-3*FREQLG)*FREQLG)&
      *FREQLG)*FREQLG
   ES=-7.342d-3+(-2.409d-15+(1.028d-30+(-4.230d-46+(1.224d-61-&
      1.351d-77*FREQ)*FREQ)*FREQ)*FREQ)*FREQ
   DO J=1,NRHOX
    H2P=EXP(-ES/TKEV(J)+FR)*xDens(J,1,1)*2.d0*xDens(J,1,2)/RHO(J)*STIM(J)
   aCont(j)=aCont(j)+H2p
   enddo
  END
!___
  SUBROUTINE HMINOP
  use atmosphere_model
  use radiation_field
  use temporary_arrays
   implicit real(8) (a-h,o-z)
   real(8)WAVELOG(1),FFTLOG(1),WAVE(1),HMINBF(1)
   DIMENSION WBF(85),BF(85),FFLOG(22,11),FF(11,22)
   DIMENSION FFBEG(11,11),FFEND(11,11),FFTT(11),WFFLOG(22)
   DIMENSION THETAFF(11),WAVEK(22)
   EQUIVALENCE (FF(1,1),FFBEG(1,1)),(FF(1,12),FFEND(1,1))
!     FROM MATHISEN (1984), AFTER WISHART(1979) AND BROAD AND REINHARDT (1976)
      DATA WBF/  18.00,  19.60,  21.40,  23.60,  26.40,  29.80,  34.30,&
         40.40,  49.10,  62.60, 111.30, 112.10, 112.67, 112.95, 113.05,&
        113.10, 113.20, 113.23, 113.50, 114.40, 121.00, 139.00, 164.00,&
        175.00, 200.00, 225.00, 250.00, 275.00, 300.00, 325.00, 350.00,&
        375.00, 400.00, 425.00, 450.00, 475.00, 500.00, 525.00, 550.00,&
        575.00, 600.00, 625.00, 650.00, 675.00, 700.00, 725.00, 750.00,&
        775.00, 800.00, 825.00, 850.00, 875.00, 900.00, 925.00, 950.00,&
        975.00,1000.00,1025.00,1050.00,1075.00,1100.00,1125.00,1150.00,&
       1175.00,1200.00,1225.00,1250.00,1275.00,1300.00,1325.00,1350.00,&
       1375.00,1400.00,1425.00,1450.00,1475.00,1500.00,1525.00,1550.00,&
       1575.00,1600.00,1610.00,1620.00,1630.00,1643.91/
      DATA BF/   0.067,  0.088,  0.117,  0.155,  0.206,  0.283,  0.414,&
         0.703,   1.24,   2.33,  11.60,  13.90,  24.30,  66.70,  95.00,&
         56.60,  20.00,  14.60,   8.50,   7.10,   5.43,   5.91,   7.29,&
         7.918,  9.453,  11.08,  12.75,  14.46,  16.19,  17.92,  19.65,&
         21.35,  23.02,  24.65,  26.24,  27.77,  29.23,  30.62,  31.94,&
         33.17,  34.32,  35.37,  36.32,  37.17,  37.91,  38.54,  39.07,&
         39.48,  39.77,  39.95,  40.01,  39.95,  39.77,  39.48,  39.06,&
         38.53,  37.89,  37.13,  36.25,  35.28,  34.19,  33.01,  31.72,&
         30.34,  28.87,  27.33,  25.71,  24.02,  22.26,  20.46,  18.62,&
         16.74,  14.85,  12.95,  11.07,  9.211,  7.407,  5.677,  4.052,&
         2.575,  1.302, 0.8697, 0.4974, 0.1989,    0. /
!    Bell and Berrington J.Phys.B,vol. 20, 801-806,1987.
      DATA WAVEK/.50,.40,.35,.30,.25,.20,.18,.16,.14,.12,.10,.09,.08,&
       .07,.06,.05,.04,.03,.02,.01,.008,.006/
      DATA THETAFF/0.5,  0.6, 0.8,  1.0,  1.2,  1.4,  1.6,  1.8,  2.0,  2.8,  3.6/
      DATA FFBEG/&
      .0178,.0222,.0308,.0402,.0498,.0596,.0695,.0795,.0896, .131, .172,&   !1823
      .0228,.0280,.0388,.0499,.0614,.0732,.0851,.0972, .110, .160, .211,&   !2278
      .0277,.0342,.0476,.0615,.0760,.0908, .105, .121, .136, .199, .262,&   !2604
      .0364,.0447,.0616,.0789,.0966, .114, .132, .150, .169, .243, .318,&   !3038
      .0520,.0633,.0859, .108, .131, .154, .178, .201, .225, .321, .418,&   !3645
      .0791,.0959, .129, .161, .194, .227, .260, .293, .327, .463, .602,&   !4557
      .0965, .117, .157, .195, .234, .272, .311, .351, .390, .549, .711,&   !5063
       .121, .146, .195, .241, .288, .334, .381, .428, .475, .667, .861,&   !5696
       .154, .188, .249, .309, .367, .424, .482, .539, .597, .830, 1.07,&   !6510
       .208, .250, .332, .409, .484, .557, .630, .702, .774, 1.06, 1.36,&   !7595
       .293, .354, .468, .576, .677, .777, .874, .969, 1.06, 1.45, 1.83/   !9113
      DATA FFEND/&
       .358, .432, .572, .702, .825, .943, 1.06, 1.17, 1.28, 1.73, 2.17,&  !10126
       .448, .539, .711, .871, 1.02, 1.16, 1.29, 1.43, 1.57, 2.09, 2.60,&  !11392
       .579, .699, .924, 1.13, 1.33, 1.51, 1.69, 1.86, 2.02, 2.67, 3.31,&  !13019
       .781, .940, 1.24, 1.52, 1.78, 2.02, 2.26, 2.48, 2.69, 3.52, 4.31,&  !15189
       1.11, 1.34, 1.77, 2.17, 2.53, 2.87, 3.20, 3.51, 3.80, 4.92, 5.97,&  !18227
       1.73, 2.08, 2.74, 3.37, 3.90, 4.50, 5.01, 5.50, 5.95, 7.59, 9.06,&  !22784
       3.04, 3.65, 4.80, 5.86, 6.86, 7.79, 8.67, 9.50, 10.3, 13.2, 15.6,&  !30378
       6.79, 8.16, 10.7, 13.1, 15.3, 17.4, 19.4, 21.2, 23.0, 29.5, 35.0,&  !45567
       27.0, 32.4, 42.6, 51.9, 60.7, 68.9, 76.8, 84.2, 91.4, 117., 140.,&  !91134
       42.3, 50.6, 66.4, 80.8, 94.5, 107., 120., 131., 142., 183., 219.,& !113918
       75.1, 90.0, 118., 144., 168., 191., 212., 234., 253., 325., 388./ !151890
      DATA ISTART/0/
      IF(ISTART == 0)THEN
       ISTART=1
       DO IWAVE=1,22
!     91.134 NUMBER TAKEN FROM BELL AND BERRINGTON
       WFFLOG(IWAVE)=LOG(91.134D0/WAVEK(IWAVE))
        DO  ITHETA=1,11
         FFLOG(IWAVE,ITHETA)=LOG(FF(ITHETA,IWAVE)*1.d-26)
        enddo
       enddo
      ENDIF

    DO  J=1,NRHOX
     THETA(J)=5040./T(J)
!     .754209 HOTOP AND LINEBERGER J.PHYS.CHEM.REF.DATA VOL 14,731-752,1985.
     XHMIN(J)=EXP(.754209d0/TKEV(J))/(2.d0*2.4148d15*T(J)*SQRT(T(J)))*&
      xDens(J,1,1)*XNE(J)
    enddo

    WAVE=2.99792458d17/FREQ
    WAVELOG=LOG(WAVE)
     DO  ITHETA=1,11
      CALL LINTER(WFFLOG,FFLOG(1,ITHETA),22,WAVELOG,FFTLOG,1)
      FFTT(ITHETA)=EXP(FFTLOG(1))/THETAFF(ITHETA)*5040.d0*1.380658d-16
     enddo
      HMINBF=0.d0
      IF(FREQ > 1.82365d14)MAXWAVE=MAP1(WBF,BF,85,WAVE,HMINBF,1)
       DO J=1,NRHOX
        CALL LINTER(THETAFF,FFTT,11,THETA(J),FFTHETA(J),1)
        HMINFF=FFTHETA(J)*xDens(J,1,1)*2.d0*XNE(J)/RHO(J)
        H=HMINBF(1)*1.d-18*(1.d0-EHVKT(J))*XHMIN(J)/RHO(J)
      aCont(j)=aCont(j)+H+HMINFF
       enddo
      END
   SUBROUTINE HRAYOP
  use atmosphere_model
  use radiation_field
  implicit real(8) (a-h,o-z)
    WAVE=2.997925d18/MIN(FREQ,2.463D15)
    WW=WAVE**2
    SIG=(5.799d-13+1.422d-6/WW+2.784d0/(WW*WW))/(WW*WW)
     DO J=1,NRHOX
      sigmaC(J)=sigmaC(j)+SIG*xDens(J,1,1)*2.d0/RHO(J)
     enddo
    END
!------
  SUBROUTINE HE1OP
  use atmosphere_model
  use radiation_field
  use temporary_arrays
  implicit real(8) (a-h,o-z)
!     REQUIRES FUNCTIONS COULX AND COULFF
   DIMENSION CHI(10),HEFREQ(10),TRANS(10),G(10)
   DIMENSION TRANS1S(10),TRANSN(27)

    DATA G/1.d0,3.d0,1.d0,9.d0,3.d0,3.d0,1.d0,9.d0,20.d0,3.d0/
    DATA HEFREQ/5.945209d15,1.152844d15,.9603331d15,.8761076d15,&
       .8147104d15,.4519048d15,.4030971d15,.3821191d15,.3660215d15,&
       .3627891d15/
    DATA CHI/0.d0,19.819d0,20.615d0,20.964d0,21.217d0,22.718d0,22.920d0,&
             23.006d0,23.073d0,23.086d0/

   RYD=109722.273d0*2.99792458d10
    DO J=1,NRHOX
     DO  N=1,10
      BOLT(J,N)=EXP(-CHI(N)/TKEV(J))*G(N)*xDens(J,2,1)/1.49725d-2
     enddo
     DO  N=4,27
      BOLTN(J,N)=EXP(-24.587d0*(1.d0-1.d0/N**2)/TKEV(J))*4.d0*N**2*&
       xDens(J,2,1)/1.49725d-2
     enddo
      FREET(J)=XNE(J)*xDens2(J,2,2)/1.49725d-2/SQRT(T(J))
      XR=xDens(J,2,1)*(4.d0/2.d0/13.595d0)*TKEV(J)/1.49725d-2
      BOLTEX(J)=EXP(-23.730d0/TKEV(J))*XR
      EXLIM(J)=EXP(-24.587d0/TKEV(J))*XR
    enddo

   FREQ3=FREQ**3
   CFREE=3.6919d8/FREQ3
   C=2.815d29/FREQ3
    DO 15 IMIN=1,10
      IF(HEFREQ(IMIN).LE.FREQ)GO TO 16
   15 CONTINUE
      IMIN=0
      GO TO 40
 16 if(imin == 1)TRANS(1)=CROSSHE(FREQ)
    if(imin == 2)TRANS(2)=HE12s3S(FREQ)
    if(imin == 3)TRANS(3)=HE12s1S(FREQ)
    if(imin == 4)TRANS(4)=HE12p3P(FREQ)
    if(imin == 5)TRANS(5)=HE12p1P(FREQ)
!     1s3s 3S
    if(imin == 6)TRANS(6)=XKARSAS(FREQ,1.236439D0,3,0)
!     1s3s 1S
    if(imin == 7)TRANS(7)=XKARSAS(FREQ,1.102898D0,3,0)
!     1s3p 3P
    if(imin == 8)TRANS(8)=XKARSAS(FREQ,1.045499D0,3,1)
!     1s3d 3D+1D
    if(imin == 9)TRANS(9)=XKARSAS(FREQ,1.001427D0,3,2)
!     1s3p 1P
    if(imin == 10)TRANS(10)=XKARSAS(FREQ,.9926D0,3,1)
!     HeII n=2
   31 ELIM=527490.06d0
      FREQHE=(ELIM-171135.000d0)*2.99792458d10
      IF(FREQ < FREQHE)GO TO 32
      ZEFF2=FREQHE/RYD
      TRANS(5)=TRANS(5)+XKARSAS(FREQ,ZEFF2,1,0)
      FREQHE=(ELIM-169087.d0)*2.99792458d10
      IF(FREQ < FREQHE)GO TO 32
      ZEFF2=FREQHE/RYD
      TRANS(4)=TRANS(4)+XKARSAS(FREQ,ZEFF2,1,0)
      FREQHE=(ELIM-166277.546d0)*2.99792458d10
      IF(FREQ < FREQHE)GO TO 32
      ZEFF2=FREQHE/RYD
      TRANS(3)=TRANS(3)+XKARSAS(FREQ,ZEFF2,1,0)
      FREQHE=(ELIM-159856.069d0)*2.99792458d10
      IF(FREQ < FREQHE)GO TO 32
      ZEFF2=FREQHE/RYD
      TRANS(2)=TRANS(2)+XKARSAS(FREQ,ZEFF2,1,0)
!     HeII n=3
   32 ELIM=588451.59d0
      FREQHE=(ELIM-186209.471d0)*2.99792458d10
      IF(FREQ < FREQHE)GO TO 40
      ZEFF2=FREQHE/RYD
      TRANS(10)=TRANS(10)+XKARSAS(FREQ,ZEFF2,1,0)
      FREQHE=(ELIM-186101.d0)*2.99792458d10
      IF(FREQ < FREQHE)GO TO 40
      ZEFF2=FREQHE/RYD
      TRANS(9)=TRANS(9)+XKARSAS(FREQ,ZEFF2,1,0)
      FREQHE=(ELIM-185564.d0)*2.99792458d10
      IF(FREQ < FREQHE)GO TO 40
      ZEFF2=FREQHE/RYD
      TRANS(8)=TRANS(8)+XKARSAS(FREQ,ZEFF2,1,0)
      FREQHE=(ELIM-184864.d0)*2.99792458d10
      IF(FREQ < FREQHE)GO TO 40
      ZEFF2=FREQHE/RYD
      TRANS(7)=TRANS(7)+XKARSAS(FREQ,ZEFF2,1,0)
      FREQHE=(ELIM-183236.d0)*2.99792458d10
      IF(FREQ < FREQHE)GO TO 40
      ZEFF2=FREQHE/RYD
      TRANS(6)=TRANS(6)+XKARSAS(FREQ,ZEFF2,1,0)
      IF(FREQ < 1.25408d+16)GO TO 40
      DO  N=4,27
       ZEFF2=4.d0-3.d0/N**2
       TRANSN(N)=XKARSAS(FREQ,ZEFF2,1,0)
      enddo
   40 DO 45 J=1,NRHOX
      EX=BOLTEX(J)
      IF(FREQ < 2.055d14)EX=EXLIM(J)/EHVKT(J)
      HE1=(EX-EXLIM(J))*C
      IF(IMIN /= 0)then
        DO  N=IMIN,10
         HE1=HE1+TRANS(N)*BOLT(J,N)
        enddo
         IF(FREQ >= 1.25408d+16)then
          DO  N=4,27
           HE1=HE1+TRANSN(N)*BOLTN(J,N)
          enddo
         endif
      endif
   45 aCont(j)=aCont(j)+(HE1+COULFF(J,1)*FREET(J)*CFREE)*STIM(J)
      RETURN
      END
!---
real(8) FUNCTION CROSSHE(FREQ)
!     Marr, G.V. and West, J.B. Atomic Data and Nuclear Data Tables,
!     vol 18, 497-508, 1976.
  use temporary_arrays
  IMPLICIT REAL(8) (A-H,O-Z)
   DIMENSION X505(92),X50(16),X20(11),X10(21)
   DATA X505/&
    7.58, 7.46, 7.33, 7.19, 7.06, 6.94, 6.81, 6.68, 6.55, 6.43,&
    6.30, 6.18, 6.05, 5.93, 5.81, 5.69, 5.57, 5.45, 5.33, 5.21,&
    5.10, 4.98, 4.87, 4.76, 4.64, 4.53, 4.42, 4.31, 4.20, 4.09,&
    4.00, 3.88, 3.78, 3.68, 3.57, 3.47, 3.37, 3.27, 3.18, 3.08,&
    2.98, 2.89, 2.80, 2.70, 2.61, 2.52, 2.44, 2.35, 2.26, 2.18,&
    2.10, 2.02, 1.94, 1.86, 1.78, 1.70, 1.63, 1.55, 1.48, 1.41,&
    1.34, 1.28, 1.21, 1.14, 1.08, 1.02, .961, .903, .847, .792,&
    .738, .687, .637, .588, .542, .497, .454, .412, .373, .335,&
    .299, .265, .233, .202, .174, .147, .123, .100,.0795,.0609,&
    .0443,.0315/
   DATA X50/.0315,.0282,.0250,.0220,.0193,.0168,.0145,.0124,.0105,&
   .00885,.00736,.00604,.00489,.00389,.00303,.00231/
   DATA X20/.00231,.00199,.00171,.00145,.00122,.00101,.000832,&
   .000673,.000535,.000417,.000318/
   DATA X10/.000318,.000274,.000235,.000200,.000168,.000139,.000115,&
   .000093,.000074,.000057,.000044,.000032,.000023,.000016,.000010,&
   .000006,.000003,.000001,.0000006,.0000003,0./
     CROSSHE=0.d0
      IF(FREQ < 5.945209d15)RETURN
      WAVE=2.99792458d18/FREQ
      IF(WAVE > 50.)THEN
       I=93.d0-(WAVE-50.d0)/5.d0
       I=MIN(92,MAX(2,I))
       CROSSHE=((WAVE-(92-I)*5-50)/5.d0*(X505(I-1)-X505(I))+X505(I))*1.d-18
      RETURN
      ENDIF
      IF(WAVE.GT.20.)THEN
      I=17.-(WAVE-20.)/2.
      I=MIN(16,MAX(2,I))
      CROSSHE=((WAVE-(16-I)*2-20)/2.d0*(X50(I-1)-X50(I))+X50(I))*1.d-18
      RETURN
      ENDIF
      IF(WAVE.GT.10.)THEN
      I=12.-(WAVE-10.)/1.
      I=MIN(11,MAX(2,I))
      CROSSHE=((WAVE-(11-I)*1-10)/1.d0*(X20(I-1)-X20(I))+X20(I))*1.d-18
      RETURN
      ENDIF
      I=22.-WAVE/.5
      I=MIN(21,MAX(2,I))
      CROSSHE=((WAVE-(21-I)*.5)/.5d0*(X10(I-1)-X10(I))+X10(I))*1.d-18
      RETURN
      END
!-----
 real(8)  FUNCTION HE111S(FREQ)
!     FOLLOWING MATHISEN
  IMPLICIT REAL(8) (A-H,O-Z)
  DIMENSION W(64),X(64)
    DATA W/&
     504.3, 501.5, 498.7, 493.3, 488.1, 480.3, 477.8, 454.0, 443.0,&
     395.0, 356.4, 348.2, 324.6, 302.0, 298.1, 275.6, 260.6, 256.2,&
     239.4, 224.6,  220.,  215,   210.,  205.,  200.,  195.,  190.,&
     185.,  180.,  175.,  170.,  165.,  160.,  155.,  150.,  145.,&
     135.,  130.,  125.,  120.,  115.,  110.,  105.,  100.,   95.,&
     90.,   85.,   80.,   75.,   70.,   65.,   60.,   55.,   50.,&
      45.,   40.,   35.,   30.,   25.,   20.,   15.,   10.,  5.,0./
    DATA X/&
     7.346, 7.317, 7.259, 7.143, 7.030, 6.857, 6.800, 6.284, 6.041,&
     4.977, 4.138, 3.961, 3.474, 3.025, 2.945, 2.522, 2.259, 2.179,&
     1.901, 1.684, 1.61 , 1.53 , 1.45 , 1.38 , 1.30 , 1.22 , 1.14 ,&
     1.08 , 1.02 , 0.961, 0.903, 0.847, 0.792, 0.738, 0.687, 0.637,&
     0.542, 0.497, 0.454, 0.412, 0.373, 0.335, 0.299, 0.265, 0.233,&
     0.202, 0.174, 0.147, 0.124, 0.103,0.0840,0.0676,0.0535,0.0414,&
    .0311,.0266,.0158,.0104,.00637,.00349,.00161,.00054,.000083,0./
      HE111S=0.d0
      IF(FREQ < 5.945209d15)RETURN
      WAVE=2.99792458d18/FREQ
      DO I=2,64
       IF(WAVE > W(I))exit
      enddo
      HE111S=((WAVE-W(I))/(W(I-1)-W(I))*(X(I-1)-X(I))+X(I))*1.d-18
      RETURN
      END
!-----
 real(8)  FUNCTION HE12S1S(FREQ)
   IMPLICIT REAL(8) (A-H,O-Z)
   DIMENSION FREQ1S(16),X1S(16)
    DATA FREQ1S/&
     15.947182,   15.913654,   15.877320,   15.837666,   15.794025,&
     15.745503,   15.690869,   15.628361,   15.555317,   15.467455,&
     15.357189,   15.289399,   15.251073,   15.209035,   15.162487,&
     14.982421/
    DATA X1S/&
     -19.635557,  -19.159345,  -18.958474,  -18.809535,  -18.676481,&
     -18.546006,  -18.410962,  -18.264821,  -18.100205,  -17.909165,&
     -17.684370,  -17.557867,  -17.490360,  -17.417876,  -17.349386,&
     -17.084441/
      HE12S1S=0.d0
      IF(FREQ < 32033.214d0*2.99792458d10)RETURN
      IF(FREQ > 2.4d0*109722.267d0*2.99792458d10)then
       WAVENO=FREQ/2.99792458d10
       EK=(WAVENO-32033.214d0)/109722.267d0
       EPS=2.d0*(EK-2.612316d0)/.00322d0
       HE12S1S=.008175d0*(484940.d0/WAVENO)**2.71d0*8.067d-18*&
               (EPS+76.21d0)**2/(1.d0+EPS**2)
       RETURN
      endif
      FREQLG=LOG10(FREQ)
      DO I=2,16
       IF(FREQLG > FREQ1S(I))exit
      enddo
      X=(FREQLG-FREQ1S(I))/(FREQ1S(I-1)-FREQ1S(I))*(X1S(I-1)-X1S(I))+X1S(I)
      HE12S1S=EXP(X*2.30258509299405d0)
      END
!-----
 real(8)  FUNCTION HE12S3S(FREQ)
  IMPLICIT REAL(8) (A-H,O-Z)
   DIMENSION FREQ3S(16),X3S(16)
    DATA FREQ3S/&
     15.956523,   15.923736,   15.888271,   15.849649,   15.807255,&
     15.760271,   15.707580,   15.647601,   15.577992,   15.495055,&
     15.392451,   15.330345,   15.295609,   15.257851,   15.216496,&
     15.061770/
    DATA X3S/&
     -18.426022,  -18.610700,  -18.593051,  -18.543304,  -18.465513,&
     -18.378707,  -18.278574,  -18.164329,  -18.033346,  -17.882435,&
     -17.705542,  -17.605584,  -17.553459,  -17.500667,  -17.451318,&
     -17.266686/
   HE12S3S=0.d0
    IF(FREQ < 38454.691d0*2.99792458d10)RETURN
    IF(FREQ > 2.4d0*109722.267d0*2.99792458d10)then
     WAVENO=FREQ/2.99792458d10
      EK=(WAVENO-38454.691d0)/109722.267d0
      EPS=2.d0*(EK-2.47898d0)/.000780d0
      HE12S3S=.01521d0*(470310.d0/WAVENO)**3.12d0*8.067d-18*&
              (EPS-122.4d0)**2/(1.d0+EPS**2)
      RETURN
    endif
      FREQLG=LOG10(FREQ)
      DO  I=2,16
       IF(FREQLG > FREQ3S(I))exit
      enddo
     X=(FREQLG-FREQ3S(I))/(FREQ3S(I-1)-FREQ3S(I))*(X3S(I-1)-X3S(I))+X3S(I)
      HE12S3S=EXP(X*2.30258509299405d0)
      END
!-----
 real(8) FUNCTION HE12P1P(FREQ)
  IMPLICIT REAL(8) (A-H,O-Z)
  DIMENSION FREQ1P(16),X1P(16)
   DATA FREQ1P/&
    15.939981,   15.905870,   15.868850,   15.828377,   15.783742,&
    15.733988,   15.677787,   15.613218,   15.537343,   15.445346,&
    15.328474,   15.255641,   15.214064,   15.168081,   15.116647,&
    14.911002/
   DATA X1P/&
    -18.798876,  -19.685922,  -20.011664,  -20.143030,  -20.091354,&
    -19.908333,  -19.656788,  -19.367745,  -19.043016,  -18.674484,&
    -18.240861,  -17.989700,  -17.852015,  -17.702677,  -17.525347,&
    -16.816344/
    HE12P1P=0.d0
     IF(FREQ < 27175.76d0*2.99792458d10)RETURN
      IF(FREQ > 2.4d0*109722.267d0*2.99792458d10)then
       WAVENO=FREQ/2.99792458d10
       EK=(WAVENO-27175.76d0)/109722.267d0
       EPS1S=2.d0*(EK-2.446534d0)/.01037d0
       EPS1D=2.d0*(EK-2.59427d0)/.00538d0
       HE12p1P=.0009487d0*(466750.d0/WAVENO)**3.69d0*8.067d-18*&
       ((EPS1S-29.30d0)**2/(1.d0+EPS1S**2)+(EPS1D+172.4d0)**2/(1.d0+EPS1D**2))
      RETURN
     endif
      FREQLG=LOG10(FREQ)
      DO  I=2,16
       IF(FREQLG > FREQ1P(I))exit
      enddo
     X=(FREQLG-FREQ1P(I))/(FREQ1P(I-1)-FREQ1P(I))*(X1P(I-1)-X1P(I))+X1P(I)
     HE12P1P=EXP(X*2.30258509299405d0)
  END
!------
 real(8) FUNCTION HE12P3P(FREQ)
  IMPLICIT REAL(8) (A-H,O-Z)
  DIMENSION FREQ3P(16),X3P(16)
   DATA FREQ3P/&
    15.943031,   15.909169,   15.872441,   15.832318,   15.788107,&
    15.738880,   15.683351,   15.619667,   15.545012,   15.454805,&
    15.340813,   15.270195,   15.230054,   15.185821,   15.136567,&
    14.942557/
   DATA X3P/&
    -19.791021,  -19.697886,  -19.591421,  -19.471855,  -19.337053,&
    -19.183958,  -19.009750,  -18.807990,  -18.570571,  -18.288361,&
    -17.943476,  -17.738737,  -17.624154,  -17.497163,  -17.403183,&
    -17.032999/
   HE12P3P=0.d0
    IF(FREQ < 29223.753d0*2.99792458d10)RETURN
      FREQLG=LOG10(FREQ)
      DO I=2,16
       IF(FREQLG > FREQ3P(I))exit
      enddo
     X=(FREQLG-FREQ3P(I))/(FREQ3P(I-1)-FREQ3P(I))*(X3P(I-1)-X3P(I))+X3P(I)
     HE12P3P=EXP(X*2.30258509299405d0)
  END
!----
  SUBROUTINE HE2OP
  use atmosphere_model
  use radiation_field
  use temporary_arrays
  implicit real(8) (a-h,o-z)
!     REQUIRES FUNCTIONS COULX AND COULFF
!     FREQUENCIES ARE 4X HYDROGEN,CHI ARE FOR ION POT=54.403
   DIMENSION CONT(9)

      DO J=1,NRHOX
        DO N=1,9
         BOLT(J,N)=EXP(-(54.403d0-54.403d0/(N*N))/TKEV(J))*2.d0*(N*N)*&
                   xDens(J,2,2)/1.49725d-2
        enddo
       FREET(J)=XNE(J)*xDens2(J,2,3)/1.49725d-2/SQRT(T(J))
       XR=xDens(J,2,2)/1.49725d-2*(2.d0/2.d0/13.595d0)*TKEV(J)
       BOLTEX(J)=EXP(-53.859d0/TKEV(J))*XR
       EXLIM(J)=EXP(-54.403d0/TKEV(J))*XR
      enddo

   DO N=1,9
    CONT(N)=XKARSAS(FREQ,4.D0,N,N)
   enddo
    FREQ3=FREQ**3
    CFREE=3.6919d8/FREQ3*4.d0
    C=2.815d29*2.d0*2.d0/FREQ3
     DO J=1,NRHOX
      EX=BOLTEX(J)
      IF(FREQ.LT.1.31522E14)EX=EXLIM(J)/EHVKT(J)
      HE2=(EX-EXLIM(J))*C
       DO N=1,9
        HE2=HE2+CONT(N)*BOLT(J,N)
       enddo
       aCont(J)=aCont(j)+(HE2+COULFF(J,2)*CFREE*FREET(J))*STIM(J)
     enddo
  END
!------
 SUBROUTINE HEMIOP
  use atmosphere_model
  use radiation_field
  implicit real(8) (a-h,o-z)
      A=3.397d-46+(-5.216d-31+7.039d-15/FREQ)/FREQ
      B=-4.116d-42+(1.067d-26+8.135d-11/FREQ)/FREQ
      C=5.081d-37+(-8.724d-23-5.659d-8/FREQ)/FREQ
       DO J=1,NRHOX
        aCont(j)=aCont(j)+(A*T(J)+B+C/T(J))*XNE(J)*xDens(J,2,1)/1.49725d-2
       enddo
  END
!------
  SUBROUTINE HERAOP
  use atmosphere_model
  use radiation_field
  implicit real(8) (a-h,o-z)
    WAVE=2.997925d18/MIN(FREQ,5.15D15)
    WW=WAVE**2
    SIG=5.484d-14/WW/WW*(1.d0+(2.44d5+5.94d10/(WW-2.90d5))/WW)**2
      DO  J=1,NRHOX
       sigmaC(J)=sigmac(j)+SIG*xDens(J,2,1)/1.49725d-2
      enddo
   END
!------
  SUBROUTINE COOLOP
  use atmosphere_model
  use radiation_field
  implicit real(8) (a-h,o-z)
!     SI1,MG1,AL1,C1,FE1
  REAL(8) MG1OP
   IF(FREQ > 3.28805E15)RETURN
    DO J=1,NRHOX
        p_C=xDens(j,6,1)/1.49725d-2
        p_Al=xDens(j,13,1)/1.49725d-2
        p_Mg=xDens(j,12,1)/1.49725d-2
        p_Si=xDens(j,14,1)/1.49725d-2
        p_Fe=xDens(j,26,1)/1.49725d-2
     aCont(J)=aCont(j)+(C1OP(J)*p_C+MG1OP(J)*p_Mg+AL1OP(J)*p_Al&
              +SI1OP(J)*p_Si+FE1OP(J)*p_Fe)*STIM(J)
     enddo
   END
!-----
 real(8) FUNCTION C1OP(J)
  use atmosphere_model
  use radiation_field
  use temporary_arrays
  implicit real(8) (a-h,o-z)
!     CROSS-SECTION TIMES THE PARTITION FUNCTION

   save x1444,x1240,x1100
    DATA FREQ1/0./

       DO K=1,NRHOX
        C1240(K)=5.d0*EXP(-1.264d0/TKEV(K))
        C1444(K)=EXP(-2.683d0/TKEV(K))
      enddo

    IF(FREQ /= FREQ1)then
      X1444=0.
      X1240=0.
      X1100=0.
      IF(FREQ >= 2.7254d15)X1100=SEATON(2.7254D15,1.219D-17,2.D0,3.317D0)
      IF(FREQ >= 2.4196d15)X1240=SEATON(2.4196D15,1.03D-17,1.5D0,2.789D0)
      IF(FREQ >= 2.0761d15)X1444=SEATON(2.0761D15,9.59D-18,1.5D0,3.501D0)
      FREQ1=FREQ
     endif
   C1OP=X1100*9.d0+X1240*C1240(J)+X1444*C1444(J)
 END
!_____
 real(8) FUNCTION SEATON(FREQ0,XSECT,POWER,A)
  use atmosphere_model
  use radiation_field
  implicit real(8) (a-h,o-z)
  SEATON=XSECT*(A+(1.d0-A)*(FREQ0/FREQ))*SQRT((FREQ0/FREQ)**&
      ( INT(2.d0*POWER+.01d0)))
  END
!-------
 real(8) FUNCTION MG1OP(J)
  use atmosphere_model
  use radiation_field
  use temporary_arrays
  implicit real(8) (a-h,o-z)
!     CROSS-SECTION TIMES THE PARTITION FUNCTION
  DIMENSION FLOG(9),FREQMG(7),PEACH(7,15),X(7),TLG(7)
  DATA PEACH/&
!     4000     5000     6000     7000     8000     9000    10000     WAVE(A)
    -42.474, -42.350, -42.109, -41.795, -41.467, -41.159, -40.883,&   !1500
    -41.808, -41.735, -41.582, -41.363, -41.115, -40.866, -40.631,&   !1550
    -41.273, -41.223, -41.114, -40.951, -40.755, -40.549, -40.347,&   !1621
    -45.583, -44.008, -42.957, -42.205, -41.639, -41.198, -40.841,&   !1622
    -44.324, -42.747, -41.694, -40.939, -40.370, -39.925, -39.566,&   !2513
    -50.969, -48.388, -46.630, -45.344, -44.355, -43.568, -42.924,&   !2514
    -50.633, -48.026, -46.220, -44.859, -43.803, -42.957, -42.264,&   !3756
    -53.028, -49.643, -47.367, -45.729, -44.491, -43.520, -42.736,&   !3757
    -51.785, -48.352, -46.050, -44.393, -43.140, -42.157, -41.363,&   !6549
    -52.285, -48.797, -46.453, -44.765, -43.486, -42.480, -41.668,&   !6550
    -52.028, -48.540, -46.196, -44.507, -43.227, -42.222, -41.408,&   !7234
    -52.384, -48.876, -46.513, -44.806, -43.509, -42.488, -41.660,&   !7235
    -52.363, -48.856, -46.493, -44.786, -43.489, -42.467, -41.639,&   !7291
    -54.704, -50.772, -48.107, -46.176, -44.707, -43.549, -42.611,&   !7292
    -54.359, -50.349, -47.643, -45.685, -44.198, -43.027, -42.418/   !9000
   DATA FREQMG/1.9341452d15,1.8488510d15,1.1925797d15,7.9804046d14,&
       4.5772110d14,4.1440977d14,4.1113514d14/
   DATA FLOG/35.23123,35.19844,35.15334,34.71490,34.31318,33.75728,&
            33.65788,33.64994,33.43947/
   DATA TLG/8.29405,8.51719,8.69951,8.85367,8.98720,9.10498,9.21034/
   DATA FREQ1/0./
   save x

      DO K=1,NRHOX
       N=MAX(MIN(6, INT(T(K)/1000.)-3),1)
       NT(K)=N
       DT(K)=(TLOG(K)-TLG(N))/(TLG(N+1)-TLG(N))
      enddo

   IF(FREQ /= FREQ1)then
    FREQ1=FREQ
      DO  N=1,7
       IF(FREQ > FREQMG(N))exit
      enddo
      N=n
    D=(FREQLG-FLOG(N))/(FLOG(N+1)-FLOG(N))
      IF(N > 2)N=2*N-2
      D1=1.d0-D
       DO IT=1,7
        X(IT)=PEACH(IT,N+1)*D+PEACH(IT,N)*D1
       enddo
    endif
     N=NT(J)
      MG1OP=EXP(X(N)*(1.d0-DT(J))+X(N+1)*DT(J))
  END
!------
 real(8) FUNCTION AL1OP(J)
  use atmosphere_model
  use radiation_field
  IMPLICIT REAL(8) (A-H,O-Z)
!     CROSS-SECTION TIMES THE PARTITION FUNCTION
      AL1OP=0.d0
      IF(FREQ >= 1.443d15)AL1OP=6.5d-17*(1.443d15/FREQ)**5*6.d0
  END
!-----
 real(8)  FUNCTION SI1OP(J)
  use atmosphere_model
  use radiation_field
  use temporary_arrays
  implicit real(8) (a-h,o-z)
!     CROSS-SECTION TIMES THE PARTITION FUNCTION
  DIMENSION FLOG(11),FREQSI(9),X(9),TLG(9)
  DIMENSION PEACH(9,19)
  DATA PEACH/&
!  4000   5000   6000   7000   8000   9000   10000  11000  12000   WAVE(A)
  38.136,38.138,38.140,38.141,38.143,38.144,38.144,38.145,38.145,&  !1200
  37.834,37.839,37.843,37.847,37.850,37.853,37.855,37.857,37.858,&  !1400
  37.898,37.898,37.897,37.897,37.897,37.896,37.895,37.895,37.894,&  !1519
  40.737,40.319,40.047,39.855,39.714,39.604,39.517,39.445,39.385,&  !1520
  40.581,40.164,39.893,39.702,39.561,39.452,39.366,39.295,39.235,&  !1676
  45.521,44.456,43.753,43.254,42.878,42.580,42.332,42.119,41.930,&  !1677
  45.520,44.455,43.752,43.251,42.871,42.569,42.315,42.094,41.896,&  !1978
  55.068,51.783,49.553,47.942,46.723,45.768,44.997,44.360,43.823,&  !1979
  53.868,50.369,48.031,46.355,45.092,44.104,43.308,42.652,42.100,&  !5379
  54.133,50.597,48.233,46.539,45.261,44.262,43.456,42.790,42.230,&  !5380
  54.051,50.514,48.150,46.454,45.176,44.175,43.368,42.702,42.141,&  !5624
  54.442,50.854,48.455,46.733,45.433,44.415,43.592,42.912,42.340,&  !5625
  54.320,50.722,48.313,46.583,45.277,44.251,43.423,42.738,42.160,&  !6260
  55.691,51.965,49.444,47.615,46.221,45.119,44.223,43.478,42.848,&  !6261
  55.661,51.933,49.412,47.582,46.188,45.085,44.189,43.445,42.813,&  !6349
  55.973,52.193,49.630,47.769,46.349,45.226,44.314,43.555,42.913,&  !6350
  55.922,52.141,49.577,47.715,46.295,45.172,44.259,43.500,42.858,&  !6491
  56.828,52.821,50.110,48.146,46.654,45.477,44.522,43.730,43.061,&  !6492
  56.657,52.653,49.944,47.983,46.491,45.315,44.360,43.569,42.901/  !6900
!    3P,1D,1S,1D,3D,3F,1D,3P
  DATA FREQSI/2.1413750E15,1.9723165E15,1.7879689E15,1.5152920E15,&
       5.5723927E14,5.3295914E14,4.7886458E14,4.7216422E14,4.6185133E14/
  DATA FLOG/35.45438,35.30022,35.21799,35.11986,34.95438,33.95402,&
       33.90947,33.80244,33.78835,33.76626,33.70518/
  DATA TLG/8.29405,8.51719,8.69951,8.85367,8.98720,9.10498,9.21034,&
       9.30565,9.39266/
   DATA FREQ1/0./
   save x

      DO K=1,NRHOX
       N=MAX(MIN(8, INT(T(K)/1000.)-3),1)
       NT(K)=N
       DT(K)=(TLOG(K)-TLG(N))/(TLG(N+1)-TLG(N))
       enddo

   IF(FREQ /= FREQ1)then
    FREQ1=FREQ
      DO  N=1,9
       IF(FREQ > FREQSI(N))exit
      enddo
      N=n
     D=(FREQLG-FLOG(N))/(FLOG(N+1)-FLOG(N))
      IF(N > 2)N=2*N-2
      D1=1.d0-D
      DO IT=1,9
       X(IT)=PEACH(IT,N+1)*D+PEACH(IT,N)*D1
      enddo
   endif
    N=NT(J)
      SI1OP=EXP(-(X(N)*(1.-DT(J))+X(N+1)*DT(J)))*9.d0
  END
!-----
 real(8) FUNCTION FE1OP(J)
  use atmosphere_model
  use radiation_field
  use temporary_arrays
  implicit real(8) (a-h,o-z)
!     CROSS-SECTION TIMES PARTITION FUNCTION
      DIMENSION G(48),E(48),WNO(48),XSECT(48)
      DATA G/25.,35.,21.,15.,9.,35.,33.,21.,27.,49.,9.,21.,27.,9.,9.,&
       25.,33.,15.,35.,3.,5.,11.,15.,13.,15.,9.,21.,15.,21.,25.,35.,&
       9.,5.,45.,27.,21.,15.,21.,15.,25.,21.,35.,5.,15.,45.,35.,55.,25./
      DATA E/500.,7500.,12500.,17500.,19000.,19500.,19500.,21000.,&
       22000.,23000.,23000.,24000.,24000.,24500.,24500.,26000.,26500.,&
       26500.,27000.,27500.,28500.,29000.,29500.,29500.,29500.,30000.,&
       31500.,31500.,33500.,33500.,34000.,34500.,34500.,35000.,35500.,&
       37000.,37000.,37000.,38500.,40000.,40000.,41000.,41000.,43000.,&
       43000.,43000.,43000.,44000./
      DATA WNO/63500.,58500.,53500.,59500.,45000.,44500.,44500.,43000.,&
       58000.,41000.,54000.,40000.,40000.,57500.,55500.,38000.,57500.,&
       57500.,37000.,54500.,53500.,55000.,34500.,34500.,34500.,34000.,&
       32500.,32500.,32500.,32500.,32000.,29500.,29500.,31000.,30500.,&
       29000.,27000.,54000.,27500.,24000.,47000.,23000.,44000.,42000.,&
       42000.,21000.,42000.,42000./
      DATA ITEMP1/0/,FREQ1/0./
      save xsect
      IF(ITEMP1 == 0)then
       ITEMP1=1
        DO K=1,NRHOX
         DO I=1,48
          BOLT2(I,K)=G(I)*EXP(-E(I)*2.99792458E10*HKT(K))
         enddo
        enddo
      endif
    IF(FREQ /= FREQ1)then
     FREQ1=FREQ
     WAVENO=FREQ/2.99792458d10
      IF(WAVENO >= 21000.d0)then
       DO I=1,48
        XSECT(I)=0.d0
         IF(WNO(I) < WAVENO)XSECT(I)= 3.d-18/&
          (1.d0+((WNO(I)+3000.d0-WAVENO)/WNO(I)/.1d0)**4)
       enddo
      endif
    endif
   FE1OP=0.d0
      IF(WAVENO < 21000.)RETURN
       DO I=1,48
        FE1OP=FE1OP+XSECT(I)*BOLT2(I,J)
       enddo
      END
!---
  SUBROUTINE LUKEOP
  use atmosphere_model
  use radiation_field
  implicit real(8) (a-h,o-z)
!     SI2,MG2,CA2,N1,O1
  REAL(8) N1OP,MG2OP
    do j=1,nRhox
     aCont(j)=aCont(j)+(N1OP(J)*xDens(J,7,1)+O1OP(J)*xDens(J,8,1)+MG2OP(J)*&
      xDens(J,12,2)+SI2OP(J)*xDens(J,14,2)+CA2OP(J)*xDens(J,20,2))*&
        STIM(J)/1.49725d-2
     enddo
  END
!-----
 real(8) FUNCTION N1OP(J)
  use atmosphere_model
  use radiation_field
  use temporary_arrays
  implicit real(8) (a-h,o-z)
!     CROSS-SECTION TIMES PARTITION FUNCTION

   DATA FREQ1,ITEMP1/0.,0/
   save x1130,x1020,x853
      IF(ITEMP1 == 0)then
        ITEMP1=1
       DO  K=1,NRHOX
        C1130(K)=6.d0*EXP(-3.575d0/TKEV(K))
        C1020(K)=10.d0*EXP(-2.384d0/TKEV(K))
      enddo
     endif
   IF(FREQ /= FREQ1)then
      X1130=0.d0
      X1020=0.d0
      X853=0.d0
      IF(FREQ >= 3.517915d15)X853=SEATON(3.517915D15,1.142D-17,2.D0,4.29D0)
      IF(FREQ >= 2.941534d15)X1020=SEATON(2.941534D15,4.41D-18,1.5D0,3.85D0)
      IF(FREQ >= 2.653317E15)X1130=SEATON(2.653317D15,4.2D-18,1.5D0,4.34D0)
      FREQ1=FREQ
   endif
    N1OP=X853*4.d0+X1020*C1020(J)+X1130*C1130(J)
  END
!-----
  real(8) FUNCTION O1OP(J)
  use atmosphere_model
  use radiation_field
  implicit real(8) (a-h,o-z)
!     FROM DEANE PETERSON AFTER PEACH
!     CROSS-SECTION TIMES PARTITION FUNCTION
    DATA FREQ1/0./
    save x911
      IF(FREQ /= FREQ1)then
       X911=0.d0
        IF(FREQ >= 3.28805d15)X911=SEATON(3.28805D15,2.94D-18,1.D0,2.66D0)
       FREQ1=FREQ
      endif
     O1OP=X911*9.d0
   END
!----
 real(8) FUNCTION MG2OP(J)
  use atmosphere_model
  use radiation_field
  use temporary_arrays
  implicit real(8) (a-h,o-z)
!     CROSS-SECTION TIMES PARTITION FUNCTION

   DATA FREQ1,ITEMP1/0.,0/
   save x824,x1169
      IF(ITEMP1 == 0)then
       ITEMP1=1
        DO K=1,NRHOX
         C1169(K)=6.d0*EXP(-4.43d0/TKEV(K))
        enddo
      endif
   IF(FREQ /= FREQ1)then
      X1169=0.d0
      X824=0.d0
       IF(FREQ >= 3.635492d15)X824=SEATON(3.635492D15,1.40D-19,4.D0,6.7D0)
       IF(FREQ >= 2.564306E15)X1169=5.11d-19*(2.564306E15/FREQ)**3
      FREQ1=FREQ
   endif
    MG2OP=X824*2.d0+X1169*C1169(J)
 END
!----
 real(8) FUNCTION SI2OP(J)
  use atmosphere_model
  use radiation_field
  use temporary_arrays
  implicit real(8) (a-h,o-z)
!     CROSS-SECTION TIMES THE PARTITION FUNCTION
   DIMENSION FLOG(9),FREQSI(7),PEACH(6,14),X(6),TLG(6)
   DATA PEACH/&
!   10000     12000     14000     16000     18000     20000       WAVE(A)
  -43.8941, -43.8941, -43.8941, -43.8941, -43.8941, -43.8941,&      !500
  -42.2444, -42.2444, -42.2444, -42.2444, -42.2444, -42.2444,&      !600
  -40.6054, -40.6054, -40.6054, -40.6054, -40.6054, -40.6054,&      !759
  -54.2389, -52.2906, -50.8799, -49.8033, -48.9485, -48.2490,&      !760
  -50.4108, -48.4892, -47.1090, -46.0672, -45.2510, -44.5933,&     !1905
  -52.0936, -50.0741, -48.5999, -47.4676, -46.5649, -45.8246,&     !1906
  -51.9548, -49.9371, -48.4647, -47.3340, -46.4333, -45.6947,&     !1975
  -54.2407, -51.7319, -49.9178, -48.5395, -47.4529, -46.5709,&     !1976
  -52.7355, -50.2218, -48.4059, -47.0267, -45.9402, -45.0592,&     !3245
  -53.5387, -50.9189, -49.0200, -47.5750, -46.4341, -45.5082,&     !3246
  -53.2417, -50.6234, -48.7252, -47.2810, -46.1410, -45.2153,&     !3576
  -53.5097, -50.8535, -48.9263, -47.4586, -46.2994, -45.3581,&     !3577
  -54.0561, -51.2365, -49.1980, -47.6497, -46.4302, -45.4414,&     !3900
  -53.8469, -51.0256, -48.9860, -47.4368, -46.2162, -45.2266/     !4200
   DATA FREQSI/4.9965417d15,3.9466738d15,1.5736321d15,1.5171539d15,&
      9.2378947d14,8.3825004d14,7.6869872d14/
!     2P,2D,2P,2D,2P
   DATA FLOG/36.32984,36.14752,35.91165,34.99216,34.95561,34.45951,&
       34.36234,34.27572,34.20161/
   DATA TLG/9.21034,9.39266,9.54681,9.68034,9.79813,9.90349/
   DATA FREQ1/0./
   save x

        DO K=1,NRHOX
         N=MAX(MIN(5, INT(T(K)/2000.d0)-4),1)
         NT(K)=N
         DT(K)=(TLOG(K)-TLG(N))/(TLG(N+1)-TLG(N))
        enddo

  IF(FREQ /= FREQ1)then
   FREQ1=FREQ
    DO  N=1,7
      IF(FREQ > FREQSI(N))exit
    enddo
     n=n
     D=(FREQLG-FLOG(N))/(FLOG(N+1)-FLOG(N))
      IF(N > 2)N=2*N-2
      IF(N == 14)N=13
      D1=1.d0-D
       DO IT=1,6
        X(IT)=PEACH(IT,N+1)*D+PEACH(IT,N)*D1
       enddo
   endif
    N=NT(J)
     SI2OP=EXP(X(N)*(1.d0-DT(J))+X(N+1)*DT(J))*6.d0
  END
!----
 real(8) FUNCTION CA2OP(J)
  use atmosphere_model
  use radiation_field
  use temporary_arrays
  implicit real(8) (a-h,o-z)
!     CROSS-SECTION TIMES THE PARTITION FUNCTION
   DATA FREQ1,ITEMP1/0.,0/
   save x120,x1218,x1044
      IF(ITEMP1 == 0)then
       ITEMP1=1
        DO K=1,NRHOX
         C1218(K)=10.d0*EXP(-1.697d0/TKEV(K))
         C1420(K)=6.d0*EXP(-3.142d0/TKEV(K))
        enddo
      endif
   IF(FREQ /= FREQ1)then
      X1420=0.d0
      X1218=0.d0
      X1044=0.d0
       IF(FREQ >= 2.870454d15)X1044=5.4d-20*(2.870454d15/FREQ)**3
       IF(FREQ >= 2.460127d15)X1218=1.64d-17*SQRT(2.460127d15/FREQ)
       IF(FREQ >= 2.110779d15)X1420=SEATON(2.110779D15,4.13D-18,3.D0,.69D0)
      FREQ1=FREQ
    endif
     CA2OP=X1044*2.d0+X1218*C1218(J)+X1420*C1420(J)
  END
!---
  SUBROUTINE HOTOP
  use atmosphere_model
  use radiation_field
  use temporary_arrays
  implicit real(8) (a-h,o-z)
   DIMENSION A(420)
   DIMENSION A1(63),A2(63),A3(63),A4(63),A5(63),A6(63),A7(42)
    EQUIVALENCE (A(1),A1(1)),(A(64),A2(1)),(A(127),A3(1))
    EQUIVALENCE (A(190),A4(1)),(A(253),A5(1)),(A(316),A6(1))
    EQUIVALENCE (A(379),A7(1))
   DATA A1/&
    4.149945d15, 6.90d-18, 1.000, 6.,  6.,  13.71,  2.,&              !6.01
    4.574341d15, 2.50d-18, 1.000, 4.,  2.,  11.96,  2.,&              !6.01
    5.220770d15, 1.08d-17, 1.000, 4., 10.,   9.28,  2.,&              !6.01
    5.222307d15, 5.35d-18, 3.769, 2.,  1.,   0.00, 16.,&              !10.00
    5.892577d15, 4.60d-18, 1.950, 6.,  6.,   0.00,  2.,&              !6.01
    6.177022d15, 3.50d-18, 1.000, 4., 12.,   5.33,  2.,&              !6.01
    6.181062d15, 6.75d-18, 3.101, 5.,  1.,   4.05,  6.,&              !7.01
    6.701879d15, 6.65d-18, 2.789, 5.,  5.,   1.90,  6.,&              !7.01
    7.158382d15, 6.65d-18, 2.860, 6.,  9.,   0.00,  6./               !7.01
   DATA A2/&
    7.284488d15, 3.43d-18, 4.174, 5.,  6.,   5.02, 11.,&              !8.01
    7.693612d15, 3.53d-18, 3.808, 5., 10.,   3.33, 11.,&              !8.01
    7.885955d15, 2.32d-18, 3.110, 5.,  6.,   5.02, 11.,&              !8.01
    8.295079d15, 3.97d-18, 3.033, 5., 10.,   3.33, 11.,&              !8.01
    8.497686d15, 7.32d-18, 3.837, 5.,  4.,   0.00, 11.,&              !8.01
    8.509966d15, 2.00d-18, 1.750, 7.,  3.,  12.69,  3.,&              !6.02
    8.572854d15, 1.68d-18, 3.751, 5.,  6.,   5.02, 11.,&              !8.01
    9.906370d15, 4.16d-18, 2.717, 3.,  6.,   0.00, 17.,&              !10.01
    1.000693d16, 2.40d-18, 1.750, 7.,  9.,   6.50,  3./              !6.02
   DATA A3/&
    1.046078d16, 4.80d-18, 1.000, 4., 10.,  12.53,  7.,&              !7.02
    1.067157d16, 2.71d-18, 2.148, 3.,  6.,   0.00, 17.,&              !10.01
    1.146734d16, 2.06d-18, 1.626, 6.,  6.,   0.00,  7.,&              !7.02
    1.156813d16, 5.20d-19, 2.126, 3.,  6.,   0.00, 17.,&              !10.01
    1.157840d16, 9.10d-19, 4.750, 4.,  1.,   0.00,  3.,&              !6.02
    1.177220d16, 5.30d-18, 1.000, 4., 12.,   7.10,  7.,&              !7.02
    1.198813d16, 3.97d-18, 2.780, 6.,  1.,   5.35, 12.,&              !8.02
    1.325920d16, 3.79d-18, 2.777, 6.,  5.,   2.51, 12.,&              !8.02
    1.327649d16, 3.65d-18, 2.014, 6.,  9.,   0.00, 12./               !8.02
   DATA A4/&
    1.361466d16, 7.00d-18, 1.000, 2.,  5.,   7.48, 12.,&              !8.02
    1.365932d16, 9.30d-19, 1.500, 7.,  6.,   8.00,  4.,&              !6.03
    1.481487d16, 1.10d-18, 1.750, 7.,  3.,  16.20,  8.,&              !7.03
    1.490032d16, 5.49d-18, 3.000, 5.,  1.,   6.91, 18.,&              !10.02
    1.533389d16, 1.80d-18, 2.277, 4.,  9.,   0.00, 18.,&              !10.02
    1.559452d16, 8.70d-19, 3.000, 6.,  2.,   0.00,  4.,&              !6.03
    1.579688d16, 4.17d-18, 2.074, 4.,  5.,   3.20, 18.,&              !10.02
    1.643205d16, 1.39d-18, 2.792, 5.,  5.,   3.20, 18.,&              !10.02
    1.656208d16, 2.50d-18, 2.346, 5.,  9.,   0.00, 18./              !10.02
   DATA A5/&
    1.671401d16, 1.30d-18, 1.750, 7.,  9.,   8.35,  8.,&              !7.03
    1.719725d16, 1.48d-18, 2.225, 5.,  9.,   0.00, 18.,&              !10.02
    1.737839d16, 2.70d-18, 1.000, 4., 10.,  15.74, 13.,&              !8.03
    1.871079d16, 1.27d-18,  .831, 6.,  6.,   0.00, 13.,&              !8.03
    1.873298d16, 9.10d-19, 3.000, 4.,  1.,   0.00,  8.,&              !7.03
    1.903597d16, 2.90d-18, 1.000, 4., 12.,   8.88, 13.,&              !8.03
    2.060738d16, 4.60d-18, 1.000, 3., 12.,  22.84, 19.,&              !10.03
    2.125492d16, 5.90d-19, 1.000, 6.,  6.,   9.99,  9.,&              !7.04
    2.162610d16, 1.69d-18, 1.937, 5.,  6.,   7.71, 19./              !10.03
   DATA A6/&
    2.226127d16, 1.69d-18, 1.841, 5., 10.,   5.08, 19.,&              !10.03
    2.251163d16, 9.30d-19, 2.455, 6.,  6.,   7.71, 19.,&              !10.03
    2.278001d16, 7.90d-19, 1.000, 6.,  9.,  10.20, 14.,&              !8.04
    2.317678d16, 1.65d-18, 2.277, 6., 10.,   5.08, 19.,&              !10.03
    2.348946d16, 3.11d-18, 1.963, 6.,  4.,   0.00, 19.,&              !10.03
    2.351911d16, 7.30d-19, 1.486, 5.,  6.,   7.71, 19.,&              !10.03
    2.366973d16, 5.00d-19, 1.000, 4.,  2.,   0.00,  9.,&              !7.04
    2.507544d16, 6.90d-19, 1.000, 6.,  3.,  19.69, 14.,&              !8.04
    2.754065d16, 7.60d-19, 1.000, 2.,  1.,   0.00, 14./               !8.04
   DATA A7/&
    2.864850d16, 1.54d-18, 2.104, 6.,  1.,   7.92, 20.,&              !10.04
    2.965598d16, 1.53d-18, 2.021, 6.,  5.,   3.76, 20.,&              !10.04
    3.054151d16, 1.40d-18, 1.471, 6.,  9.,   0.00, 20.,&              !10.04
    3.085141d16, 2.80d-18, 1.000, 4.,  5.,  11.01, 20.,&              !10.04
    3.339687d16, 3.60d-19, 1.000, 6.,  2.,   0.00, 15.,&              !8.05
    3.818757d16, 4.90d-19, 1.145, 6.,  6.,   0.00, 21./               !10.05
   DATA NUM/60/
!      IF(FREQ.LT.A(1))RETURN
!     FREE-FREE
  DO J=1,NRHOX
   FREE=COULFF(J,1)*1.d0**2*(xDens2(J,6,2)+xDens2(J,7,2)+xDens2(J,8,2)+&
     xDens2(J,10,2)+xDens2(J,12,2)+xDens2(J,14,2)+xDens2(J,16,2)+xDens2(J,26,2))+&
        COULFF(J,2)*2.d0**2*(xDens2(J,6,3)+xDens2(J,7,3)+xDens2(J,8,3)+&
     xDens2(J,10,3)+xDens2(J,12,3)+xDens2(J,14,3)+xDens2(J,16,3)+xDens2(J,26,3))+&
        COULFF(J,3)*3.d0**2*(xDens2(J,6,4)+xDens2(J,7,4)+xDens2(J,8,4)+&
     xDens2(J,10,4)+xDens2(J,12,4)+xDens2(J,14,4)+xDens2(J,16,4)+xDens2(J,26,4))+&
        COULFF(J,4)*4.d0**2*(xDens2(J,6,5)+xDens2(J,7,5)+xDens2(J,8,5)+&
     xDens2(J,10,5)+xDens2(J,12,5)+xDens2(J,14,5)+xDens2(J,16,5)+xDens2(J,26,5))+&
        COULFF(J,5)*5.d0**2*(xDens2(J,6,6)+xDens2(J,7,6)+xDens2(J,8,6)+&
     xDens(J,10,6)+xDens2(J,12,6)+xDens2(J,14,6)+xDens2(J,16,6))
   ahot(J)=FREE*3.6919d8/FREQ**3*XNE(J)/SQRT(T(J))/1.49725d-2
  enddo
   L=-6
    DO I=1,NUM
     L=L+7
     IF(FREQ < A(L))cycle
      XSECT=A(L+1)*(A(L+2)+(A(L)/FREQ)-A(L+2)*(A(L)/FREQ))*&
        SQRT((A(L)/FREQ)** INT(A(L+3)))
      ID=A(L+6)
      if (ID <= 4) then
           m = ID
           k = 6
      else
        if (ID <= 9) then
           m = ID - (21 - 16) + 1
           k = 7
        else
          if (ID <= 15) then
           m = ID - (21 - 15) - 3
           k = 8
          else
           m = ID - (21 - 15) - 9
           k = 10
          endif
        endif
      endif
    DO J=1,NRHOX
      XX=XSECT*xDens(J,k,m)*A(L+4)/1.49725d-2
      IF(XX > AHOT(J)/100.d0)AHOT(J)=AHOT(J)+XX/EXP(A(L+5)/TKEV(J))
    enddo
   enddo
      DO J=1,NRHOX
       AHOT(J)=AHOT(J)*STIM(J)
       Acont(j)=Acont(j)+Ahot(j)
      enddo
 END
!---
  SUBROUTINE ELECOP
  use atmosphere_model
  use radiation_field
  implicit real(8) (a-h,o-z)
     DO J=1,NRHOX
      sigmaC(j)=sigmaC(j)+0.6653d-24*XNE(J)/RHO(J)
     enddo
  END
