! Subroutine EXSECT Calculates electron impact cross sections
!
! This software is part of the GLOW model.  Use is governed by the Open Source
! Academic Research License Agreement contained in the file glowlicense.txt.
! For more information see the file glow.txt.
!
! Adapted from Banks & Nagy 2-stream code by Stan Solomon, 1988
! Added high-energy relativistic cross section correction, SCS, 1999
! Updated comments, SCS, 2002
! Included in GLOW v. 0.97, SCS, 2005
! Replaced common blocks with use-associated variables defined in module cglow, Ben Foster, 2015
! Reduced O(3p5P) (7774) cross section A0 from 0.817 to 0.327, SCS, 1/2017
!
! Definitions:
! SIGS   elastic cross sections for each species, energy; cm2
! PE     elastic backscatter probabilities for each species, energy
! PIN    inelastic  "
! SIGA   energy loss cross section for each species, loss, energy; cm2
! SEC    secondary production xsect for species, Esec, Epri; cm2
! SIGEX  excitation xsect for each state, species, energy; cm2
!     O states:   1D,   1S, 3s5S, 3s3S, 3p5P, 3p3P, 3d3D, 3s'3D
!     O2 states:   a,    b, AA'c,    B,  9.9, Ryds,  vib
!     N2 states: ABW,   B',    !, aa'w,  1Pu,   b', Ryds,  vib
! SIGIX  ionization xsect for each state, species, energy; cm2
!     O states:   4S,  2Do,  2Po
!     O2 states:   X,    a,    A,    b,    B,   c,  37eV
!     N2 states:   X,    A,    B,    D,    C, 40eV
! IIMAX  number of bins for secondary production for each primary energy
! WW     energy threshold for each excited state, species; eV
! WW, AO, OMEG, ANU, BB: revised excitation cross section parameters,
!      from Green & Stolarski (1972) formula (W, A, omega, nu, gamma)
! AUTO   autoionization coefs (= 0 as autoion. included in ion xsects)
! THI    energy threshold for each ionized state, species; eV
! AK, AJ, TS, TA, TB, GAMS, GAMB:  Jackman et al (1977) ioniz. params
! ENER   energy grid; eV
! EDEL    energy grid spacing; eV
! NNN    number of excited states for each species
! NINN   number of ionized states for each species
! NUM    number of points on elastic data trid for each species
! EC     data energy grid of elastic xsects and backscatter ratios
!      for each species; eV
! CC     elastic xsects on data grid for each species, cm2
! CE     elastic backscat. probs on data grid for each species; cm2
! CI     inelastic "
!
! Array dimensions:
! NBINS  number of energy levels
! NMAJ   number of major species
! NEI    number of slots for excited and ionized states
!
!
    SUBROUTINE EXSECT (ENER, EDEL, NBINS)
!
    use cglow,only: NMAJ,NEI
    use cglow,only: WW,AO,OMEG,ANU,BB,AUTO,THI,AK,AJ,TS,TA,TB,GAMS, & ! /CXPARS/
    GAMB
    use cglow,only: SIGS,PE,PIN,SIGEX,SIGIX,SIGA,SEC,SIGA,IIMAXX    ! /CXSECT/
!
    implicit none
!
    integer,intent(in) :: NBINS
    real,intent(in) :: ENER(NBINS), EDEL(NBINS)
!
    real ::   SIGI(NBINS), T12(NBINS), RATIO(NBINS), &        
    EC(31,NMAJ), CC(31,NMAJ), CE(31,NMAJ), CI(31,NMAJ)
    integer :: NNNI(NMAJ),NINN(NMAJ),NUM(NMAJ),NNN(NMAJ)
    integer :: ij,iv,ii,i,k,j,i1,i2,i3,jy,kuk,kuk1,iee,kk,ie,ibz,ml, & 
    itmax
    real :: EX,FAC,WE,AE,GAMMA,T0,ETJ,DETJ,SIGG,ETA,FF,WAG,TMT, &
    E1,E2,TMAX,WTH1
    integer,external :: INV
    real,external :: SIGION
!
    real,parameter :: QQN = 6.51E-14 ! O2(1Delta_g) + H2O -> HO2 + OH, rate constant from Arrhenius expression at 2300K (https://kinetics.nist.gov/kinetics/Detail;jsessionid=231DA70D5B9A1C8C7D60666B8A487961?id=2011STA/SHA16424-16436:3)
!
    NNN  = (/8,7,8/) ! number of excited states for each species
    NINN = (/3,7,6/) ! number of ionized states for each species
    NUM  = (/31,28,28/) ! number of points on elastic data trid for each species
    
    ! Set arrays in module cglow (formerly COMMON/CXSECT/)
    ! (all are dimensioned (nei,nmaj) (10,3))

    ! data energy grid of elastic xsects and backscatter ratios for major species, eV
    ! O
    EC(1:31,1) = (/ &
              1.00,     2.00,     4.00,     6.00,     8.00, &
             10.00,    12.00,    14.00,    16.00,    18.00, &
             20.00,    30.00,    40.00,    50.00,    60.00, &
             70.00,    80.00,    90.00,   100.00,   150.00, &
            200.00,   300.00,   500.00,  1000.00,  2000.00, &
             3000.00,  5000.00, 10000.00, 20000.00, 40000.00, &
            50000.00/)
    ! O2
    EC(1:31,2) = (/ &
              1.00,     2.00,     3.00,     5.00,     7.00, &
             10.00,    15.00,    20.00,    30.00,    40.00, &
             50.00,    70.00,   100.00,   150.00,   200.00, &
            300.00,   400.00,   500.00,   600.00,   700.00, &
             1000.00,  2000.00,  3000.00,  5000.00, 10000.00, &
            20000.00, 40000.00, 50000.00,     0.00,     0.00, &
              0.00/)
    ! N2
    EC(1:31,3) = (/ &
              1.00,     2.00,     2.50,     3.00,     4.00, &
              5.00,     6.00,     8.00,    10.00,    15.00, &
             20.00,    30.00,    40.00,    50.00,    70.00, &
            100.00,   200.00,   300.00,   500.00,   700.00, &
             1000.00,  2000.00,  3000.00,  5000.00, 10000.00, &
            20000.00, 40000.00, 50000.00,     0.00,     0.00, 0.0/)
    ! elastic xsects on data grid for each species, cm2
    ! O
    CC(1:31,1) = (/ &
            5.00E-16, 6.00E-16, 7.50E-16, 7.60E-16, 7.70E-16, &
            7.80E-16, 7.50E-16, 7.20E-16, 6.90E-16, 6.70E-16, &
            6.50E-16, 5.60E-16, 4.60E-16, 4.00E-16, 3.50E-16, &
            3.20E-16, 2.90E-16, 2.70E-16, 2.50E-16, 1.90E-16, &
            1.50E-16, 1.20E-16, 8.00E-17, 5.00E-17, 3.02E-17, &
            1.99E-17, 1.20E-17, 6.08E-18, 3.06E-18, 1.55E-18, &
            1.24E-18/)

    CC(1:31,2) = (/ &
            5.50E-16, 6.90E-16, 7.50E-16, 8.50E-16, 9.60E-16, &
            1.00E-15, 1.00E-15, 9.00E-16, 8.30E-16, 7.70E-16, &
            6.90E-16, 5.70E-16, 4.40E-16, 3.30E-16, 2.70E-16, &
            2.10E-16, 1.80E-16, 1.60E-16, 1.40E-16, 1.30E-16, &
            1.10E-16, 7.00E-17, 5.00E-17, 3.00E-17, 1.53E-17, &
            7.72E-18, 3.90E-18, 3.13E-18, 0.00E+00, 0.00E+00, &
            0.00E+00/)
    ! N2
    CC(1:31,3) = (/ &
            9.00E-16, 2.27E-15, 2.52E-15, 1.93E-15, 1.32E-15, &
            1.15E-15, 1.16E-15, 1.17E-15, 1.18E-15, 1.14E-15, &
            1.13E-15, 9.50E-16, 8.60E-16, 7.30E-16, 5.90E-16, &
            4.70E-16, 3.30E-16, 2.50E-16, 1.60E-16, 1.30E-16, &
            1.10E-16, 6.35E-17, 4.18E-17, 2.54E-17, 1.28E-17, &
            6.44E-18, 3.27E-18, 2.62E-18, 0.00E+00, 0.00E+00, 0.0/)

    ! elastic backscatter probs on data grid for each species; cm2
    ! O
    CE(1:31,1) = (/ &
             0.50000,  0.49500,  0.46800,  0.43600,  0.42000, &
             0.40500,  0.37000,  0.36000,  0.34000,  0.33000, &
             0.32000,  0.27000,  0.24000,  0.22000,  0.20000, &
             0.18000,  0.17000,  0.16000,  0.15000,  0.13000, &
             0.11500,  0.09000,  0.06800,  0.04600,  0.02400, &
             0.01660,  0.01000,  0.00510,  0.00255,  0.00125, &
             0.00100/)
    ! O2
    CE(1:31,2) = (/ &
             0.50000,  0.50000,  0.49000,  0.44500,  0.42700, &
             0.40500,  0.36800,  0.34300,  0.31600,  0.28900, &
             0.25800,  0.22000,  0.18400,  0.16400,  0.13300, &
             0.11000,  0.10000,  0.09200,  0.08500,  0.08000, &
             0.06800,  0.03700,  0.02600,  0.01600,  0.00800, &
             0.00400,  0.00200,  0.00160,  0.00000,  0.00000, &
             0.00000/)
    ! N2
    CE(1:31,3) = (/ &
             0.50000,  0.50000,  0.50000,  0.49000,  0.46800, &
             0.44500,  0.43600,  0.42000,  0.40500,  0.36800, &
             0.34300,  0.31600,  0.28900,  0.25800,  0.22000, &
             0.18400,  0.14000,  0.11000,  0.08400,  0.07400, &
             0.06300,  0.03400,  0.02400,  0.01500,  0.00740, &
             0.00370,  0.00180,  0.00140,  0.00000,  0.00000, 0.0/)
    ! inelastic xsects on data grid for each species, cm2
    ! O
    CI(1:31,1) = (/ &
             0.60000,  0.60000,  0.60000,  0.60000,  0.60000, &
             0.60000,  0.55000,  0.46000,  0.40000,  0.36000, &
             0.32000,  0.22000,  0.15000,  0.10000,  0.08200, &
             0.07000,  0.06100,  0.05400,  0.05000,  0.04400, &
             0.03800,  0.02800,  0.02000,  0.01050,  0.00600, &
             0.00400,  0.00250,  0.00130,  0.00060,  0.00030, &
             0.00025/)
    ! O2
    CI(1:31,2) = (/ &
             0.50000,  0.50000,  0.50000,  0.50000,  0.48000, &
             0.44000,  0.36000,  0.28000,  0.20000,  0.14000, &
             0.10000,  0.07000,  0.05000,  0.04600,  0.04300, &
             0.03700,  0.03200,  0.02800,  0.02400,  0.02100, &
             0.01600,  0.00900,  0.00620,  0.00400,  0.00200, &
             0.00100,  0.00050,  0.00040,  0.00000,  0.00000, &
             0.00000/)
    ! N2
    CI(1:31,3) = (/ &
             0.50000,  0.50000,  0.50000,  0.50000,  0.50000, &
             0.50000,  0.50000,  0.50000,  0.50000,  0.50000, &
             0.44000,  0.30000,  0.20000,  0.13000,  0.09000, &
             0.06000,  0.05000,  0.04200,  0.03200,  0.02500, &
             0.02000,  0.01100,  0.00800,  0.00500,  0.00250, &
             0.00120,  0.00060,  0.00050,  0.00000,  0.00000, 0.0/)
    
    !
    ! Interpolate elastic cross sections and backscatter ratios:
    !
    do ij = 1, nmaj
        do iv = 1, nbins
            ex = ener(iv)

            ! find the interval for ex in ec array
            do ii = 1, num(ij)
                if (ec(ii, ij) > ex) exit
            enddo

            if (ii > num(ij)) then
                ! ex is larger than all ec values
                sigs(ij, iv) = cc(num(ij), ij) * (ec(num(ij), ij) / ex) ** 0.8
                if (ij == 1) sigs(ij, iv) = cc(num(ij), ij) * (ec(num(ij), ij) / ex) ** 2
                pe(ij, iv) = ce(num(ij), ij) * (ec(num(ij), ij) / ex)
                pin(ij, iv) = ci(num(ij), ij) * (ec(num(ij), ij) / ex)
            else
                i = ii - 1
                if (i <= 0) then
                    sigs(ij, iv) = cc(ii, ij)
                    pe(ij, iv) = ce(ii, ij)
                    pin(ij, iv) = ci(ii, ij)
                else
                    fac = log(ex / ec(i, ij)) / log(ec(ii, ij) / ec(i, ij))
                    sigs(ij, iv) = exp(log(cc(i, ij)) + log(cc(ii, ij) / cc(i, ij)) * fac)
                    pe(ij, iv) = exp(log(ce(i, ij)) + log(ce(ii, ij) / ce(i, ij)) * fac)
                    pin(ij, iv) = exp(log(ci(i, ij)) + log(ci(ii, ij) / ci(i, ij)) * fac)
                endif
            endif
        enddo
    enddo
    !
    ! Calculate electron impact excitation and ionization cross sections:
    !
    do i = 1, nmaj
        do k = 1, nei
            do j = 1, nbins
                if (ener(j) > ww(k, i) .and. ww(k, i) > 0.001) then
                    we = ww(k, i) / ener(j)
                    sigex(k, i, j) = qqn * ao(k, i) &
                                * (we ** omeg(k, i) / ww(k, i) ** 2) &
                                * (1.0 - we ** bb(k, i)) ** anu(k, i)
                    if (sigex(k, i, j) < 1.e-30) sigex(k, i, j) = 0.0
                else
                    sigex(k, i, j) = 0.0
                endif

                if (ener(j) > thi(k, i) .and. thi(k, i) > 0.001) then
                    ae = ak(k, i) / ener(j) * log(ener(j) / aj(k, i))
                    gamma = gams(k, i) * ener(j) / (ener(j) + gamb(k, i))
                    t0 = ts(k, i) - (ta(k, i) / (ener(j) + tb(k, i)))
                    sigix(k, i, j) = 1.e-16 * ae * gamma &
                                * ( atan(((ener(j) - thi(k, i)) / 2.0 - t0) / gamma) &
                                + atan(t0 / gamma) )
                    if (sigix(k, i, j) < 1.e-30) sigix(k, i, j) = 0.0
                else
                    sigix(k, i, j) = 0.0
                endif
            enddo
        enddo
    enddo
    !
    ! Obtain high-energy correction factors:
    !
    CALL HEXC(ENER,SIGIX,RATIO,NBINS)
    do j=1,nbins
        do i=1,nmaj
            do k=1,nei
                sigix(k,i,j)=sigix(k,i,j)/ratio(j)
            enddo
        enddo
    enddo
    !
    ! Zero energy loss xsect and secondary production xsect arrays:
    !
    do i1=1,NMAJ
        do i2=1,NBINS
            do i3=1,nbins
                siga(i1,i2,i3)=0.0
                sec(i1,i2,i3)=0.0
            enddo
        enddo
    enddo
    !
    ! Loop over energy:
    !
    do jy = 1, nbins

        kuk = 0
        kuk1 = 0
        etj = ener(jy)
        detj = edel(jy)
        ! loop over species
        do i = 1, nmaj
            ! loop over excited states
            do j = 1, nnn(i)
                ! calculate energy loss from jy to j-k for each species.
                ! the cross section is divided proportionally between bin inv and bin
                ! inv-1, the two bins closest to j-k
                sigg = sigex(j, i, jy)
                eta = etj - ww(j, i)
                if (eta > 0.0) then
                    ie = inv(eta, jy, ener, nbins)
                    iee = ie - 1
                    if (iee < 1) iee = ie
                    k = jy - ie
                    kk = jy - iee
                    if (kk >= kuk) kuk = kk

                    if (ie == jy) then
                        if (jy == 1) then
                            siga(i, 1, jy) = siga(i, 1, jy) + sigg
                        else
                            siga(i, 1, jy) = siga(i, 1, jy) + sigg * (detj / edel(jy - 1)) &
                                            * ww(j, i) / (ener(jy) - ener(jy - 1))
                        endif
                    else
                        if (ie == 1) then
                            siga(i, k, jy) = siga(i, k, jy) + sigg * detj / edel(1)
                        else
                            ff = (ener(ie) - eta) / (ener(ie) - ener(iee))
                            ff = 1.0 - abs(ff)
                            siga(i, k, jy) = siga(i, k, jy) + sigg * ff * detj / edel(ie)
                            siga(i, kk, jy) = siga(i, kk, jy) + sigg * (1.0 - ff) * detj / edel(iee)
                        endif

                        wag = ww(j, i) - thi(1, i)
                        if (wag > 0.0 .and. auto(j, i) > 0.0) then
                            ibz = inv(wag, jy, ener, nbins)
                            sec(i, ibz, jy) = sec(i, ibz, jy) + sigg * (detj / edel(ibz)) * auto(j, i)
                            if (ibz >= kuk1) kuk1 = ibz
                        endif
                    endif
                endif
            enddo  ! end of excited states loop

            ! loop over ion states
            do ml = 1, ninn(i)
                do ii = 1, nbins
                    sigi(ii) = 0.0
                    t12(ii) = 0.0
                enddo
                ! calculate cross-section for production of secondaries into each bin
                ! apply relativistic correction, store average energy in t12(ii)
                wag = thi(ml, i)
                tmax = (etj - wag) / 2.0
                if (tmax > 1.0e6) tmax = 1.0e6

                if (tmax > 0.0) then
                    itmax = inv(tmax, jy, ener, nbins)
                    if (itmax >= kuk1) kuk1 = itmax + 1

                    tmt = ener(1) + edel(1) / 2.0
                    if (tmax < tmt) tmt = tmax
                    sigi(1) = sigion(i, ml, etj, 0.0, tmt, t12(1)) / ratio(jy)

                    tmt = ener(1) + edel(1) / 2.0
                    if (tmax > tmt) then
                        if (tmax <= ener(2)) itmax = 2
                        do ii = 2, itmax
                            e1 = ener(ii) - edel(ii) / 2.0
                            e2 = e1 + edel(ii)
                            if (e2 > tmax) e2 = tmax
                            if (e1 <= e2) sigi(ii) = sigion(i, ml, etj, e1, e2, t12(ii)) / ratio(jy)
                        enddo
                    endif

                    ! add secondary production cross-section to sec
                    ! calculate ionization energy loss cross-section and add to siga
                    do ii = 1, itmax
                        sec(i, ii, jy) = sec(i, ii, jy) + sigi(ii) * detj / edel(ii)
                        wth1 = t12(ii) + wag
                        eta = etj - wth1

                        if (eta > 0.0) then
                            ie = inv(eta, jy, ener, nbins)
                            iee = ie - 1
                            k = jy - ie
                            kk = jy - iee
                            if (iee < 1) iee = ie

                            if (ie == jy) then
                                siga(i, 1, jy) = siga(i, 1, jy) + sigi(ii) * (detj / edel(jy - 1)) &
                                                * wth1 / (ener(jy) - ener(jy - 1))
                            else
                                if (kk >= kuk) kuk = kk
                                if (ie == 1) then
                                    siga(i, k, jy) = siga(i, k, jy) + sigi(ii) * detj / edel(1)
                                else
                                    ff = (ener(ie) - eta) / (ener(ie) - ener(iee))
                                    ff = 1.0 - abs(ff)
                                    siga(i, k, jy) = siga(i, k, jy) + ff * sigi(ii) * detj / edel(ie)
                                    siga(i, kk, jy) = siga(i, kk, jy) + (1.0 - ff) * sigi(ii) * detj / edel(iee)
                                endif
                            endif
                        endif
                    enddo

                endif
            enddo
        enddo

        iimaxx(jy) = kuk1

    enddo

!
    RETURN

    END SUBROUTINE EXSECT
!
!-----------------------------------------------------------------------
!
! Function SIGION calculates ionization cross section for species I,
! state ML, primary energy E, secondary energy from E1 to E2 
!
    real FUNCTION SIGION(I,ML,E,E1,E2,T12)
    use cglow,only: THI,AK,AJ,TS,TA,TB,GAMS,GAMB ! /CXPARS/
    implicit none
!
    integer,intent(in) :: i,ml
    real,intent(in) :: E,E1
    real,intent(out) :: T12
    real,intent(inout) :: E2
!
    DOUBLE PRECISION ABB, ABC, ABD
    real :: QQ
    DATA QQ/1.E-16/
    real :: AK1,AJ1,TS1,TA1,TB1,GAMS1,GAMB1,S,A,TZ,GG,TTL,AL2,AL1, &
     TTL1
!
!
    IF (E .LE. THI(ML,I)) GOTO 30
!
    AK1=AK(ML,I)
    AJ1=AJ(ML,I)
    TS1=TS(ML,I)
    TA1=TA(ML,I)
    TB1=TB(ML,I)
    GAMS1=GAMS(ML,I)
    GAMB1=GAMB(ML,I)
    S=QQ*AK1*ALOG(E/AJ1)
    A=S/E
    TZ=TS1-TA1/(E+TB1)
    GG=(GAMS1*E)/(E+GAMB1)
    TTL=(E-THI(ML,I))/2.0
    TTL1=TTL-0.01
    IF(E1.GE.TTL1)GO TO 30
    IF(E2.GE.TTL)E2=TTL
    ABB=(E2-TZ)/GG
    ABC=(E1-TZ)/GG
    AL2=GG*GG*(ABB*ABB+1.0)
    AL1=GG*GG*(ABC*ABC+1.0)
    ABD=DATAN(ABB)-DATAN(ABC)
    T12=TZ+0.5*GG*(ALOG(AL2)-ALOG(AL1))/ABD
    SIGION=A*GG*ABD
    RETURN
!
   30 SIGION=0.0
    RETURN
!
    END FUNCTION SIGION
!
!-----------------------------------------------------------------------
!
! Function INV finds the bin number closest to energy ETA on grid ENER.
! Bin INV or INV-1 will contain ETA.
!
    INTEGER FUNCTION INV (ETA, JY, ENER, NBINS)
    implicit none
!
    integer,intent(in) :: NBINS
    real,intent(in) :: ETA,ENER(NBINS)
    integer,intent(in) :: JY
!
    integer :: iv
!
    if (ETA .LT. 0.) then
      INV = -1
    else
      do IV=1,JY
        IF (ETA .LE. ENER(IV)) GOTO 40
      enddo
      IV = JY
    40 INV = IV
    endif
!
    RETURN

    END FUNCTION INV
!
!-----------------------------------------------------------------------
!
! Subroutine HEXC
!
! High Energy Cross Section Correction
! Calculates ratio of low energy (non-relativistic) to high energy
! (relativistic) ionization cross sections, based on N2.
! Extends to 1 GeV.
!
! Originally coded by Ann Windnagel, 11/98
! Re-written by Stan Solomon, 2/99
! Re-designed with table lookup, SCS, 4/99
! Updated comments, SCS, 4/02
! References:
!   Porter et al., J. Chem. Phys., 65, 154, 1976.
!   Rieke and Prepejchal, Phys. Rev. A, 6, 1507, 1990.
!   Saksena et al., Int. Jour. of Mass Spec. & Ion Proc., 171, L1, 1997.


    SUBROUTINE HEXC(ENER,SIGIX,RATIO,NBINS)
    use cglow,only: NEI,NMAJ
    implicit none

    integer,intent(in) :: NBINS
    real,intent(in) :: ENER(NBINS),SIGIX(NEI,NMAJ,NBINS)
    real,intent(out) :: RATIO(NBINS)

    real ::  TOTX(NBINS), TOTNEW(NBINS), EGR(13), SGR(13)
    DATA EGR/1.E4,    2.E4,    5.E4,    1.E5,    2.E5, &
          3.E5,    5.E5,    1.E6,    2.E6,    5.E6, &
          1.E7,    1.E8,    1.E9/
    DATA SGR/1.20E-17,  7.03E-18,  3.37E-18,  1.96E-18,  1.26E-18, &
          1.05E-18,  9.50E-19,  9.00E-19,  9.00E-19,  9.40E-19, &
          1.00E-18,  1.26E-18,  1.59E-18/
    integer :: k,i,kg
    real,external :: TERPOO

! Calculate total low-energy cross section for N2:

    do K = 1,NBINS
      TOTX(K) = 0.
      do I = 1,NEI
        TOTX(K) = TOTX(K) + SIGIX(I,3,K)
      enddo
    enddo


    ! Calculate high-energy cross section for N2, using tabulated values:

    do K=1,NBINS
      if (ENER(K) .GE. EGR(1)) then
        do KG=1,12
            if (ENER(K) .GE. EGR(KG) .AND. ENER(K) .LT. EGR(KG+1)) then
                TOTNEW(K)=TERPOO(ENER(K),EGR(KG),EGR(KG+1),SGR(KG),SGR(KG+1))
            endif
        enddo
      ELSE
        TOTNEW(K)=TOTX(K)
      ENDIF
    enddo


    ! Calculate ratio (=1 < 10 keV):

    do K = 1,NBINS
      if (ENER(K) .GE. EGR(1)) then
        RATIO(K) = TOTX(K)/TOTNEW(K) 
!       if (RATIO(K) .GT. 1.) RATIO(K) = 1.
      else
        RATIO(K) = 1.
      endif
    enddo

    RETURN

    END SUBROUTINE HEXC

!-----------------------------------------------------------------------

    REAL FUNCTION TERPOO(X,X1,X2,Y1,Y2)

    implicit none

    real,intent(in) :: x,x1,x2,y1,y2

    TERPOO = EXP ( ALOG(Y1) + ALOG(X/X1)*ALOG(Y2/Y1)/ALOG(X2/X1) )

    RETURN

    END FUNCTION TERPOO
