module atmo
   character(1024) :: data_dir

   real :: snoem_zin(16)          ! altitude grid
   real :: snoem_mlatin(33)       ! magnetic latitude grid
   real :: snoem_no_mean(33,16)   ! mean nitric oxide distribution
   real :: snoem_eofs(33,16,3)    ! empirical orthogonal functions

contains

   subroutine atmo_init(direct)
      ! Initialize SNOEM model data from file
      implicit none
      character(len=*), intent(in) :: direct ! directory containing data file
      character(len=1024) :: filepath
      integer :: j,k,n
      data_dir = trim(direct)
      filepath = trim(data_dir)//'snoem_eof.dat'
      open(unit=1,file=filepath,status='old',action='read')
      read(1,*) (snoem_zin(k),k=1,16)
      read(1,*) (snoem_mlatin(j),j=1,33)
      read(1,*) ((snoem_no_mean(j,k),j=1,33),k=1,16)
      read(1,*) (((snoem_eofs(j,k,n),j=1,33),k=1,16),n=1,3)
      close(unit=1)
   end subroutine atmo_init

   subroutine snoem(doy, kp, f107, z, mlat, nozm)
      ! Subroutine SNOEM calculates nitric oxide zonal mean altitude profile
      ! as function of magnetic latitude for specified day of year, Kp, and F10.7.
      !
      ! The NOEM empirical model is based on data from the SNOE ultraviolet
      ! spectrometer during 1998-2000, using empirical orthogonal function analysis.
      ! Altitude range is from 100 to 150 km.
      ! Reference:
      ! Marsh et al., JGR, 109, A07301, doi:10.1029/2003JA010199, 2004.
      !
      ! Adapted by Stan Solomon, 5/2014, from IDL and F90 code supplied by Dan Marsh.
      ! Minor revisions to make compatible with gfortran, SCS, 4/2017

      implicit none

      integer,intent(in) :: doy
      real,intent(in) :: kp, f107
      real,intent(out) :: z(16), mlat(33), nozm(33,16)

      real :: theta0                 ! day number in degrees
      real :: dec                    ! solar declination angle
      real :: m1, m2, m3             ! coefficients for first 3 eofs
      real, parameter :: pi=3.1415926536
      integer :: j, k

      !... calculate coefficients (m1 to m3) for eofs based on geophysical parameters

      !... eof1 - kp

      m1 =  kp * 0.689254 - 1.53366

      !... eof2 - declination

      theta0 = 2.*pi * float(doy - 1) / 365.

      dec = 0.006918 &
         - 0.399912 * cos(theta0)   + 0.070257 * sin(theta0) &
         - 0.006758 * cos(2*theta0) + 0.000907 * sin(2*theta0) &
         - 0.002697 * cos(3*theta0) + 0.001480 * sin(3*theta0)

      dec = dec * 180./pi

      m2 = -0.31978 + dec*0.097309 + dec**2*0.00048979 - dec**3*0.00010360

      !... eof3 - f107

      m3 =  alog10(f107) * 6.35777 - 13.8163

      !... zonal mean distrib. is sum of mean and eofs

      do k=1,16
         do j=1,33
            nozm(j,k) = snoem_no_mean(j,k)-m1*snoem_eofs(j,k,1)+m2*snoem_eofs(j,k,2)-m3*snoem_eofs(j,k,3)
         end do
      end do

      z(:) = snoem_zin(:)
      mlat(:) = snoem_mlatin(:)

      return

   end subroutine snoem

   subroutine snoemint(idate,mlat,f107,ap,jmax,z,ztn,zno)
      ! Subroutine SNOEMINT gets NO estimate from the NOEM emperical model and
      ! INTerpolates it onto an altitude grid.  Extrapolation is done above 150
      ! km assuming a scale height approximation, and below 100 km
      ! assuming a constant number density profile.
      !
      ! Stan Solomon, 12/2014
      ! Refactored to f90, scs, 6/2016
      ! Fixed bug in Kp (estimated from Ap) so that xkp >= 0, scs, 1/2017
      !
      ! Input:
      !   IDATE  Date in yyddd or yyyyddd format
      !   MLAT   Magnetic latitude in degrees
      !   F107   10.7 cm radio flux index
      !   AP     Ap index
      !   JMAX   Number of points in altitude grid
      !   Z      Altitude grid in km
      !   ZTN    Temperature at Z in K
      ! Output:
      !   ZNO    Nitric oxide density at Z in cm-3

      implicit none

      integer,intent(in) :: idate, jmax
      real,intent(in) :: mlat, f107, ap, z(jmax), ztn(jmax)
      real,intent(out) :: zno(jmax)

      real,parameter :: pi=3.1415926536

      integer :: iday, klat1, klat2, kz1, kz2, j
      real :: zg(16), xmlatno(33), zmno(33,16), zmnoi(16)
      real :: xkp, rat, h

      ! Find magnetic latitude:

      ! call geomag(0,glong,glat,xmlong,xmlat)

      ! Get zonal mean no profiles:

      iday=mod(idate, 1000)
      xkp=1.75*alog(0.4*ap)
      if (xkp < 0.0) xkp=0.0
      call snoem(iday,xkp,f107,zg,xmlatno,zmno)

      ! Interpolate altitude profile at magnetic latitude:

      klat1=ifix(mlat+80.)/5+1
      klat2=klat1+1
      if (klat1 < 1) klat1=1
      if (klat1 > 33) klat1=33
      if (klat2 < 1) klat1=1
      if (klat2 > 33) klat2=33
      rat=mlat/5.0-ifix(mlat)/5

      do j=1,16
         zmnoi(j) = alog(zmno(klat1,j)*(1.-rat)+zmno(klat2,j)*rat)
      end do

      ! Interpolate onto altitude grid:
      ! Use constant value below 100 km and scale height assumption above 150 km:

      h=0.03*ztn(jmax)
      do j=1,jmax
         if (z(j) <= 100.) zno(j)=exp(zmnoi(16))
         if (z(j) > 100. .and. z(j) <= 150.) then
            kz2=ifix((150.-z(j))*.3)+1
            kz1=kz2+1
            zno(j)=exp(zmnoi(kz1) + (zmnoi(kz2)-zmnoi(kz1)) * (z(j)-zg(kz1)) / (zg(kz2)-zg(kz1)))
         endif
         if (z(j) > 150.) zno(j)=exp(zmnoi(1)+(150.-z(j))/h)
      end do

      return

   end

   subroutine mzgrid (&
      jmax,nex,idate,ut,glat,glong,mlat,&
      stl,f107a,f107,f107p,ap,z, &
      jfin, oarr, &
      zo,zo2,zn2,zns,znd,zno,ztn,zun,zvn,ze,zti,zte,zxden)
      ! Subroutine MZGRID gets fields from eMpirical models on default GLOW altitude grid.
      !
      ! This software is part of the GLOW model.  Use is governed by the Open Source
      ! Academic Research License Agreement contained in the file glowlicense.txt.
      ! For more information see the file glow.txt.
      !
      ! Stan Solomon, 12/15, 1/16
      ! Extracted from glowdriver.f90 into separate file mzgrid.f90, SCS, 12/16
      !
      ! Neutral densities from NRLMSISE-00 model (a.k.a. MSIS2K).
      ! Nitric oxide densities from NOEM model (via snoem.f and snoemint.f).
      ! Electron densities from IRI-90.
      ! Inputs:
      !         jmax   Number of altitude levels used by GLOW (should = 102 for MSIS/IRI/NOEM runs)
      !         nex    Number of ionized/excited species (for array zxden)
      !         idate  Date in yyyyddd or yyddd format
      !         ut     Universal time, seconds)
      !         glat   Latitude, degrees
      !         glong  Longitude, degrees
      !         mlat   Magnetic latitude, degrees
      !         stl    Local solar time, hours
      !         f107a  F10.7 index 81-day centered average
      !         f107   F10.7 index of day
      !         f107p  F107. index of previous day
      !         ap     Ap index daily value
      !         z      Geographic altitude grid (km)
      !         iri90_dir  Directory containing IRI input files (set in namelist inputs)
      !
      ! Outputs:
      !         zo     O number density, cm-3
      !         zo2    O2    "
      !         zn2    N2    "
      !         zns    N(4S) "
      !         n2d    N(2D) "      (set to zero since this is calculated by GLOW)
      !         no     NO    "
      !         ztn    Tn, K
      !         zun    Zonal wind velocity, cm/s, currently = 0 (not used by GLOW)
      !         zvn    Meridional wind velocity, cm/s currently = 0 (not used by GLOW)
      !         zti    Ti, K (Ion temperature is at least neutral temperature)
      !         zte    Te, K
      !         zxden  Array of ionized/excited species density, cm-3, must be dimensioned (nex,jmax)
      !                O+(2P), O+(2D), O+(4S)[x], N+, N2+, O2+[x], NO+[x], N2(A), N(2P), N(2D), O(1S), O(1D)
      implicit none

      integer,intent(in) :: jmax,nex,idate,jfin(12)
      real,intent(in) :: ut,glat,glong,stl,f107a,f107,f107p,ap,z(jmax),mlat
      real,intent(inout) :: oarr(30)
      real,intent(out) :: zo(jmax),zo2(jmax),zn2(jmax),zns(jmax),znd(jmax), &
         zno(jmax),ztn(jmax),zti(jmax),zte(jmax),zun(jmax),zvn(jmax),ze(jmax),zxden(nex,jmax)

      integer :: j,jmag,iday,mmdd
      real :: rz12, d(8), t(2), sw(25)
      logical :: jf(12)
      character(1024) :: iri90_dir
      real,allocatable :: outf(:,:)              ! iri output (11,jmax)
      data sw/25*1./

      iri90_dir = trim(data_dir)//'iri90'

      allocate(outf(11,jmax))

      ! Call MSIS-2K to get neutral densities and temperature:
      !
      ! call alt_grid(jmax, 60., 0.5, 4., z)
      call tselec(sw)

      do j=1,jmax ! levels
         call gtd7(idate,ut,z(j),glat,glong,stl,f107a,f107p,ap,48,d,t) ! mass number, 0 is temperature only, 48 is all, 17 is anomalous O only
         zo(j) = d(2)
         zn2(j) = d(3)
         zo2(j) = d(4)
         zns(j) = d(8)
         ztn(j) = t(2)
         znd(j)  = 0. ! calculated by glow
      enddo
      !
      ! Call SNOEMINT to obtain NO profile from the Nitric Oxide Empirical Model (NOEM):
      !
      call snoemint(idate,mlat,f107,ap,jmax,z,ztn,zno)
      !
      ! Call International Reference Ionosphere-1990 subroutine to get
      ! electron density, electron temperature, and ion temperature:
      ! The directory iri90_dir is the location of the ccirnn.asc and ursinn.asc files.
      !
      do j=1,12
         if (jfin(j) == 0) then
            jf(j) = .false.
         else
            jf(j) = .true.
         endif
      enddo

      jmag = 0
      rz12 = -f107a
      iday = mod(idate,1000)
      mmdd = -iday
      outf = 0.

      call iri90(jf,jmag,glat,glong,rz12,mmdd,stl,z,jmax, iri90_dir,outf,oarr)

      do j=1,jmax
         ze(j) = outf(1,j) / 1.E6
         if (ze(j) < 100.) ze(j) = 100.
         zti(j) = outf(3,j)
         if (zti(j) < ztn(j)) zti(j) = ztn(j)
         zte(j) = outf(4,j)
         if (zte(j) < ztn(j)) zte(j) = ztn(j)
         zxden(3,j) = ze(j) * outf(5,j)/100.
         zxden(6,j) = ze(j) * outf(8,j)/100.
         zxden(7,j) = ze(j) * outf(9,j)/100.
      enddo
!
! Until implementation of an empirical neutral wind model, winds are set to zero:
!
      zun(:) = 0.
      zvn(:) = 0.

   end subroutine mzgrid


end module
