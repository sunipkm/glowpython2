module noem
   character(1024) :: data_dir

   real :: snoem_zin(16)          ! altitude grid
   real :: snoem_mlatin(33)       ! magnetic latitude grid
   real :: snoem_no_mean(33,16)   ! mean nitric oxide distribution
   real :: snoem_eofs(33,16,3)    ! empirical orthogonal functions

contains
   subroutine snoem_init(direct)
      ! Initialize SNOEM model data from file
      implicit none
      character(len=*), intent(in) :: direct ! directory containing data file
      character(len=1024) :: filepath
      integer :: j,k,n
      data_dir = trim(direct)
      filepath = trim(data_dir)//'/'//'snoem_eof.dat'
      open(unit=1,file=filepath,status='old',action='read')
      read(1,*) (snoem_zin(k),k=1,16)
      read(1,*) (snoem_mlatin(j),j=1,33)
      read(1,*) ((snoem_no_mean(j,k),j=1,33),k=1,16)
      read(1,*) (((snoem_eofs(j,k,n),j=1,33),k=1,16),n=1,3)
      close(unit=1)
   end subroutine snoem_init

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
end module
