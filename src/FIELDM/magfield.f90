subroutine dipangle_pogo68(year,lat,lon,alt,dip)
   ! Calculate dip angle using POGO68 model
   implicit none
   real, intent(in)  :: year ! Decimal year: Year + (Day of Year)/365 or 366
   real, intent(in)  :: lat  ! latitude
   real, intent(in)  :: lon  ! longitude
   real, intent(in)  :: alt  ! altitude
   real, intent(out) :: dip  ! dip angle (output)
   real              :: x,y,z,f,dec,smodip
   call FIELDM(lat,lon,alt,x,y,z,f,dip,dec,smodip)
end subroutine

subroutine diparray_pogo68(year,lat,lon,alt,dip,nalt)
   ! Calculate dip angle array using POGO68 model
   implicit none
   integer, intent(in) :: nalt ! number of altitude points
   real, intent(in)    :: year ! Decimal year: Year + (Day of Year)/365 or 366
   real, intent(in)    :: lat  ! latitude
   real, intent(in)    :: lon  ! longitude
   real, intent(in)    :: alt(nalt) ! altitude array
   real, intent(inout) :: dip(nalt) ! dip angle array (output)
   real                :: x,y,z,f,dec,smodip
   integer             :: i
   do i=1,nalt
      call FIELDM(lat,lon,alt(i),x,y,z,f,dip(i),dec,smodip)
   enddo
end subroutine

subroutine bmagarray_pogo68(year,lat,lon,alt,bmag,nalt)
   ! Calculate magnetic field strength array using POGO68 model
   implicit none
   integer, intent(in) :: nalt       ! number of altitude points
   real, intent(in)    :: year       ! Decimal year: Year + (Day of Year)/365 or 366
   real, intent(in)    :: lat        ! latitude
   real, intent(in)    :: lon        ! longitude
   real, intent(in)    :: alt(nalt)  ! altitude array
   real, intent(inout) :: bmag(nalt) ! magnetic field strength array, Tesla (output)
   real                :: x,y,z,dip,dec,smodip
   integer             :: i
   do i=1,nalt
      call FIELDM(lat,lon,alt(i),x,y,z,bmag(i),dip,dec,smodip)
      bmag(i) = bmag(i)*1.0e-4
   enddo
end subroutine

subroutine mlatlon_pogo68(year,lat,lon,alt,mlat,mlon)
   ! Calculate magnetic latitude and longitude using POGO68 model
   implicit none
   real, intent(in)  :: year ! Decimal year: Year + (Day of Year)/365 or 366
   real, intent(in)  :: lat  ! geographic latitude
   real, intent(in)  :: lon  ! geographic longitude
   real, intent(in)  :: alt  ! altitude
   real, intent(out) :: mlat ! magnetic latitude (output)
   real, intent(out) :: mlon ! magnetic longitude (output)
   call GEOMAG(0, lon, lat, mlon, mlat)
end subroutine
