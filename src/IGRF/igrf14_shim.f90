module glowigrf
    real :: byear = 0.0
    character(1024) :: data_dir

contains
    subroutine igrf_init(direct)
        implicit none
        COMMON /CONST/UMR,PI &
        /IGRF1/ERA,AQUAD,BQUAD,DIMO &
        /findRLAT/FLON,RYEAR &
        /iounit/konsol,mess
        character(len=*), intent(in) :: direct ! directory containing data file
        real :: umr, pi, era, aquad, bquad, dimo, flon, ryear
        real :: eexc, erequ, erpol
        integer :: konsol, mess
        ! Initialize constants
        pi=ATAN(1.0)*4.
        UMR=pi/180.
        ERA=6371.2
        EREQU=6378.16
        ERPOL=6356.775
        AQUAD=EREQU*EREQU
        BQUAD=ERPOL*ERPOL
        EEXC=0.01675
        dimo=0.311653
        mess=0 ! no messages
        konsol=6 ! standard output unit
        ! Set data directory
        data_dir = trim(direct)
    end subroutine igrf_init

    subroutine dipangle_igrf(year,lat,lon,alt,dip)
        ! Calculate dip angle using IGRF model
        implicit none
        real, intent(in)  :: year ! Decimal year: Year + (Day of Year)/365 or 366
        real, intent(in)  :: lat ! latitude
        real, intent(in)  :: lon ! longitude
        real, intent(in)  :: alt ! altitude
        real, intent(out) :: dip ! dip angle (output)
        real              :: dec,dipl,ymodip
        if (byear .ne. year) then
            print *, 'Loading IGRF data for year ', year
            call FELDCOF(year,data_dir)
            print *, 'IGRF data loaded for year ', year
            byear = year
        end if
        call igrf_dip(lat,lon,year,alt,dec,dip,dipl,ymodip)
    end subroutine

    subroutine diparray_igrf(year,lat,lon,alt,dip,nalt)
        ! Calculate dip angle array using IGRF model
        implicit none
        integer, intent(in) :: nalt ! number of altitude points
        real, intent(in)    :: year ! Decimal year: Year + (Day of Year)/365 or 366
        real, intent(in)    :: lat  ! latitude
        real, intent(in)    :: lon  ! longitude
        real, intent(in)    :: alt(nalt) ! altitude array
        real, intent(inout) :: dip(nalt) ! dip angle array (output)
        real                :: dec,dipl,ymodip
        integer             :: i
        if (byear .ne. year) then
            call FELDCOF(year,data_dir)
            byear = year
        end if
        do i=1,nalt
            call igrf_dip(lat,lon,year,alt(i),dec,dip(i),dipl,ymodip)
        enddo
    end subroutine

    subroutine bmagarray_igrf(year,lat,lon,alt,bmag,nalt)
        ! Calculate magnetic field strength array using IGRF model
        implicit none
        integer, intent(in) :: nalt       ! number of altitude points
        real, intent(in)    :: year       ! Decimal year: Year + (Day of Year)/365 or 366
        real, intent(in)    :: lat        ! latitude
        real, intent(in)    :: lon        ! longitude
        real, intent(in)    :: alt(nalt)  ! altitude array
        real, intent(inout) :: bmag(nalt) ! magnetic field strength array (output)
        real                :: x,dip
        integer             :: i,icode
        if (byear .ne. year) then
            call FELDCOF(year,data_dir)
            byear = year
        end if
        do i=1,nalt
            print *, 'Calculating Bmag for alt=', alt(i)
            call igrf_sub(lat,lon,year,alt(i),x,icode,dip,bmag(i))
            bmag(i) = bmag(i)*1.0e-4
        enddo
    end subroutine

    subroutine mlatlon_igrf(year,lat,lon,alt,mlat,mlon)
        ! Calculate magnetic latitude and longitude using IGRF model
        implicit none
        real, intent(in)  :: year ! Decimal year: Year + (Day of Year)/365 or 366
        real, intent(in)  :: lat  ! geographic latitude
        real, intent(in)  :: lon  ! geographic longitude
        real, intent(in)  :: alt  ! altitude
        real, intent(out) :: mlat ! magnetic latitude (output)
        real, intent(out) :: mlon ! magnetic longitude (output)
        real              :: dat(11,4),pla(4),plo(4)
        dat(:,:) = 0.0
        dat(1,1) = lat
        dat(1,2) = lon
        call GEOCGM01(1,int(year),alt,dat,pla,plo)
        mlat = dat(1,3)
        mlon = dat(1,4)
    end subroutine
end module