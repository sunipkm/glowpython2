program igrftest
    ! COMMON /CONST/UMR,PI &
    ! /IGRF1/ERA,AQUAD,BQUAD,DIMO &
    ! /findRLAT/FLON,RYEAR &
    ! /iounit/konsol,mess
    ! pi=ATAN(1.0)*4.
    ! UMR=pi/180.
    ! ERA=6371.2
    ! EREQU=6378.16
    ! ERPOL=6356.775
    ! AQUAD=EREQU*EREQU
    ! BQUAD=ERPOL*ERPOL
    ! EEXC=0.01675
    ! dimo=0.311653
    ! mess=0
    ! konsol=6
    ! call FELDCOF(2019.21643,'/home/sunip/Codes/python/glowpython2/src/IGRF/data')
    ! ! CALL GEODIP(IYEAR,LATI,LONGI,MLAT,MLONG,JMAG)
    ! call igrf_dip(0.0,0.0,2019.21643,300.0,dec,dip,magbr,modip)
    ! call igrf_sub(0.0, 0.0, 2019.21643, 110.0, x, icode, dip, bmag)
    ! print *, 'Magnetic field strength at lat=45, lon=45, alt=60 km in 2018.0 is', bmag, 'Gauss'
    use glowigrf, only: igrf_init, dipangle_igrf, diparray_igrf, bmagarray_igrf, mlatlon_igrf
    implicit none
    integer, parameter :: nalt = 5
    real :: alt(nalt), dip_array(nalt), bmag_array(nalt), dip
    real :: mlat, mlon
    integer :: i
    call igrf_init('/home/sunip/Codes/python/glowpython2/src/IGRF/data')
    call dipangle_igrf(2019.21643, 45.0, 45.0, 300.0, dip)
    print *, 'Dip angle at lat=45, lon=45, alt=300 km in 2019.21643 is', dip, 'degrees'
    alt = [100.0, 200.0, 300.0, 400.0, 500.0]
    call diparray_igrf(2019.21643, 45.0, 45.0, alt, dip_array, nalt)
    print *, 'Dip angle array at lat=45, lon=45 in 2019.21643:'
    do i = 1, nalt
        print *, '  Altitude =', alt(i), 'km, Dip angle =', dip_array(i), 'degrees'
    end do
    call bmagarray_igrf(2019.21643, 45.0, 45.0, alt, bmag_array, nalt)
    print *, 'Magnetic field strength array at lat=45, lon=45 in 2019.21643:'
    do i = 1, nalt
        print *, '  Altitude =', alt(i), 'km, Bmag =', bmag_array(i), 'Tesla'
    end do
    call mlatlon_igrf(2019.21643, 45.0, 45.0, 300.0, mlat, mlon)
    print *, 'Magnetic latitude and longitude at lat=45, lon=45, alt=300 km in 2019.21643 are:', mlat, 'degrees,', mlon, 'degrees'
end program igrftest