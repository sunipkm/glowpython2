program fieldmtest
    implicit none
    real :: dip, alt(5), dip_array(5), bmag_array(5), mlat, mlon
    integer :: i
    call dipangle_pogo68(2020.12485, 0.0, 0.0, 0.0, dip)
    print *, 'Dip angle at equator, sea level in 2020:', dip
    alt = [0.0, 100.0, 200.0, 300.0, 400.0]
    call diparray_pogo68(2020.12485, 0.0, 0.0, alt, dip_array, 5)
    print *, 'Dip angles at various altitudes:'
    do i = 1, 5
        print *, 'Altitude:', alt(i), 'Dip angle:', dip_array(i)
    end do
    call bmagarray_pogo68(2020.12485, 0.0, 0.0, alt, bmag_array, 5)
    print *, 'Magnetic field strengths at various altitudes:'
    do i = 1, 5
        print *, 'Altitude:', alt(i), 'Bmag:', bmag_array(i)
    end do
    call mlatlon_pogo68(2020.12485, 45.0, 45.0, 1.0, mlat, mlon)
    print *, 'Magnetic latitude and longitude for (45N, 45E):', mlat, mlon
end program fieldmtest