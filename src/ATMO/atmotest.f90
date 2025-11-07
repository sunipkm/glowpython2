program atmotest
    implicit none
    logical :: jf(12)
    real :: d(10,3), t(3), exot
    real :: outf(11,3), oarr(30)
    integer :: j, sw(25)
    sw(:) = 1
    call msis00_init('data', sw)
    call msis00_eval(2020124, 43200.0, [100.0, 200.0, 300.0], 45.0, 45.0, 150.0, 110.0, 4.0, d, t, exot, 3)
    print *, 'Densities (cm^-3):'
    print *, '    O:', (d(1,j), j=1,3)
    print *, '    N2:', (d(2,j), j=1,3)
    print *, '    O2:', (d(3,j), j=1,3)
    print *, 'Temperatures (K):', (t(j), j=1,3)
    print *, 'Exospheric Temperature (K):', exot
    jf(:) = .true.
    call iri90_eval(jf,0,45.0,45.0,124,43200.0,150.0,[100.0, 200.0, 300.0],'data',outf,oarr,3)
    print *, 'Electron Densities (cm^-3):', (outf(1,j), j=1,3)
    print *, 'Ion Temperatures (K):', (outf(3,j), j=1,3)
    print *, 'Electron Temperatures (K):', (outf(4,j), j=1,3)
    print *, 'nmF2 (cm^-3):', oarr(1) / 1.e6, ' hmF2 (km):', oarr(2)
end program atmotest