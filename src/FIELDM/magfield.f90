subroutine dipangle_pogo86(lat,lon,alt,dip)
    implicit none
    real, intent(in)  :: lat, lon, alt
    real, intent(out) :: dip
    real              :: x,y,z,f,dec,smodip
    call FIELDM(lat,lon,alt,x,y,z,f,dip,dec,smodip)
end subroutine

subroutine diparray_pogo86(lat,lon,alt,dip,nalt)
    implicit none
    integer, intent(in) :: nalt
    real, intent(in)    :: lat, lon, alt(nalt)
    real, intent(inout) :: dip(nalt)
    real                :: x,y,z,f,dec,smodip
    integer             :: i
    do i=1,nalt
        call FIELDM(lat,lon,alt(i),x,y,z,f,dip(i),dec,smodip)
    enddo
end subroutine

subroutine bmagarray_pogo86(lat,lon,alt,bmag,nalt)
    implicit none
    integer, intent(in) :: nalt
    real, intent(in)    :: lat, lon, alt(nalt)
    real, intent(inout) :: bmag(nalt)
    real                :: x,y,z,dip,dec,smodip
    integer             :: i
    do i=1,nalt
        call FIELDM(lat,lon,alt(i),x,y,z,bmag(i),dip,dec,smodip)
        bmag(i) = bmag(i)*1.0e-4
    enddo
end subroutine