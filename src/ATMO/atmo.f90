subroutine msis00_init(data_dir, sw)
   ! Initialize MSIS model data from file
   implicit none
   integer, save :: snoem_initialized = 0
   character(len=*), intent(in) :: data_dir ! directory containing data files
   integer, intent(in) :: sw(25) ! MSIS switch array

   if (snoem_initialized .ne. 0) then
      call snoem_init(data_dir)
      snoem_initialized = 1
   end if

   call tselec(sw)
end subroutine msis00_init

subroutine msis00_eval(iyd,sec,alt,glat,glong,f107a,f107,ap,&
   d,t,exot,nalt)
   implicit none
   ! MSIS Legacy subroutine arguments
   integer, intent(in)         :: iyd
   real(4), intent(in)         :: sec
   real(4), intent(in)         :: alt(nalt)
   real(4), intent(in)         :: glat
   real(4), intent(in)         :: glong
   real(4), intent(in)         :: f107a
   real(4), intent(in)         :: f107
   real(4), intent(in)         :: ap(7)
   integer, intent(in)         :: nalt
   real(4), intent(inout)      :: d(10,nalt), t(nalt)
   real(4), intent(out)        :: exot
   integer                     :: i, j
   real(4)                     :: tmpd(9), tmpt(2), stl, mlon, mlat
   ! Calculate local solar time
   stl = sec/3600. + glong/15.
   if (stl < 0.) stl = stl + 24.
   if (stl >= 24.) stl = stl - 24.
   ! Call MSIS to get neutral densities and temperature
   do i=1,nalt
      call GTD7(iyd, sec, alt(i), glat, glong, stl, f107a, f107, &
         ap, 48, tmpd, tmpt)
      t(i)=tmpt(2)
      do j=1,9
         d(j,i)=tmpd(j)
      end do
   end do
   exot = tmpt(1)
   ! Calculate magnetic latitude
   call GEOMAG(0,glong,glat,mlon,mlat)
   ! Call SNOEMINT to obtain NO profile from the Nitric Oxide Empirical Model (NOEM)
   call snoemint(iyd, mlat, f107, ap, nalt, alt, t, d(10,:))
end subroutine msis00_eval

subroutine iri90_eval(jf,jmag,glat,glong,mmdd,sec,f107a,z,jmax, &
   iri90_dir,outf,oarr)
   ! jf,jmag,alat,alon,iyyy,mmdd,dhour,zkm,nzkm,outf,oarr
   implicit none

   logical,intent(in) :: jf(12)
   integer,intent(in) :: jmag, jmax, mmdd
   real,intent(in) :: glat, glong, f107a, z(jmax), sec
   character(len=*), intent(in) :: iri90_dir
   real,intent(inout) :: oarr(30)
   real,intent(out) :: outf(11,jmax)

   integer :: iday, j
   real :: rz12, stl

   stl = sec/3600. + glong/15.
   if (stl < 0.) stl = stl + 24.
   if (stl >= 24.) stl = stl - 24.
   rz12 = -f107a
   iday = -mod(mmdd,1000)

   call iri90(jf, jmag, glat, glong, rz12, iday, stl, z, jmax, iri90_dir, outf, oarr)

   do j=1,jmax
      outf(1,j) = outf(1,j) / 1.E6 ! electron density, m-3 to cm-3
      if (outf(1,j) < 100.) outf(1,j) = 100. ! minimum electron density
      outf(5:11,j) = outf(5:11,j) / 100. ! ion species densities
   enddo

end subroutine iri90_eval
