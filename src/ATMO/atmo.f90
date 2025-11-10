subroutine msis00_init(data_dir, sw)
   use noem, only: snoem_init
   ! Initialize MSIS model data from file
   implicit none
   integer, save :: snoem_initialized = 0
   character(len=*), intent(in) :: data_dir ! directory containing data files
   real, intent(in), optional :: sw(25) ! MSIS switch array
   real :: usesw(25)
   if (present(sw)) then
      usesw(:) = sw(:)
   else
      usesw(:) = 1
   end if
   if (snoem_initialized .eq. 0) then
      call snoem_init(data_dir)
      snoem_initialized = 1
   end if

   call tselec(usesw)
end subroutine msis00_init

subroutine msis00_eval(iyd,sec,alt,glat,glong,f107a,f107,ap,&
   d,t,exot,nalt)
   use noem, only: snoemint
   implicit none
   ! MSIS Legacy subroutine arguments
   integer, intent(in)      :: iyd
   real, intent(in)         :: sec
   real, intent(in)         :: alt(nalt)
   real, intent(in)         :: glat
   real, intent(in)         :: glong
   real, intent(in)         :: f107a
   real, intent(in)         :: f107
   real, intent(in)         :: ap(7)
   integer, intent(in)      :: nalt
   real, intent(inout)      :: d(10,nalt), t(nalt)
   real, intent(out)        :: exot
   integer                  :: i, j
   real                     :: tmpd(9), tmpt(2), stl, mlon, mlat, no(nalt)
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
   call snoemint(iyd, mlat, f107, ap(1), nalt, alt, t, no)
   d(10,:) = no(:)
end subroutine msis00_eval

subroutine iri90_eval(jf,jmag,glat,glong,mmdd,sec,f107a,z, &
   data_dir,outf,oarr,jmax)
   ! jf,jmag,alat,alon,iyyy,mmdd,dhour,zkm,nzkm,outf,oarr
   implicit none

   logical,intent(in) :: jf(12)
   integer,intent(in) :: jmag, jmax, mmdd
   real,intent(in) :: glat, glong, f107a, z(jmax), sec
   character(len=*), intent(in) :: data_dir
   real,intent(inout) :: oarr(30)
   real,intent(inout) :: outf(11,jmax)
   character(len=1024) :: iri90_dir
   
   integer :: iday, i, j
   real :: rz12, stl
   
   iri90_dir = trim(data_dir)//'/iri90/'
   stl = sec/3600. + glong/15.
   if (stl < 0.) stl = stl + 24.
   if (stl >= 24.) stl = stl - 24.
   rz12 = -f107a
   iday = -mod(mmdd,1000)

   call iri90(jf, jmag, glat, glong, rz12, iday, stl, z, jmax, iri90_dir, outf, oarr)

   do j=1,jmax
      if (outf(11,j).lt.0.0) outf(11,j) = 0.0
      if (outf(4,j).lt.outf(3,j)) outf(4,j) = outf(3,j) ! Te >= Ti
      do i=5,11
         outf(i,j) = outf(i,j) * outf(1,j) / 100.0 ! percentage to cm-3
      enddo
   enddo
end subroutine iri90_eval
