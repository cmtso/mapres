!________________________________________________________
subroutine mapit(sigma,sat,porosity,m,n,mfname,pfname,nv)
  implicit none
  character*40, intent(in) :: mfname,pfname
  integer, intent(in) ::nv
  real*8, intent(in) :: m,n
  real*8, dimension(nv),intent(in) :: sigma,sat,porosity
  
!python directives
!f2py intent(in) nv
!f2py intent(in) sigma
!f2py intent(in) sat
!f2py intent(in) porosity
!f2py intent(in) m
!f2py intent(in) n
!f2py intent(in) mfname
!f2py intent(in) pfname

  integer*8 :: i
  integer*8 :: cpos
  integer :: pfnx,pfny,pfnz,nelem
  integer, dimension(:), allocatable :: rw,v
  real, dimension(:), allocatable :: w,sige4d
  
  !read in the map
  open(10,file=trim(mfname),status='old',form='unformatted',action='read')
  read(10) pfnx
  read(10) pfny
  read(10) pfnz
  read(10) nelem
  read(10) cpos
  allocate(rw(cpos),v(cpos),w(cpos))
  read(10) rw
  read(10) v
  read(10) w
  close(10)

  allocate(sige4d(nelem))
  sige4d = 0
  do i=1,cpos
     sige4d(rw(i))=sige4d(rw(i)) + (porosity(v(i))**m)*sigma(v(i))*sat(v(i))*n
  end do
 
  write(*,*) "mapit: writing sigma file: "//trim(pfname)//'.sig'
  open(10,file=trim(pfname)//'.sig',status='replace',action='write')
  write(10,*) nelem, "1"
  do i=1,nelem
     write(10,*) sige4d(i)
  end do
  close(10)
end subroutine mapit


!________________________________________________________
subroutine map_waxsmit(fc,sat,por,temp,bsig_fnam,pet_fnam,mp_fnam,tm,nv)
  implicit none
  character*80, intent(in) :: bsig_fnam,pet_fnam,mp_fnam
  integer, intent(in) ::nv
  real, intent(in) :: tm
  real, dimension(nv), intent(in) :: fc,sat,por,temp
  
  
!python directives
!f2py intent(hide), depend(fc) :: nv=shape(fc,0)
!f2py intent(in) fc
!f2py intent(in) sat
!f2py intent(in) por
!f2py intent(in) temp
!f2py intent(in) bsig_fname
!f2py intent(in) pet_fname
!f2py intent(in) mp_fname
!f2py intent(in) tm


  integer*8 :: i
  integer*8 :: cpos
  integer :: pfnx,pfny,pfnz,nelem,netest
  integer, dimension(:), allocatable :: rw,v
  real, dimension(:), allocatable :: w,sig_e4d,bsig
  real, dimension(:,:), allocatable :: petro
  real :: a,B,Qv,c,m,t,Tt,ifc,isat,ipor
  character*20 :: str
  
!read in the map matrix
  open(10,file=trim(mp_fnam),status='old',form='unformatted',action='read')
  read(10) pfnx
  read(10) pfny
  read(10) pfnz
  read(10) nelem
  read(10) cpos
  allocate(rw(cpos),v(cpos),w(cpos))
  read(10) rw
  read(10) v
  read(10) w
  close(10)

!read the uninterpolated conductivity values
write(*,*) "Reading uninterpolated conductivity: ",trim(bsig_fnam)
open(10,file=trim(bsig_fnam),status='old',action='read')
read(10,*) netest
if(nelem .ne. netest) then
  write(*,*) "The number of elements in the map matrix (",nelem,")"
  write(*,*) "does not match the number of elements (",netest,")" 
  write(*,*) "in the baseline sigma file ",trim(bsig_fnam)
  write(*,*) "aborting ..."
  close(10)
  return
end if
allocate(bsig(nelem))
do i=1,nelem
   read(10,*) bsig(i)
end do
close(10)

!read the Waxman-Smit petrophysical parmaters
write(*,*) "Reading Waxman_Smit petro parameters: ",trim(pet_fnam)
open(10,file=trim(pet_fnam),status='old',action='read')
read(10,*) netest
if(nelem .ne. netest) then
  write(*,*) "The number of elements in the map matrix (",nelem,")"
  write(*,*) "does not match the number of elements (",netest,")" 
  write(*,*) "in the petro file ",trim(pet_fnam)
  write(*,*) "aborting ..."
  close(10)
  return
end if
allocate(petro(nelem,6))
do i=1,nelem
   read(10,*) petro(i,:)
end do
close(10)

  !allocate and init the bulk conductivity vector
  allocate(sig_e4d(nelem))
 sig_e4d = 0

  !interpolate the conductivities via Waxman-Smit eqn
  write(*,*) 'Interpolating to e4d mesh'
  do i=1,cpos
     a = petro(rw(i),1)
     B = petro(rw(i),2)
     Qv = petro(rw(i),3)
     c = petro(rw(i),4)
     m = petro(rw(i),5)
     t = petro(rw(i),6)
     ifc = fc(v(i))
     ipor = por(v(i))
     isat = sat(v(i))
     Tt = temp(v(i))
     
     if(isat>0) then
     	sig_e4d(rw(i)) = sig_e4d(rw(i)) + w(i)*((ifc/a+B*Qv/isat)*(ipor**c)*(isat**m))
     end if
  end do
  
  !set the uniterpolated conductvities
  do i=1,nelem
     if(sig_e4d(i) == 0) then
     	sig_e4d(i) = bsig(i)
     end if
  end do
  
  !write the conductivity file for this time step
  write(str,*) tm
  write(*,*) 'Writing conductivity for time ',tm,' to file: '//trim(str)//'.sig'
  write(*,*)
  write(*,*)
  
  open(10,file=adjustl(trim(adjustl(str)))//'.sig',status='replace',action='write')
  write(10,*) nelem,1
  do i=1,nelem
     write(10,*) sig_e4d(i)
  end do
  close(10)
  
end subroutine map_waxsmit
