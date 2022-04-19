!------------------------------------------------------------------------------


module skyvault

integer :: ndir    ! Number of directions used for sky discretisation
real, allocatable :: hmoy(:) ! Elevation angle of direction jdir (jdir=1,ndir)
real, allocatable :: azmoy(:) ! Azimuth angle of direction jdir (jdir=1,ndir)
real, allocatable :: omega(:) ! Solid angle associated to direction jdir (jdir=1,ndir)
real, allocatable :: pc(:)  ! Relative contribution of direction jdir (jdir=1,ndir) to incident diffuse radiation

contains
 subroutine sv_read(fname)

 use constant_values    ! for number pi

!  integer :: lectureciel  ! option for input file format of skyvault discretisation data

  real :: sumomega   ! sum of solid angle values (should be 2 pi after reading all direction data)
  real :: sumpc    ! sum of solid angle contribution to diffuse incident radiation (should be 1 after reading all direction data)



  character*28 fname
  !write(*,*) 'nom skyvault',fname

!  Lecture des fichiers d'entree  : Reading the input files

  !write(*,*) 'Creating the sky vault ...'
  !write(*,*)

!  select case (lectureciel)

!  Case (1)
!  Lecture des directions (h et az) des ndir directions, ainsi que de l'angle solide et la contribution au diffus incident
  !fname='skyvault.'//spec
  open (1,file=fname)
  read(1,*) ndir    ! Lecture du nombre de directions ndir

  call sv_destroy   ! deallocate allocatable arrays if already allocated

  allocate(hmoy(ndir))  ! Allocation des tableaux
  allocate(azmoy(ndir))
  allocate(omega(ndir))
  allocate(pc(ndir))

  sumomega = 0.
  sumpc=0.
  do jdir=1,ndir
   read(1,*) hmoy(jdir), azmoy(jdir), omega(jdir), pc(jdir)  ! angles in degrees
   hmoy(jdir)=hmoy(jdir)*pi/180.  ! conversion to radians
   azmoy(jdir)=azmoy(jdir)*pi/180.
   sumomega=sumomega+omega(jdir)
   sumpc=sumpc+pc(jdir)
  end do

!  end select
  close (1)

  !write(*,*) '  Total diffuse irradiance of sky vault: ',sumpc,'. SHOULD BE 1 '
  !write(*,*) '  Total solid angle of sky vault: ',sumomega,'. SHOULD BE ',2.*pi

 end subroutine sv_read

   subroutine sv_destroy

  if (allocated(hmoy))   deallocate(hmoy)
  if (allocated(azmoy))  deallocate(azmoy)
  if (allocated(omega))  deallocate(omega)
  if (allocated(pc))   deallocate(pc)

   end subroutine sv_destroy


end module skyvault

!------------------------------------------------------------------------------
