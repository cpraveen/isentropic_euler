module comvar
   implicit none

   integer :: nx, ny
   real    :: xmin, xmax, ymin, ymax, dx, dy, dt
   integer :: itmax
   integer :: itsave
   real    :: cfl
   real    :: final_time

   real    :: mu, gas_const

   integer :: nrk
   real    :: ark(3)

   real :: gamma = 1.4
   real :: M_PI = 4.0*atan(1.0)

   integer :: fileid_sol, fileid_omg

   integer :: fluxtype
   integer :: iroe=1, irusanov=2

   integer :: limtype
   integer :: ford=0, muscl3=1, mmod=2

   integer :: no=0, yes=1

   integer :: xperiod, yperiod
   integer,parameter :: nvar=3

end module comvar
