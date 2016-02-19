subroutine init_cond_isen(rho, vex, vey, pre)
   use comvar
   implicit none

   real    :: rho(-1:nx+2, -1:ny+2)
   real    :: vex(-1:nx+2, -1:ny+2)
   real    :: vey(-1:nx+2, -1:ny+2)
   real    :: pre(-1:nx+2, -1:ny+2)

   integer :: i, j
   real    :: x, y, r2, Temp, beta, c1, c2, mach_inf
   real    :: x0, y0, theta

   final_time = 100.0

   xmin =-5.0
   xmax = 5.0
   ymin =-5.0
   ymax = 5.0

   dx = (xmax - xmin)/nx
   dy = (ymax - ymin)/ny

   x0 = 0.0
   y0 = 0.0
   mach_inf = 0.5
   theta = 0.0
   beta = 5.0

   theta= theta*M_PI/180.0
   c1= (gamma-1.0)*beta**2/(8.0*gamma*M_PI**2)
   c2= beta/(2.0*M_PI)

   do i=1,nx
      do j=1,ny
         x = xmin + (i-1)*dx + 0.5*dx
         y = ymin + (j-1)*dy + 0.5*dy

         r2   = (x-x0)**2 + (y-y0)**2

         Temp     = 1.0 - c1*exp(1.0-r2)
         rho(i,j) = Temp**(1.0/(gamma-1.0))
         vex(i,j) = mach_inf * cos(theta) - c2*(y-y0)*exp(0.5*(1.0-r2))
         vey(i,j) = mach_inf * sin(theta) + c2*(x-x0)*exp(0.5*(1.0-r2))
         pre(i,j) = rho(i,j)**gamma
      enddo
   enddo

   write(*,*)'Set isentropic vortex as initial condition'

end subroutine init_cond_isen
