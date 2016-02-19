program main

   use comvar

   implicit none

   real, dimension(:), allocatable :: rho, vex, vey, pre, co0, co1, res
   real, dimension(:), allocatable :: omg

   nx = 100
   ny = 100

   ! Default value, modified in initial condition function
   final_time = 10.0

   itmax = 50000
   itsave= 50

   ! Material properties
   gas_const = 1.0
   mu        = 0.0

   ! periodicity conditions
   xperiod = yes
   yperiod = yes

   ! options: iroe, irusanov
   fluxtype = irusanov

   ! limiter: ford, muscl3, mmod
   limtype = muscl3

   ! file id for saving solution
   fileid_sol = 0
   fileid_omg = 0

   cfl = 0.5

   nrk    = 3

   if(nrk.eq.2)then
      ark(1) = 0.0
      ark(2) = 0.5
   else if(nrk.eq.3)then
      ark(1) = 0.0
      ark(2) = 3.0/4.0
      ark(3) = 1.0/3.0
   endif

   allocate( rho( (nx+4)*(ny+4) ) )
   allocate( vex( (nx+4)*(ny+4) ) )
   allocate( vey( (nx+4)*(ny+4) ) )
   allocate( pre( (nx+4)*(ny+4) ) )
   allocate( omg( (nx+1)*(ny+1) ) )
   allocate( co0( nvar*(nx+4)*(ny+4) ) )
   allocate( co1( nvar*(nx+4)*(ny+4) ) )
   allocate( res(nvar*(nx+2)*(ny+2)) )

   call solveFVM(rho, vex, vey, pre, omg, co0, co1, res)

end program main
