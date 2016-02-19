subroutine solveFVM(rho, vex, vey, pre, omg, co0, co1, res)

   use comvar

   implicit none

   real :: rho(-1:nx+2, -1:ny+2)
   real :: vex(-1:nx+2, -1:ny+2)
   real :: vey(-1:nx+2, -1:ny+2)
   real :: pre(-1:nx+2, -1:ny+2)
   real :: omg( 1:nx+1,  1:ny+1)
   real :: co0(nvar, -1:nx+2, -1:ny+2)
   real :: co1(nvar, -1:nx+2, -1:ny+2)
   real :: res(nvar,0:nx+1,0:ny+1)

   integer :: it, i, j, rks
   real    :: lambda
   real    :: xflux(nvar), yflux(nvar)
   real    :: time
   real    :: resid(nvar), resid1(nvar)
   real    :: ke, entropy
   logical :: tostop


   ! set initial condition
   call init_cond(rho, vex, vey, pre)
   call prim2cons(rho, vex, vey, co1)
   call periodic(co1)
   call cons2prim(co1, rho, vex, vey, pre)
   call saveprim(0.0, rho, vex, vey, pre)
   call vorticity(vex, vey, omg)
   call savevort(0.0, omg)
   call timestep(rho, vex, vey, pre)

   time   = 0.0
   it     = 0

   do while(time < final_time .and. it < itmax)

      ! Exactly match final time
      tostop = .false.
      if(time + dt > final_time)then
         dt = final_time - time
         tostop = .true.
      endif

      lambda = dt/dx/dy

      co0(:,:,:) = co1(:,:,:)

      do rks=1,nrk

         res = 0.0

         ! x fluxes
         do i=0,nx
            do j=1,ny
               call numflux_x(co1(:,i-1,j), co1(:,i,j), co1(:,i+1,j), &
                              co1(:,i+2,j), xflux)
               res(:,i,j)   = res(:,i,j)   + dy*xflux(:)
               res(:,i+1,j) = res(:,i+1,j) - dy*xflux(:)
            enddo
         enddo

         ! y fluxes
         do j=0,ny
            do i=1,nx
               call numflux_y(co1(:,i,j-1), co1(:,i,j), co1(:,i,j+1), &
                              co1(:,i,j+2), yflux)
               res(:,i,j)   = res(:,i,j)   + dx*yflux(:)
               res(:,i,j+1) = res(:,i,j+1) - dx*yflux(:)
            enddo
         enddo

         ! update conserved variables
         resid = 0.0
         do i=1,nx
            do j=1,ny
               co1(:,i,j) = ark(rks)*co0(:,i,j) + &
                            (1.0-ark(rks))*(co1(:,i,j) - lambda*res(:,i,j))
               resid = resid + res(:,i,j)**2
            enddo
         enddo
         resid = sqrt(resid)

         call periodic(co1)

      enddo ! Rk stage loop

      it = it + 1
      if(it==1)then
         resid1 = 1.0
         if(resid(1) > 0.0) resid1(1) = resid(1)
         if(resid(2) > 0.0) resid1(2) = resid(2)
         if(resid(3) > 0.0) resid1(3) = resid(3)
      endif

      call cons2prim(co1,rho,vex,vey,pre)
      call global_quantities(rho,vex,vey,pre,ke,entropy)

      time = time + dt
      write(*,'(I6,F10.2,7E12.4)')it,time,resid(:)/resid1(:),ke,entropy

      if(mod(it,itsave)==0 .or. it==itmax .or. tostop)then
         call cons2prim(co1,rho,vex,vey,pre)
         call saveprim(time, rho, vex, vey, pre)
         call vorticity(vex, vey, omg)
         call savevort(time, omg)
      endif

   enddo ! time iteration loop


end subroutine solveFVM
