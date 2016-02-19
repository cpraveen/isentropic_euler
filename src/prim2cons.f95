subroutine prim2cons(rho, vex, vey, con)
   use comvar
   implicit none

   real :: rho(-1:nx+2, -1:ny+2)
   real :: vex(-1:nx+2, -1:ny+2)
   real :: vey(-1:nx+2, -1:ny+2)
   real :: con(nvar, -1:nx+2, -1:ny+2)

   integer :: i, j

   do i=-1,nx+2
      do j=-1,ny+2
         con(1,i,j) = rho(i,j)
         con(2,i,j) = rho(i,j)*vex(i,j)
         con(3,i,j) = rho(i,j)*vey(i,j)
      enddo
   enddo

end subroutine prim2cons
