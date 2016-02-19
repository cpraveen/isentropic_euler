subroutine cons2prim(con, rho, vex, vey, pre)
   use comvar
   implicit none

   real :: con(nvar, -1:nx+2, -1:ny+2)
   real :: rho(-1:nx+2, -1:ny+2)
   real :: vex(-1:nx+2, -1:ny+2)
   real :: vey(-1:nx+2, -1:ny+2)
   real :: pre(-1:nx+2, -1:ny+2)

   integer :: i, j

   do i=-1,nx+2
      do j=-1,ny+2
         rho(i,j) = con(1,i,j)
         vex(i,j) = con(2,i,j)/con(1,i,j)
         vey(i,j) = con(3,i,j)/con(1,i,j)
         pre(i,j) = rho(i,j)**gamma
      enddo
   enddo

end subroutine cons2prim
