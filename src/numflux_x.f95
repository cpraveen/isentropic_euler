subroutine numflux_x(conjm1, conj, conjp1, conjp2, flux)
   use comvar
   implicit none

   real :: conjm1(nvar), conj(nvar), conjp1(nvar), conjp2(nvar)
   real :: flux(nvar)

   real :: conl(nvar), conr(nvar)

   ! reconstructed states
   call reconstruct(conjm1, conj, conjp1, conjp2, conl, conr)

   if(fluxtype == iroe)then
      call roe_flux(1.0, 0.0, conl, conr, flux)
   else if(fluxtype == irusanov)then
      call rusanov_flux(1.0, 0.0, conl, conr, flux)
   else
      write(*,*)'Uknown flux type fluxtype =', fluxtype
      stop
   endif

end subroutine numflux_x
