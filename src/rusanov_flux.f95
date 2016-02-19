subroutine rusanov_flux(sx, sy, conl, conr, flux)
   use comvar
   implicit none

   real :: sx, sy, conl(nvar), conr(nvar), flux(nvar)

   real :: rhol, vexl, veyl, prel, sonl, vnl, laml
   real :: rhor, vexr, veyr, prer, sonr, vnr, lamr
   real :: lam

   ! left state
   rhol = conl(1)
   vexl = conl(2)/conl(1)
   veyl = conl(3)/conl(1)
   prel = rhol**gamma
   sonl = sqrt(gamma*prel/rhol)
   vnl  = vexl*sx + veyl*sy
   laml = max( abs(vnl - sonl), abs(vnl), abs(vnl + sonl) )

   ! right state
   rhor = conr(1)
   vexr = conr(2)/conr(1)
   veyr = conr(3)/conr(1)
   prer = rhor**gamma
   sonr = sqrt(gamma*prer/rhor)
   vnr  = vexr*sx + veyr*sy
   lamr = max( abs(vnr - sonr), abs(vnr), abs(vnr + sonr) )

   lam  = max(laml, lamr)

   flux(1) = 0.5*(rhol*vnl + rhor*vnr - lam*(conr(1)-conl(1)))
   flux(2) = 0.5*(prel*sx + rhol*vexl*vnl + prer*sx + rhor*vexr*vnr - &
             lam*(conr(2)-conl(2)))
   flux(3) = 0.5*(prel*sy + rhol*veyl*vnl + prer*sy + rhor*veyr*vnr - &
             lam*(conr(3)-conl(3)))

end subroutine rusanov_flux
