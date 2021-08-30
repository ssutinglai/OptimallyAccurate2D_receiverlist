
subroutine calf2(nx,nz,it,t,ist,isx,isz,dt,dx,dz,rho,f0,t0,fx,fz )

  implicit none
  real*8 pi
  parameter ( pi=3.1415926535897932d0 )
  integer nx,nz,it,ist,isx,isz
  real*8 t,dt,dx,dz,rho,f0,t0,fx(nx+1,nz+1),fz(nx+1,nz+1)
  real*8 b,a,factor

   factor=1.d3
    if ( it.le.ist ) then

     a = pi*pi*f0*f0*(t-t0)*(t-t0)

     ! Ricker source time function (second derivative of a Gaussian)

     fx(isx,isz) = -factor * (1.d0 - 2.d0*a)*exp(-a); !(04.04.2018)

     fx(isx,isz) = fx(isx,isz) * dt * dt / rho !%! Closed for matching with Specfem2D (08.03.2018)
     fz(isx,isz) = fx(isx,isz)
     if ( (it.eq.0).or.(it.eq.ist) ) then
        fx(isx,isz) = fx(isx,isz) / 2.d0
        fz(isx,isz) = fz(isx,isz) / 2.d0
     endif

  else
     fx(isx,isz) = 0.d0
     fz(isx,isz) = 0.d0

  endif


  !NF for point source
  fx(isx,isz)=0.d0



  return
end subroutine calf2



subroutine calf( nx,nz,it,t,ist,isx,isz,dt,dx,dz,rho,tp,ts,fx,fz )
  implicit none
  real*8 pi
  parameter ( pi=3.1415926535897932d0 )
  integer nx,nz,it,ist,isx,isz
  real*8 t,dt,dx,dz,rho,tp,ts,fx(nx+1,nz+1),fz(nx+1,nz+1)
  real*8 b
  
  if ( it.le.ist ) then
     b = pi * ( t - ts ) / tp
     fx(isx,isz) &
          = dsqrt(pi) / 2.d0 * (b*b-0.5d0) * dexp(-b*b) &
          / ( dx * dz )
     fx(isx,isz) = fx(isx,isz) * dt * dt / rho
     fz(isx,isz) = fx(isx,isz)
     if ( (it.eq.0).or.(it.eq.ist) ) then
        fx(isx,isz) = fx(isx,isz) / 2.d0
        fz(isx,isz) = fz(isx,isz) / 2.d0
     endif
  else
     fx(isx,isz) = 0.d0
     fz(isx,isz) = 0.d0
  endif

  ! NF for point source
  !! NOT change fz and fx into different direction 28/02/2018
  fx(isx,isz)=0.d0
  fz(isx,isz)=fz(isx,isz)*1.d4   !Change for matching with specfem2d
  !fz(isx,isz)=fz(isx,isz)

!%! Plot source wavelet!


  return
end subroutine calf
