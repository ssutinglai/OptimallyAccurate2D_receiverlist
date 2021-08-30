program SincInterpolation !( nx,nz,rho,lam,mu,dx,dz,dt )
  implicit none
  integer nx,nz
  integer ix,iz
  integer jx,jz
  integer mx,mz,nx,nz ! for trial functions
  integer mxmax,mxmin,mzmax,mzmin ! the max and min values of mx,mz,nx,nz 
  double precision xm,zm
  integer ndis,ngrid,npTF
  logical, parameter :: sincfunction = .FALSE.
  double precision x,xx,z,zz
  !double precision rho(nx+1,nz+1),lam(nx+1,nz+1),mu(nx+1,nz+1)
  double precision, allocatable :: lam(:,:),mu(:,:),rho(:,:)
  double precision dx,dz,dt,smalldx,smalldz
  double precision dx2,dz2,dxdz,dt2
  double precision, allocatable :: phix(:),phiz(:) ! non-zero values for phix, phiz
  double precision, allocatable :: phixderiv(:),phizderiv(:) ! and theirs derivatives
  double precision, allocatable :: H1(:,:),H2(:,:)
  double precision, parameter :: pi = 3.141592653589793238462643383
 

  dt = 1.d0
  dx = 1.d0
  dz = 1.d0

  
  !npTF defines points in scheme (3, 5, 7)
  npTF = 3
  ngrid = (npTF-1)/2
  ndis = 100


  smalldx = dx/dble(ndis)
  smalldz = dz/dble(ndis)


  mxmin = -npTF+1
  mxmax = npTF-1

  mzmin = -npTF+1
  mzmax = npTF-1

  allocate(phix(-ngrid*ndis:ngrid*ndis))
  allocate(phiz(-ngrid*ndis:ngrid*ndis))
  allocate(phixderiv(-ngrid*ndis:ngrid*ndis))
  allocate(phizderiv(-ngrid*ndis:ngrid*ndis))
  allocate(lam(-ngrid*ndis:ngrid*ndis,-ngrid*ndis:ngrid*ndis))
  allocate(H1(mxmin:mxmax,mzmin:mzmax))
  allocate(H2(mxmin:mxmax,mzmin:mzmax))
  

  lam =1.d0

  dt2 = dt * dt
  dx2 = dx * dx
  dz2 = dz * dz
  dxdz = dx * dz
  
  !xm = 0.d0
  !zm = 0.d0


  ! I put nx, nz to be the centre (0,0)

  nx = 0
  nz = 0

  ! then mx, mz have mxmin:mxmax and mzmin:mzmax


  ! Initialising vectors

  phix = 0.d0
  phiz = 0.d0
  
  phixderiv = 0.d0
  phizderiv = 0.d0


  H1 = 0.d0
  H2 = 0.d0

  !trialfunction decides on sinc(true) or linear(false) interpolation
  !sincfunction = .true.



  ! Trial function calculations for sinc/spline
  
  if (sincfunction) then
     
     do ix=-ngrid*ndis, ngrid*ndis

        xx =dble(ix/ndis)*dx
        phix(ix) = sin(pi*xx)/(pi*xx)
        phixderiv(ix) = pi*cos(pi*xx)/(pi*xx) - sin(pi*xx)/(pi*xx*xx)
        phiz(ix) = phix(ix)
        phizderiv(ix) =phixderiv(ix)
        
          

        open(unit=8,file="phix.dat",form="formatted"&
             ,status="replace",action="write")
        
        write(8,*)'phix', phix(ix,iz)
        
        close(8)
        
        
     enddo
     
     
     
  else
     

     ! B-splines (for the moment only for 3 points so it's not correct for 5-, 7- point schemes)

     do ix=-ngrid*ndis,0
        
        x=dble(ix/ndis)*dx
        phix(ix) = (x+dx)/dx
        phiz(ix) = phix(ix)
        phixderiv(ix) = 1.d0/dx
        phizderiv(ix) = phixderiv(ix)
        
     enddo
     
     do ix=0,ngrid*ndis
      
        x=dble(ix/ndis)*dx
        phix(ix) = (-x+dx)/dx
        phiz(ix) = phix(ix)
        phixderiv(ix) = -1.d0/dx
        phizderiv(ix) = phixderiv(ix)
        
     enddo
  endif

  

  nx = 0
  nz = 0


  do mx = mxmin,mxmax
     do mz = mzmin,mzmax
        
        H1(mx,mz) = 0.d0
        
        do ix = -ngrid*ndis,ngrid*ndis
           if(((ix-ndis).ge.-ngrid*ndis).and.((ix-ndis).le.ngrid*ndis))
           H1(mx,mz) = H1(mx,mz) + phix(ix-ndis)*phix(ix)*
        enddo
     enddo
  enddo


     
     






end program SincInterpolation
  
  
