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
  double precision dx,dz,dt
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

  mxmin = -npTF+1
  mxmax = npTF-1

  mzmin = -npTF+1
  mzmax = npTF-1

  allocate(phix(-ngrid*ndis:ngrid*ndis))
  allocate(phiz(-ngrid*ndis:ngrid*ndis))
  allocate(phixderiv(-ngrid*ndis:ngrid*ndis))
  allocate(phizderiv(-ngrid*ndis:ngrid*ndis))
  allocate(lam(-ngrid*ndis:ngrid*ndis,-ngrid*ndis:ngrid*ndis))
  allocate(H1(1,1))
  allocate(H2(1,1))
  

  lam =1.d0

  dt2 = dt * dt
  dx2 = dx * dx
  dz2 = dz * dz
  dxdz = dx * dz
  
  !xm = 0.d0
  !zm = 0.d0


  ! I put mx, mz to be the centre (0,0)

  mx = 0
  mz = 0

  ! then nx, nz have mxmin:mxmax and mzmin:mzmax


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


  

  do jx=-ngrid, ngrid
        do jz=-ngrid, ngrid
           
           H1(m,n)=(phix(jx,jz)*phix(jx,jz))+ &
                (phiz(jx,jz)*phiz(jx,jz))
           
           
           H2(m,n)=lam(jx,jz)*(phixderiv(jx,jz)+phizderiv(jx,jz))&
               *(phixderiv(jx,jz)+phizderiv(jx,jz))

           print *,'jx,jz',jx,jz,'H1', H1(m,n), 'H2', H2(m,n)
           
        end do
     end do
     






end program SincInterpolation
  
  
