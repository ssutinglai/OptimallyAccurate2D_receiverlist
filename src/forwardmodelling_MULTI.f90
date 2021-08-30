subroutine forwardmodelling
  use parameters
  use paramFWI
  implicit none

  
  ! record size
  recl_size=(nx+1-rmargin(1)-lmargin(1))*(nz+1-rmargin(2)-lmargin(2))*kind(0e0)
  recl_size_syn=(maxnt+1)*(nReceiver+1)*kind(0e0)


  
  ! Smoothed version of CONV/OPT operators


  ! NF is not using cales for circle free surface (KATAYAMA)

  if(0.eq.1) then
  call cales( nx,nz,rho,lam,mu,dt,dx,dz, &
       e1, e2, e3, e4, e5, e6, e7, e8, &
       e13,e14,e15,e16,e17,e18,e19,e20, &
       f1, f2, f3, f4, f5, f6, f7, f8, &
       f13,f14,f15,f16,f17,f18,f19,f20 ) 
   endif

  ! NF uses free surface circle


  call cales_circle_free( nx,nz,rho,lam,mu,dt,dx,dz, &
       e1, e2, e3, e4, e5, e6, e7, e8, &
       e13,e14,e15,e16,e17,e18,e19,e20, &
       f1, f2, f3, f4, f5, f6, f7, f8, &
       f13,f14,f15,f16,f17,f18,f19,f20 )

  ! discontinuities
  
  ee12 = 0.d0
  ee34 = 0.d0
  ee56 = 0.d0
  ee65 = 0.d0
  ee78 = 0.d0
  ee87 = 0.d0

  ff12 = 0.d0
  ff34 = 0.d0
  ff56 = 0.d0
  ff65 = 0.d0
  ff78 = 0.d0
  ff87 = 0.d0
  
  if(nDiscon.ne.0) then
 
     ! changing dscr by putting lmargin(1) and (2)
     tmpvaluex=dble(lmargin(1))*dx
     tmpvaluez=dble(lmargin(2))*dz
     do ix=1,nDiscon
        do iz=1,lengthDiscon
           dscr(1,iz,ix)=dscr(1,iz,ix)+tmpvaluex
           dscr(2,iz,ix)=dscr(2,iz,ix)+tmpvaluez
        enddo
     enddo
          
     call cales_discon( nx,nz,rho,lam,mu,dt,dx,dz,e1, e2, e3, e4, e5, e6, e7, e8,&
     e13,e14,e15,e16,e17,e18,e19,e20, &
     f1, f2, f3, f4, f5, f6, f7, f8, &
     f13,f14,f15,f16,f17,f18,f19,f20, & 
     ! hereafter are new variables for cales_discon
     ee12,ee34,ee56,ee65,ee78,ee87, &
     ff12,ff34,ff56,ff65,ff78,ff87, &
     markers,nDiscon,lengthDiscon,dscr)

  endif
  
  if(lengthFreeSurface.ne.0) then
          
     call cales_free( maxnx,nx,nz,rho,lam,mu,dt,dx,dz,e1, e2, e3, e4, e5, e6, e7, e8,&
     e13,e14,e15,e16,e17,e18,e19,e20, &
     f1, f2, f3, f4, f5, f6, f7, f8, &
     f13,f14,f15,f16,f17,f18,f19,f20, & 
     ! hereafter are new variables for cales_discon
     ee12,ee34,ee56,ee65,ee78,ee87, &
     ff12,ff34,ff56,ff65,ff78,ff87, &
     zerodisplacement,lengthFreeSurface,free)

     
  endif



  if(lengthFreeSurface.ne.0) then
     ! changing free by putting lmargin(1) and (2)
     tmpvaluex=dble(lmargin(1))*dx
     tmpvaluez=dble(lmargin(2))*dz
     
     do ix=1,lengthFreeSurface
           free(1,ix)=free(1,ix)+tmpvaluex
           free(2,ix)=free(2,ix)+tmpvaluez           
     enddo
  endif

  ! for Cerjan absorbing boundary


  weightBC=1.d0
     
  call compNRBCpre(weightBC(1:nx+1,1:nz+1),CerjanRate,lmargin,rmargin,nx+1,nz+1)
  




  do iSource = 1, nSource

     isx=iisx(iSource)+lmargin(1)
     isz=iisz(iSource)+lmargin(2)

     ist=nt/4

!     write(*,*)isx,isz

     ALPHA_MAX_PML = 2.d0*PI*(f0/2.d0) ! from Festa and Vilotte
     
     ! Initializing the data

     ux=0.d0
     uz=0.d0
     ux1=0.d0
     uz1=0.d0
     ux2=0.d0
     uz2=0.d0
     work=0.d0
     

     
 
     ! R. Courant et K. O. Friedrichs et H. Lewy (1928)
     cp=maxval(vp)
     Courant_number = cp * dt * sqrt(1.d0/dx**2 + 1.d0/dz**2)
     print *, 'Courant number is', Courant_number
     
     

     fx=0.d0
     fz=0.d0
     


     
     t=0.d0
     time(0)=t
     do it=0,nt
       
        call calf2( nx,nz,it,t,ist,isx,isz,dt,dx,dz,rho(isx,isz),f0,t0,fx,fz )
        t=t+dt
        !write(13,*) t, fx(isx,isz),fz(isx,isz)
        
     enddo
     !print *, maxnz,it,t,ist,isx,isz,dt,dx,dz,rho(isx,isz),f0,t0
     !stop

     
     t = 0.d0

     do ir = 1,nReceiver
        synx(0,ir)=ux(nrx(ir),nrz(ir))
        synz(0,ir)=uz(nrx(ir),nrz(ir))
!        synp(0,ir)=ux(nrx(ir),nrz(ir))
!        syns(0,ir)=uz(nrx(ir),nrz(ir))             !%!05/04/2018 Commented
     enddo



     do it=0,nt
      
        call calf2( nx,nz,it,t,ist,isx,isz,dt,dx,dz,rho(isx,isz),f0,t0,fx,fz )
        ! evaluating the next step
        
        !if(nDiscon.eq.0) then
!           call calstep( maxnz,nx,nz, &
!                e1, e2, e3, e4, e5, e6, e7, e8, &
!                e13,e14,e15,e16,e17,e18,e19,e20, &
!                f1, f2, f3, f4, f5, f6, f7, f8, &
!                f13,f14,f15,f16,f17,f18,f19,f20, &
!                ux,uz,ux1,ux2,uz1,uz2,isx,isz,fx,fz, &
!                work(1,1), work(1,5), work(1,9),work(1,13), &
!                work(1,17),work(1,18),work(1,20),work(1,21), &
!                work(1,23),work(1,24),work(1,28),work(1,29), optimise)
       ! endif
        !else

        ! NF uses cales_circle_free but normally you should use calstep_discon or calstep
! !%! Ssu-Ting close it for test
           call calstep_circle_free( nx,nz, &
                e1, e2, e3, e4, e5, e6, e7, e8, &
                e13,e14,e15,e16,e17,e18,e19,e20, &
                f1, f2, f3, f4, f5, f6, f7, f8, &
                f13,f14,f15,f16,f17,f18,f19,f20, &
                ux,uz,ux1,ux2,uz1,uz2,isx,isz,fx,fz, &
                work(1,1), work(1,5), work(1,9),work(1,13), &
                work(1,17),work(1,18),work(1,20),work(1,21), &
                work(1,23),work(1,24),work(1,28),work(1,29), optimise, &
                ! Hereafter are new variables for cales_discon
                ee12,ee34,ee56,ee65,ee78,ee87, &
                ff12,ff34,ff56,ff65,ff78,ff87)


        !endif


           if(lengthFreeSurface.ne.0) then
              do iz=1,nz+1
                 do ix=1,nx+1
                    if(zerodisplacement(ix,iz).eq.1) then
                       ux(ix,iz) = 0.d0
                       uz(ix,iz) = 0.d0
                    endif
                 enddo
              enddo
           endif


           ! increment of t
        t = t + dt
        time(it)=t
        !write(14,*) real(t),real(ux(nrx,nrz)),real(uz(nrx,nrz))
        do ir = 1,nReceiver                                  !%! Recording as seismic records
           synx(it,ir)=ux(nrx(ir),nrz(ir))
           synz(it,ir)=uz(nrx(ir),nrz(ir))
        enddo

     
        
        ! applying Cerjan boundary
        
        do iz=1,nz+1
           do ix=1,nx+1
              uz(ix,iz)=uz(ix,iz)*weightBC(ix,iz)
              ux1(ix,iz)=ux1(ix,iz)*weightBC(ix,iz)
              uz1(ix,iz)=uz1(ix,iz)*weightBC(ix,iz)
              ux2(ix,iz)=ux2(ix,iz)*weightBC(ix,iz)
              uz2(ix,iz)=uz2(ix,iz)*weightBC(ix,iz)
           enddo
        enddo
        
        
        
        ! calculating strains
        
        
        !if(writingStrain.and.(mod(it,IT_DISPLAY).eq.0)) then
           singleStrainDiagonal=0.e0
           singleStrainShear=0.e0
           tmpsingleStrain=0.e0
           call calStrainDiagonal(nx,nz,ux,uz,lmargin,rmargin,singleStrainDiagonal)
           call calStrainShear(nx,nz,ux,uz,lmargin,rmargin,singleStrainShear)

           do ir=1,nReceiver
              synp(it,ir)=singleStrainDiagonal(nrx(ir),nrz(ir))
              syns(it,ir)=singleStrainShear(nrx(ir),nrz(ir))
           enddo

!%! Close for nor print out strain (04.04.2018)
!           if(optimise) then
!              write(outfile,'("strainD",I5,".",I5,".",I5,".OPT_dat") ') it,isx-lmargin(1),isz-lmargin(2)
!           else
!              write(outfile,'("strainD",I5,".",I5,".",I5,".CON_dat") ') it,isx-lmargin(1),isz-lmargin(2)
!           endif
!           do j=1,24
!              if(outfile(j:j).eq.' ') outfile(j:j)='0'
!           enddo
!           
!           outfile = './strains/'//trim(modelname)//'/'//outfile
!           open(1,file=outfile,form='unformatted',access='direct',recl=recl_size)
!
!           tmpsingleStrain(1:nx+1-lmargin(1)-rmargin(1),1:nz+1-lmargin(2)-rmargin(2)) = &
!                singleStrainDiagonal(lmargin(1)+1:nx+1-rmargin(1),lmargin(2)+1:nz+1-rmargin(2))
!           write(1,rec=1)  tmpsingleStrain
!           close(1,status='keep')

           tmpsingleStrain=0.e0
!%! Close for nor print out strain (04.04.2018)        
!           if(optimise) then
!              write(outfile,'("strainS",I5,".",I5,".",I5,".OPT_dat") ') it,isx-lmargin(1),isz-lmargin(2)
!           else
!              write(outfile,'("strainS",I5,".",I5,".",I5,".CON_dat") ') it,isx-lmargin(1),isz-lmargin(2)
!           endif
!           do j=1,24
!              if(outfile(j:j).eq.' ') outfile(j:j)='0'
!           enddo
!           
!           outfile = './strains/'//trim(modelname)//'/'//outfile
!           open(1,file=outfile,form='unformatted',access='direct',recl=recl_size)
!
!           tmpsingleStrain(1:nx+1-lmargin(1)-rmargin(1),1:nz+1-lmargin(2)-rmargin(2)) = &
!                singleStrainShear(lmargin(1)+1:nx+1-rmargin(1),lmargin(2)+1:nz+1-rmargin(2))
!           write(1,rec=1)  tmpsingleStrain
!           close(1,status='keep')
!!


           
        !endif
!write(*,*)iisx,iisz

        
        !write(*,*) it, ' of ', nt
        if(mod(it,IT_DISPLAY) == 0)then
           !
           !head=0
           !head(58) = nx
           !head(59) = dz * 1E3
           !snapux=0.e0
           !snapuz=0.e0
           !snapux(1:nx,1:nz) = ux(1:nx,1:nz)
           !snapuz(1:nx,1:nz) = uz(1:nx,1:nz)
           !write(routine,'(a12,i5.5,a9)') './snapshots/',it,'snapUx.su'
           !open(21,file=routine,access='stream')
           !do j = 1,nx,1
           !   write(21) head,(real(snapux(k,j)),k=1,nz)
           !enddo
           !close(21)
           !write(routine,'(a12,i5.5,a9)') './snapshots/',it,'snapUz.su'
           !open(21,file=routine,access='stream')
           !do j = 1,nx,1
           !   write(21) head,(real(snapuz(k,j)),k=1,nz)
           !enddo
           !close(21)
           
           
           
           !call create_color_image(ux(1:nx+1,1:nz+1),nx+1,nz+1,it,isx,isz,ix_rec,iz_rec,1,0, &
           !     dummylog,dummylog,dummylog,dummylog,1)
           !call create_color_image(ux(1:nx+1,1:nz+1),nx+1,nz+1,it,isx,isz,ix_rec,iz_rec,1,&
           !    NPOINTS_PML,USE_PML_XMIN,USE_PML_XMAX,USE_PML_YMIN,USE_PML_YMAX,1)

           if(videoornot) then
              call create_color_image(ux(1:nx+1,1:nz+1),nx+1,nz+1,it,nSource, &   !%! Change for x-component plot
                   nrx(1:nReceiver),nrz(1:nReceiver),nReceiver,&
                   NPOINTS_PML,USE_PML_XMIN,USE_PML_XMAX,USE_PML_YMIN,USE_PML_YMAX,2,lmargin)
              
              
           endif

!           if(videoornot) then
!              call create_color_image(ux(1:nx+1,1:nz+1),nx+1,nz+1,it,isx,isz, &   !%! Change for x-component plot
!                   nrx(1:nReceiver),nrz(1:nReceiver),nReceiver, &
!                   NPOINTS_PML,USE_PML_XMIN,USE_PML_XMAX,USE_PML_YMIN,USE_PML_YMAX,2)
!
!
!           endif





           !if(optimise) then
           !   write(outfile,'("video",I5,".",I5,".",I5,".OPT_UX") ') it,isx-lmargin(1),isz-lmargin(2)
           !else
           !   write(outfile,'("video",I5,".",I5,".",I5,".CON_UX") ') it,isx-lmargin(1),isz-lmargin(2)
           !endif
           !do j=1,24
           !   if(outfile(j:j).eq.' ') outfile(j:j)='0'
           !enddo
           
           !outfile = './synthetics/'//trim(modelname)//'/'//outfile
           !open(1,file=outfile,form='unformatted',access='direct',recl=recl_size)
           !video(1:nx+1-lmargin(1)-rmargin(1),1:nz+1-lmargin(2)-rmargin(2))= &
           !     ux(lmargin(1)+1:nx+1-rmargin(1),lmargin(2)+1:nz+1-rmargin(2))
           !write(1,rec=1)  video(1:nx+1-lmargin(1)-rmargin(1),1:nz+1-lmargin(2)-rmargin(2))
           !close(1,status='keep')
           
           
           
           !if(optimise) then
           !   write(outfile,'("video",I5,".",I5,".",I5,".OPT_UX") ') it,isx-lmargin(1),isz-lmargin(2)
           !else
           !   write(outfile,'("video",I5,".",I5,".",I5,".CON_UX") ') it,isx-lmargin(1),isz-lmargin(2)
           !endif
           !do j=1,24
           !   if(outfile(j:j).eq.' ') outfile(j:j)='0'
           !enddo
           
           !outfile = './synthetics/'//trim(modelname)//'/'//outfile
           !open(1,file=outfile,form='unformatted',access='direct',recl=recl_size)
           !video(1:nx+1-lmargin(1)-rmargin(1),1:nz+1-lmargin(2)-rmargin(2))= &
           !     ux(lmargin(1)+1:nx+1-rmargin(1),lmargin(2)+1:nz+1-rmargin(2))
           !write(1,rec=1)  video(1:nx+1-lmargin(1)-rmargin(1),1:nz+1-lmargin(2)-rmargin(2))
           !close(1,status='keep')
           
        endif

        
        
        !call compNRBC2(ux(1:nx+1,1:nz+1),ux1(1:nx+1,1:nz+1),ux2(1:nx+1,1:nz+1), &
        !     uz(1:nx+1,1:nz+1),uz1(1:nx+1,1:nz+1),uz2(1:nx+1,1:nz+1), CerjanRate, lmargin, rmargin,nx+1,nz+1)
        
     enddo



     !write(18,*) singleStrainDiagonal(:,:)


     if(videoornot) then
        
        if(optimise) then
           write(outfile,'("video",".",I5,".",I5,".OPT.mp4") ') isx-lmargin(1),isz-lmargin(2)
        else
           write(outfile,'("video",".",I5,".",I5,".CON.mp4") ') isx-lmargin(1),isz-lmargin(2)
        endif
        do j=1,24
           if(outfile(j:j).eq.' ') outfile(j:j)='0'
        enddo
        
        outfile = './videos/'//trim(modelname)//'/'//trim(outfile)
        
      
        commandline="ffmpeg -framerate 5 -pattern_type glob -i 'snapshots/*.png' -c:v libx264 -pix_fmt yuv420p -vf 'scale=trunc(iw/"
        commandline=trim(commandline)//"2)*2:trunc(ih/2"//")*2' "
        commandline=trim(commandline)//" "//trim(outfile)
        !print *, commandline
        call system(commandline)

        
     
     endif
     
     

  
  
     if(iterationIndex.eq.0) then
        if(optimise) then
           write(outfile,'(I5,".",I5,".OPT_UX") ')  &
                isx-lmargin(1),isz-lmargin(2)
        else
           write(outfile,'(I5,".",I5,".CON_UX") ') &
                isx-lmargin(1),isz-lmargin(2)
        endif
        
        do j=1,12
           if(outfile(j:j).eq.' ') outfile(j:j)='0'
        enddo
        
        outfile = './synthetics/'//trim(modelname)//'/'//outfile
     else
        if(optimise) then
           write(outfile,'(I5,".",I5,".OPT_UX.it",I3.3) ')  &
                isx-lmargin(1),isz-lmargin(2),iterationIndex
        else
           write(outfile,'(I5,".",I5,".CON_UX,it",I3.3) ') &
                isx-lmargin(1),isz-lmargin(2),iterationIndex
        endif
        
        do j=1,12
           if(outfile(j:j).eq.' ') outfile(j:j)='0'
        enddo
        
        outfile = './synthetics/'//trim(modelname)//'/'//outfile
        
     endif
        
        


     open(1,file=outfile,form='unformatted',access='direct',recl=recl_size_syn)
     write(1,rec=1) synx(0:maxnt,1:nReceiver)
     close(1)

     
     if(iterationIndex.eq.0) then
        
        if(optimise) then
           write(outfile,'(I5,".",I5,".OPT_UZ") ') &
                isx-lmargin(1),isz-lmargin(2)

        else
           write(outfile,'(I5,".",I5,".CON_UZ") ') &
                isx-lmargin(1),isz-lmargin(2)
        endif
        
        do j=1,12
           if(outfile(j:j).eq.' ') outfile(j:j)='0'
        enddo
        
        outfile = './synthetics/'//trim(modelname)//'/'//outfile
        
     else

          
        if(optimise) then
           write(outfile,'(I5,".",I5,".OPT_UZ.it",I3.3) ') &
                isx-lmargin(1),isz-lmargin(2),iterationIndex
        else
           write(outfile,'(I5,".",I5,".CON_UZ.it",I3.3) ') &
                isx-lmargin(1),isz-lmargin(2),iterationIndex
        endif
        
        do j=1,12
           if(outfile(j:j).eq.' ') outfile(j:j)='0'
        enddo
        
        outfile = './synthetics/'//trim(modelname)//'/'//outfile
     endif
        

     !open(1, file=outfile,status='unknown',form='formatted')
     !do it=0,nt
     !   write (1,*) time(it), synz(it,ir)
     !enddo
     !close(1)
     
     open(1,file=outfile,form='unformatted',access='direct',recl=recl_size_syn)
     write(1,rec=1) synz(0:maxnt,1:nReceiver)
     close(1)
     !open(1,file=outfile,form='unformatted',access='direct',recl=kind(0e0)*(nt+1))
        !write(1,rec=1) synz(0:nt,ir)
        !close(1,status='keep')


    !%! For records in different receivers in the vertical line in the middle of the circle (01/03/2018 Ssu-Ting)

    if(iterationIndex.eq.0) then
      do ir=1,nReceiver
       if(optimise)then
          write(outfile,'(I5,".",I5,".OPT_UX") ') &
                nrx(ir)-lmargin(1),nrz(ir)-lmargin(2)
        else
          write(outfile,'(I5,".",I5,".CON_UX") ') &
                nrx(ir)-lmargin(1),nrz(ir)-lmargin(2)
       endif

       do j=1,12
          if(outfile(j:j).eq.' ') outfile(j:j)='0'
       enddo

       outfile = './synthetics/'//trim(modelname)//'/'//outfile

       open(1,file=outfile,form='unformatted',access='direct',recl=recl_size_syn)
       write(1,rec=1) synx(0:maxnt,ir)
       close(1)

!      else
!
!       if(optimise) then
!         write(outfile,'(I5,".",I5,".OPT_UX.it",I3.3) ') &
!              nrx(ir)-lmargin(1),nrz(ir)-lmargin(2),iterationIndex
!       else
!         write(outfile,'(I5,".",I5,".CON_UX.it",I3.3) ') &
!              nrx(ir)-lmargin(1),nrz(ir)-lmargin(2),iterationIndex
!       endif
!
!       do j=1,12
!          if(outfile(j:j).eq.' ') outfile(j:j)='0'
!       enddo
!
!       outfile = './synthetics/'//trim(modelname)//'/'//outfile
!       open(1,file=outfile,form='unformatted',access='direct',recl=recl_size_syn)
!       write(1,rec=1) synx(0:maxnt,ir)
!       close(1)

      enddo
     endif


    if(iterationIndex.eq.0) then
      do ir=1,nReceiver
        if(optimise)then
          write(outfile,'(I5,".",I5,".OPT_UZ") ') &
                nrx(ir)-lmargin(1),nrz(ir)-lmargin(2)
          else
          write(outfile,'(I5,".",I5,".CON_UZ") ') &
                nrx(ir)-lmargin(1),nrz(ir)-lmargin(2)
        endif

        do j=1,12
           if(outfile(j:j).eq.' ') outfile(j:j)='0'
        enddo

        outfile = './synthetics/'//trim(modelname)//'/'//outfile

        open(1,file=outfile,form='unformatted',access='direct',recl=recl_size_syn)
        write(1,rec=1) synz(0:maxnt,ir)
        close(1)

!     else
!
!        if(optimise) then
!          write(outfile,'(I5,".",I5,".OPT_UZ.it",I3.3) ') &
!                nrx(ir)-lmargin(1),nrz(ir)-lmargin(2),iterationIndex
!          else
!          write(outfile,'(I5,".",I5,".CON_UZ.it",I3.3) ') &
!                nrx(ir)-lmargin(1),nrz(ir)-lmargin(2),iterationIndex
!          endif
!
!        do j=1,12
!        if(outfile(j:j).eq.' ') outfile(j:j)='0'
!        enddo
!
!        outfile = './synthetics/'//trim(modelname)//'/'//outfile
!        open(1,file=outfile,form='unformatted',access='direct',recl=recl_size_syn)
!        write(1,rec=1) synz(0:maxnt,ir)
!        close(1)

      enddo
    endif

!!%!! For Plotting Strain for UP & US

    !%! For UP

    if(iterationIndex.eq.0) then
    if(optimise) then
    write(outfile,'(I5,".",I5,".OPT_UP") ')  &
    isx-lmargin(1),isz-lmargin(2)
    else
    write(outfile,'(I5,".",I5,".CON_UP") ') &
    isx-lmargin(1),isz-lmargin(2)
    endif

    do j=1,12
    if(outfile(j:j).eq.' ') outfile(j:j)='0'
    enddo

    outfile = './synthetics/'//trim(modelname)//'/'//outfile
    else
    if(optimise) then
    write(outfile,'(I5,".",I5,".OPT_UP.it",I3.3) ')  &
    isx-lmargin(1),isz-lmargin(2),iterationIndex
    else
    write(outfile,'(I5,".",I5,".CON_UP,it",I3.3) ') &
    isx-lmargin(1),isz-lmargin(2),iterationIndex
    endif

    do j=1,12
    if(outfile(j:j).eq.' ') outfile(j:j)='0'
    enddo

    outfile = './synthetics/'//trim(modelname)//'/'//outfile

    endif




    open(1,file=outfile,form='unformatted',access='direct',recl=recl_size_syn)
    write(1,rec=1) synp(0:maxnt,1:nReceiver)
    close(1)

    !%! For US

    if(iterationIndex.eq.0) then
    if(optimise) then
    write(outfile,'(I5,".",I5,".OPT_US") ')  &
    isx-lmargin(1),isz-lmargin(2)
    else
    write(outfile,'(I5,".",I5,".CON_US") ') &
    isx-lmargin(1),isz-lmargin(2)
    endif

    do j=1,12
    if(outfile(j:j).eq.' ') outfile(j:j)='0'
    enddo

    outfile = './synthetics/'//trim(modelname)//'/'//outfile
    else
    if(optimise) then
    write(outfile,'(I5,".",I5,".OPT_US.it",I3.3) ')  &
    isx-lmargin(1),isz-lmargin(2),iterationIndex
    else
    write(outfile,'(I5,".",I5,".CON_US,it",I3.3) ') &
    isx-lmargin(1),isz-lmargin(2),iterationIndex
    endif

    do j=1,12
    if(outfile(j:j).eq.' ') outfile(j:j)='0'
    enddo

    outfile = './synthetics/'//trim(modelname)//'/'//outfile

    endif




open(1,file=outfile,form='unformatted',access='direct',recl=recl_size_syn)
write(1,rec=1) syns(0:maxnt,1:nReceiver)
close(1)




  enddo
end subroutine forwardmodelling
