# OptimallyAccurate2D_receiverlist

This is a 2D finite difference code (Opt2D) for modeling wave propagtion written by Nobuaki FUJI and modeified by Ssu-Ting LAI.
"OptimallyAccurate2D_receiverlist" version of Opt2D code is specially designed for Matlab users for generating models and plotting the waveforms using Matlab.

For the users,

1. There are mainly two parts of the code to take care: 1) models and 2) file of parameter setting.

2. models/modelgenerator.m is an example Matlab file for generating Opt2D Vp, Vs and rho files input. Also, the file models/modelgenerator.m will also generate a list of receivers (inffile/receiver_depth.txt) on the top of the models.

3. inffile/para.inf is the most important file for controlling the parameters of wave propagation. The parameters in the para.inf file are list below:
   (Line 1).  MarmousiOPT => The name of the project.
   (Line 2).  ../models/dunepilat.vp  => The address to the Vp file.
   (Line 3).  ../models/dunepilat.vs  => The address to the Vs file.
   (Line 4).  ../models/dunepilat.rho => The address to the who file.
   (Line 5).  F optimise or conventional => If T, the code will process in optimise mode. If F, the code use conventional mode.
   (Line 6).  T video or not  => If T, the code will make all the snapshots into a video. (You may need to download "imagemagick" in your computer to do it.)
   (Line 7).  T writingStrain => If T, the code will output strain files into inffile/strains.
   (Line 8).  10 IT_DISPLAY   => The frame rate of the video.  
   (Line 9).  550 10 1  iSourceStart, iSourceInverval, nSource  => The starting point of the first source (X-component)/ Interval between the sources/ The number of the sources you want to assign.
   (Line10).  56 izSourceStart                                  => The starting point of the first source (Z-component)
   (Line11).  110 1 10 iReceiverStart, iReceiverInterval, nReceiver => The starting point of the first receiver (X-component)/ Interval between the receivers/ The number of the receivers you want to assign. However, YOU DON'T NEED TO TAKE CARE OF "iReceiverStart" AND "iReceiverInterval" IN THIS VERSION. YOU ONLY NEED TO PUT THE NUMBER OF RECEIVERS, WHICH IS "nReceiver". THE COORDINATES OF THE RECEIVERS ARE GIVEN BY "inffile/receiver_depth.txt".
   (Line12).  110 izReceiverStart                                   => The starting point of the first receivers (Z-component)
   (Line13).  6000 1012 250 nt nx nz                                => The number of time steps(nt)/The numbers of points in X-direction(nx)/The numbers of points in Z-direction(nz)--When you design your model, it is better to make extra length in Z-direction.
   (Line14).  1.d-4 8.d-4 8.d-4 dt,dx,dz                            => The time step(dt)/The step of X-direction in kilometers(dx)/The step of Z-direction in kilometers(dz)       
   (Line15).  1.d2 0 f0, t0                                         => The dominant frequency of Ricker wavelet(f0)/ Time 0(t0)
   (Line16).  1  1 i1Source, i2Source ! not necessary for FWI
   (Line17).  ~/OptimallyAccurate2D/Marmousi/synthetics/MarmousiOPT/tmp/
   (Line18).  OPT_UX
   (Line19).  OPT_UZ
   (Line20).  0
   (Line21).  1.e-3 1.e-3
   
4. A easy guide for running the code:
   (1) cd ./models
   (2) run matlab code modelgenerator.m and create dunepilat.vp, dunepilat.vs, dunepilat.rho and inffile/receiver_depth.txt 
   (3) cd ../inffile
   (4) modify para.inf to the parameters you want
   (5) ../bin/FWIfull.x < para.inf
   (6) look at the results in inffile/snapshots, inffile/synthetics and inffile/videos/MarmousiOPT
   (7) use Matlab file inffile/synthetics/binary_read.m to plot the waveforms
