clear all

for iReceiver=1;
    iReceiver
    filename=strcat('dunepilat.vp');
    filename2=strcat('modelVP.txt');
    
    fileID1=fopen(filename,'r');
    AA=fread(fileID1,'single'); %For reading kind(1.e0)
    fclose(fileID1);
    
    %
    fileID = fopen(filename2,'w');
    fprintf(fileID,'%e\n',AA');
    fclose(fileID);
    iReceiver
end

for iReceiver=1;
    iReceiver
    filename=strcat('dunepilat.vs');
    filename2=strcat('modelVS.txt');
    
    fileID1=fopen(filename,'r');
    AA=fread(fileID1,'single'); %For reading kind(1.e0)
    fclose(fileID1);
    
    %
    fileID = fopen(filename2,'w');
    fprintf(fileID,'%e\n',AA');
    fclose(fileID);
    iReceiver
end

for iReceiver=1;
    iReceiver
    filename=strcat('dunepilat.rho');
    filename2=strcat('modelRHO.txt');
    
    fileID1=fopen(filename,'r');
    AA=fread(fileID1,'single'); %For reading kind(1.e0)
    fclose(fileID1);
    
    %
    fileID = fopen(filename2,'w');
    fprintf(fileID,'%e\n',AA');
    fclose(fileID);
    iReceiver
end