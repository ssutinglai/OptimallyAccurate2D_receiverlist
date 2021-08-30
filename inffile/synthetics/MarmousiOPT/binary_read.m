clear all;clc;close all
nt=6000;
nr = 400;
a = load('/Users/alice/Desktop/pilat_dune/MODEL/receiver_depth.txt');
syn = zeros(nt,nr);

for iReceiver=1:nr
%     iReceiver
    num5digit=sprintf('%05d',a(iReceiver,1));
    num2digit=sprintf('%05d',a(iReceiver,2));
    filename=strcat(num5digit,'.',num2digit,'.CON_UX');
    
    fileID1=fopen(filename,'r');
    AA=fread(fileID1,'single'); %For reading kind(1.e0)
    seismogram=AA(1:nt);
    fclose(fileID1);
    syn(:,iReceiver) = seismogram(1:end);
%     A=[time',seismogram];
%     
%     %
%     fileID = fopen(filename2,'w');
%     fprintf(fileID,'%e %e\n',A');
%     fclose(fileID);
%     iReceiver
end
perc= 0.1;
syn1 = (syn);
max_syn=max(max(abs(syn)));
for i=1:1:nt
    for j=1:1:nr
        if(abs(syn(i,j))>perc*max_syn)
            syn1(i,j) = syn(i,j)*perc;
        end
    end
end
figure()
    imagesc(syn1);   colormap(gray)
    