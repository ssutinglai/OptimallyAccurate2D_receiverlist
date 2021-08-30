set(0,'defaultfigurecolor','w');
set(0,'defaultaxesfontname','Times new roman')
set(0,'defaultaxesfontsize',18)

%%
Lx = 810;
Lz = 600;
dx = 0.8;
dz = 0.8;
v1 = 800;
v2 = 1750;
v3 = 800;
v4 = 1750;
v5 = 800;
%%
nx = ceil(Lx/dx);
nz = ceil(Lz/dz)+1;

%%
vel = 10+zeros(nz,nx);  %VP
vels = 0+zeros(nz,nx);  %VS 
rho  = 50+zeros(nz,nx);  %rho
%% layer 1
p1 = [110 0];
p2 = [85 190];
p3 = [80 210];
p4 = [45 440];
p5 = [40 450];
p6 = [0 700];
p7 = [110 810];
layer1 = zeros(1, nx);
%
k1 = (p1(1)-p2(1))/(p1(2)-p2(2));
ii = 0;
x  = p1(2);
while(x<=p2(2))
    ii = ii+1;
    layer1(ii) = p1(1)+k1*(x-p1(2));
    x = x+dx;
end
%
k2 = (p2(1)-p3(1))/(p2(2)-p3(2));
x  = p2(2);
ii = ii-1;
while(x<=p3(2))
    ii = ii+1;
    layer1(ii) = p2(1)+k2*(x-p2(2));
    x = x+dx;
end
%
k3 = (p3(1)-p4(1))/(p3(2)-p4(2));
x  = p3(2);
ii = ii-1;
while(x<=p4(2))
    ii = ii+1;
    layer1(ii) = p3(1)+k3*(x-p3(2));
    x = x+dx;
end
%
k4 = (p4(1)-p5(1))/(p4(2)-p5(2));
x  = p4(2);
ii = ii-1;
while(x<=p5(2))
    ii = ii+1;
    layer1(ii) = p4(1)+k4*(x-p4(2));
    x = x+dx;
end
%
k5 = (p5(1)-p6(1))/(p5(2)-p6(2));
x  = p5(2);
ii = ii-1;
while(x<=p6(2))
    ii = ii+1;
    layer1(ii) = p5(1)+k5*(x-p5(2));
    x = x+dx;
end
%
k6 = (p6(1)-p7(1))/(p6(2)-p7(2));
x  = p6(2);
ii = ii-1;
while(x<=p7(2))
    ii = ii+1;
    layer1(ii) = p6(1)+k6*(x-p6(2));
    x = x+dx;
end
% figure()
%     plot(layer1)

%% layer 2
p1 = [110 0];
p2 = [85 190];
p3 = [80 210];
p4 = [45 440];
p5 = [40 450];
p6 = [30 550];
p7 = [90 710];
p8 = [90 790];
p9 = [110 810];

layer2 = zeros(1, nx);
%
k1 = (p1(1)-p2(1))/(p1(2)-p2(2));
ii = 0;
x  = p1(2);
while(x<=p2(2))
    ii = ii+1;
    layer2(ii) = p1(1)+k1*(x-p1(2));
    x = x+dx;
end
%
k2 = (p2(1)-p3(1))/(p2(2)-p3(2));
x  = p2(2);
ii = ii-1;
while(x<=p3(2))
    ii = ii+1;
    layer2(ii) = p2(1)+k2*(x-p2(2));
    x = x+dx;
end
%
k3 = (p3(1)-p4(1))/(p3(2)-p4(2));
x  = p3(2);
ii = ii-1;
while(x<=p4(2))
    ii = ii+1;
    layer2(ii) = p3(1)+k3*(x-p3(2));
    x = x+dx;
end
%
k4 = (p4(1)-p5(1))/(p4(2)-p5(2));
x  = p4(2);
ii = ii-1;
while(x<=p5(2))
    ii = ii+1;
    layer2(ii) = p4(1)+k4*(x-p4(2));
    x = x+dx;
end
%
k5 = (p5(1)-p6(1))/(p5(2)-p6(2));
x  = p5(2);
ii = ii-1;
while(x<=p6(2))
    ii = ii+1;
    layer2(ii) = p5(1)+k5*(x-p5(2));
    x = x+dx;
end
%
k6 = (p6(1)-p7(1))/(p6(2)-p7(2));
x  = p6(2);
ii = ii-1;
while(x<=p7(2))
    ii = ii+1;
    layer2(ii) = p6(1)+k6*(x-p6(2));
    x = x+dx;
end
%
k7 = (p7(1)-p8(1))/(p7(2)-p8(2));
x  = p7(2);
ii = ii-1;
while(x<=p8(2))
    ii = ii+1;
    layer2(ii) = p7(1)+k7*(x-p7(2));
    x = x+dx;
end
%
k8 = (p8(1)-p9(1))/(p8(2)-p9(2));
x  = p8(2);
ii = ii-1;
while(x<=p9(2))
    ii = ii+1;
    layer2(ii) = p8(1)+k8*(x-p8(2));
    x = x+dx;
end
% figure()
%     plot(layer2)

%% layer 3
p1 = [110 0];
p2 = [85 190];
p3 = [80 210];
p4 = [45 440];
p5 = [35 550];
p6 = [95 710];
p7 = [95 795];
p8 = [110 810];
layer3 = zeros(1, nx);
%
k1 = (p1(1)-p2(1))/(p1(2)-p2(2));
ii = 0;
x  = p1(2);
while(x<=p2(2))
    ii = ii+1;
    layer3(ii) = p1(1)+k1*(x-p1(2));
    x = x+dx;
end
%
k2 = (p2(1)-p3(1))/(p2(2)-p3(2));
x  = p2(2);
ii = ii-1;
while(x<=p3(2))
    ii = ii+1;
    layer3(ii) = p2(1)+k2*(x-p2(2));
    x = x+dx;
end
%
k3 = (p3(1)-p4(1))/(p3(2)-p4(2));
x  = p3(2);
ii = ii-1;
while(x<=p4(2))
    ii = ii+1;
    layer3(ii) = p3(1)+k3*(x-p3(2));
    x = x+dx;
end
%
k4 = (p4(1)-p5(1))/(p4(2)-p5(2));
x  = p4(2);
ii = ii-1;
while(x<=p5(2))
    ii = ii+1;
    layer3(ii) = p4(1)+k4*(x-p4(2));
    x = x+dx;
end
%
k5 = (p5(1)-p6(1))/(p5(2)-p6(2));
x  = p5(2);
ii = ii-1;
while(x<=p6(2))
    ii = ii+1;
    layer3(ii) = p5(1)+k5*(x-p5(2));
    x = x+dx;
end
%
k6 = (p6(1)-p7(1))/(p6(2)-p7(2));
x  = p6(2);
ii = ii-1;
while(x<=p7(2))
    ii = ii+1;
    layer3(ii) = p6(1)+k6*(x-p6(2));
    x = x+dx;
end
%
k7 = (p7(1)-p8(1))/(p7(2)-p8(2));
x  = p7(2);
ii = ii-1;
while(x<=p8(2))
    ii = ii+1;
    layer3(ii) = p7(1)+k7*(x-p7(2));
    x = x+dx;
end
% figure()
%     plot(layer3)

%% layer 4
p1 = [110 0];
p2 = [85 190];
p3 = [80 210];
p4 = [65 410];
p5 = [85 480];
p6 = [100 800];
p7 = [110 810];
layer4 = zeros(1, nx);
%
k1 = (p1(1)-p2(1))/(p1(2)-p2(2));
ii = 0;
x  = p1(2);
while(x<=p2(2))
    ii = ii+1;
    layer4(ii) = p1(1)+k1*(x-p1(2));
    x = x+dx;
end
%
k2 = (p2(1)-p3(1))/(p2(2)-p3(2));
x  = p2(2);
ii = ii-1;
while(x<=p3(2))
    ii = ii+1;
    layer4(ii) = p2(1)+k2*(x-p2(2));
    x = x+dx;
end
%
k3 = (p3(1)-p4(1))/(p3(2)-p4(2));
x  = p3(2);
ii = ii-1;
while(x<=p4(2))
    ii = ii+1;
    layer4(ii) = p3(1)+k3*(x-p3(2));
    x = x+dx;
end
%
k4 = (p4(1)-p5(1))/(p4(2)-p5(2));
x  = p4(2);
ii = ii-1;
while(x<=p5(2))
    ii = ii+1;
    layer4(ii) = p4(1)+k4*(x-p4(2));
    x = x+dx;
end
%
k5 = (p5(1)-p6(1))/(p5(2)-p6(2));
x  = p5(2);
ii = ii-1;
while(x<=p6(2))
    ii = ii+1;
    layer4(ii) = p5(1)+k5*(x-p5(2));
    x = x+dx;
end
%
k6 = (p6(1)-p7(1))/(p6(2)-p7(2));
x  = p6(2);
ii = ii-1;
while(x<=p7(2))
    ii = ii+1;
    layer4(ii) = p6(1)+k6*(x-p6(2));
    x = x+dx;
end
% figure()
%     plot(layer4)   
    
%% layer 5
p1 = [110 0];
p2 = [85 190];
p3 = [70 410];
p4 = [93 480];
p5 = [105 805];
p6 = [110 810];
layer5 = zeros(1, nx);
%
k1 = (p1(1)-p2(1))/(p1(2)-p2(2));
ii = 0;
x  = p1(2);
while(x<=p2(2))
    ii = ii+1;
    layer5(ii) = p1(1)+k1*(x-p1(2));
    x = x+dx;
end
%
k2 = (p2(1)-p3(1))/(p2(2)-p3(2));
x  = p2(2);
ii = ii-1;
while(x<=p3(2))
    ii = ii+1;
    layer5(ii) = p2(1)+k2*(x-p2(2));
    x = x+dx;
end
%
k3 = (p3(1)-p4(1))/(p3(2)-p4(2));
x  = p3(2);
ii = ii-1;
while(x<=p4(2))
    ii = ii+1;
    layer5(ii) = p3(1)+k3*(x-p3(2));
    x = x+dx;
end
%
k4 = (p4(1)-p5(1))/(p4(2)-p5(2));
x  = p4(2);
ii = ii-1;
while(x<=p5(2))
    ii = ii+1;
    layer5(ii) = p4(1)+k4*(x-p4(2));
    x = x+dx;
end
%
k5 = (p5(1)-p6(1))/(p5(2)-p6(2));
x  = p5(2);
ii = ii-1;
while(x<=p6(2))
    ii = ii+1;
    layer5(ii) = p5(1)+k5*(x-p5(2));
    x = x+dx;
end
%

% figure()
%     plot(layer5)       
%%
for i=1:1:nx
    vel(ceil(layer1(i)/dz)+1:ceil(layer2(i)/dz),i) = v1;
    vels(ceil(layer1(i)/dz)+1:ceil(layer2(i)/dz),i) = 300;
    rho(ceil(layer1(i)/dz)+1:ceil(layer2(i)/dz),i) = 1600;
    
end
for i=1:1:nx
    vel(ceil(layer2(i)/dz)+1:nz,i) = v2;
    vels(ceil(layer2(i)/dz)+1:nz,i) = 500;
    rho(ceil(layer2(i)/dz)+1:nz,i) = 2000;
end
for i=1:1:nx
    vel(ceil(layer3(i)/dz)+1:nz,i) = v3;
    vels(ceil(layer3(i)/dz)+1:nz,i) = 300;
    rho(ceil(layer3(i)/dz)+1:nz,i) = 1600;
end
for i=1:1:nx
    vel(ceil(layer4(i)/dz)+1:nz,i) = v4;
    vels(ceil(layer4(i)/dz)+1:nz,i) = 500;
    rho(ceil(layer4(i)/dz)+1:nz,i) = 2000;
end
for i=1:1:nx
    vel(ceil(layer5(i)/dz)+1:nz,i) = v5;
    vels(ceil(layer5(i)/dz)+1:nz,i) = 300;
    rho(ceil(layer5(i)/dz)+1:nz,i) = 1600;
end
figure()
    imagesc(vel); 
%     axis([0 405 0 00])
figure()
    imagesc(vels); 
%     axis([0 405 0 100])
figure()
    imagesc(rho); 
%     axis([0 405 0 100])
    
vel = vel(1:end-1,1:end);
vels = vels(1:end-1,1:end);
rho = rho(1:end-1,1:end);
%% Output the model file as Opt2D format
fid=fopen('dunepilat.vp','w');
fwrite(fid,vel','single');
fclose(fid);
fid=fopen('dunepilat.vs','w');
fwrite(fid,vels','single');
fclose(fid);
fid=fopen('dunepilat.rho','w');
fwrite(fid,rho','single');
fclose(fid);

%% Output receiver file as Opt2D format
nreceiver = 10;    %number of the receivers 500
dreceiver = 100;      %distances between the receivers 2
receiver0 = 2; 
receiver = zeros(nreceiver,2);
for i = 1:1:nreceiver
    receiver(i,1)= (i-1)*dreceiver+receiver0; 
    receiver(i,2)= floor(layer1(receiver(i,1))/dz)+1;
end
fid = fopen('../inffile/receiver_depth.txt','w');
for i=1:1:nreceiver-1
    fprintf(fid,'%d\t%d\n',receiver(i,1),receiver(i,2));
end
fprintf(fid,'%d\t%d',receiver(nreceiver,1),receiver(nreceiver,2));
fclose(fid);
