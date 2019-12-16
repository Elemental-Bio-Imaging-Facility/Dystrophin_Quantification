clear
close all
clc

%% Import
delimiterIn = ',';
headerlinesIn = 1;
Directory='D:\';
mkdir ([Directory], 'OutDirectory');
OutDirectory=([Directory 'OutDirectory\']);
mkdir ([Directory], 'SliceDirectory');
SliceDirectory=([Directory 'SliceDirectory\']);

%% Files and initial conditions for loop
HFiles=dir([Directory '*horz*.csv']);
VFiles=dir([Directory '*vert*.csv']);
lengthV=length(VFiles);
lengthH=length(HFiles);
layer=0; %layers for counting import layers of 3D matrix

%% Constants
Speed=15;%in um/s
Scantime=0.125;%in seconds
SpotSize=15;%in um
gradient=1;%from calibraiton curve
intercept=0;%from calibraiton curve
Frequency=20;%Hz
M=SpotSize/(Speed*Scantime); %Magnification factor (ratio of speed and scan time)
N=floor(M/2);
if mod(N,2)==0
    N=N+1;
end
%%N=floor(M/2)+1; %Padding offset
Layers=2; %layer of ablation offset
TheoryRes=Speed/Frequency;%Physical resolution of raster regime
ExpRes=Speed/(M*2);%Physical resolution of raster regime


%% Parameters
dataH = importdata([Directory HFiles(1).name],delimiterIn,headerlinesIn);
dataV = importdata([Directory VFiles(1).name],delimiterIn,headerlinesIn);

la1=size(dataH.data,1);
wa=size(dataH.data,2);
lb1=size(dataV.data,1);
wb=size(dataV.data,2);
j=round(((1.2/Scantime)*10));%removed from top (M*laser warmup time or just speed times 10s)DONT FORGET BEGINNING IS ALWAYS LONGER normally about 12*M
k=round((la1-(wa*M))-j);%removed from bottom (M times the lines of the other matrix - large number) leftright
l=round((lb1-(wb*M))-j);%removed from bottom (M times the lines of the other matrix - large number) updown 
A1=dataH.data;
A1(la1-l+1:la1, :) = []; %remove k rows from the bottom
A1(1:j, :) = []; %remove j rows from the top
updowncut1=A1;
sizex=size(updowncut1,1)*Layers+N; % size of y dimension including upsample, kron and padding
sizey=size(updowncut1,2)*M*Layers+N;  % size of x dimension including upsample, kron and padding
lengthz=length(dir([Directory '*.csv'])); % size of z dimension
FirstLayer=kron(dataH.data,ones(1,M)); %First layer converted to dimensions

%% Allocation
RM=zeros(sizey,sizex,lengthz);%to include padding and kron

%% Loop of all data manipulation

for ifp=1:lengthV
    fileH = importdata([Directory HFiles(ifp).name],delimiterIn,headerlinesIn);
    fileV = importdata([Directory VFiles(ifp).name],delimiterIn,headerlinesIn);
    dataH=fileH.data;
    dataV=fileV.data;
    
    dataH(la1-l+1:la1, :) = []; %remove k rows from the bottom
    dataH(1:j, :) = []; %remove j rows from the top
    
    dataH=(dataH-intercept)/gradient; %convert CPS to ppm,ppb
    dataH=kron(dataH,ones(1,M)); %Duplicate M number columns
    dataH=upsample(dataH,Layers); %Fill Layer number zero columns
    dataH=transpose(dataH); %Rotate to upsample
    dataH=upsample(dataH,Layers); %Fill M number zero columns
    
    dataH=padarray(dataH,[N N],0,'post'); %Pad last row and column with zeros for 3D
    
    dataV(lb1-k+1:lb1, :) = []; %remove k rows from the bottom
    dataV(1:j, :) = []; %remove j rows from the top
    dataV=(dataV-intercept)/gradient; %convert CPS to ppm,ppb
    dataV=kron(dataV,ones(1,M)); %Duplicate M number columns
    dataV=upsample(dataV,Layers); %Fill Layer number zero columns
    dataV=transpose(dataV); %Rotate to match
    dataV=upsample(dataV,Layers); %Fill M number zero columns
    dataV=transpose(dataV); %Rotate to match
    
    dataV=padarray(dataV,[N N],0,'pre'); %Pad first row and column with zeros
    
    layer=layer+1;
    
    RM(:,:,layer) = dataH;
    
    layer=layer+1;
    
    RM(:,:,layer) = dataV;
    
    
end

if lengthV<lengthH
    
    dataH = importdata([Directory HFiles(lengthH).name],delimiterIn,headerlinesIn);
      
    dataH(la1-l+1:la1, :) = []; %remove k rows from the bottom
    dataH(1:j, :) = []; %remove j rows from the top
    
    dataH=(dataH-intercept)/gradient; %convert CPS to ppm,ppb
    dataH=kron(dataH,ones(1,M)); %Duplicate M number columns
    dataH=upsample(dataH,Layers); %Fill Layer number zero columns
    dataH=transpose(dataH); %Rotate to upsample
    dataH=upsample(dataH,Layers); %Fill M number zero columns
    
    dataH=padarray(dataH,[N N],0,'post'); %Pad last row and column with zeros for 3D
    
    layer=layer+1;
    
    RM(:,:,layer) = dataH;
    
end
    
RM(RM == 0) = NaN;%replace zeroes with NaN

%% Trilinear interpolation
FM = ones(size(RM));
MM = ones(size(RM));
MM = MM.*RM;
[NcomX, NcomY, NcomZ] = size(RM);

for iz=1:NcomZ
    for ic=1:NcomY
        for ir=1:NcomX
            if isnan(MM(ir,ic,iz))
                cf=0;
                avenum=0;
                if (iz>1)&&(~isnan(MM(ir,ic,iz-1)))
                    cf=cf+MM(ir,ic,iz-1);
                    avenum=avenum+1;
                end
                if (ir-1>0)&&(iz>1)&&(~isnan(MM(ir-1,ic,iz-1)))
                    cf=cf+MM(ir-1,ic,iz-1);
                    avenum=avenum+1;
                end
                if (ir+1<=NcomX)&&(iz>1)&&(~isnan(MM(ir+1,ic,iz-1)))
                    cf=cf+MM(ir+1,ic,iz-1);
                    avenum=avenum+1;
                end
                if (ic-1>0)&&(iz>1)&&(~isnan(MM(ir,ic-1,iz-1)))
                    cf=cf+MM(ir,ic-1,iz-1);
                    avenum=avenum+1;
                end
                if (ic+1<=NcomY)&&(iz>1)&&(~isnan(MM(ir,ic+1,iz-1)))
                    cf=cf+MM(ir,ic+1,iz-1);
                    avenum=avenum+1;
                end
                if (ir-1>0)&&(ic-1>0)&&(iz>1)&&(~isnan(MM(ir-1,ic-1,iz-1)))
                    cf=cf+MM(ir-1,ic-1,iz-1);
                    avenum=avenum+1;
                end
                if (ir+1<=NcomX)&&(ic+1<=NcomY)&&(iz>1)&&(~isnan(MM(ir+1,ic+1,iz-1)))
                    cf=cf+MM(ir+1,ic+1,iz-1);
                    avenum=avenum+1;
                end
                if (ir+1<=NcomX)&&(ic-1>0)&&(iz>1)&&(~isnan(MM(ir+1,ic-1,iz-1)))
                    cf=cf+MM(ir+1,ic-1,iz-1);
                    avenum=avenum+1;
                end
                if (ir-1>0)&&(ic+1<=NcomY)&&(iz>1)&&(~isnan(MM(ir-1,ic+1,iz-1)))
                    cf=cf+MM(ir-1,ic+1,iz-1);
                    avenum=avenum+1;
                end
                %from layer-1 to layer
                if (ir-1>0)&&(~isnan(MM(ir-1,ic,iz)))
                    cf=cf+MM(ir-1,ic,iz);
                    avenum=avenum+1;
                end
                if (ir+1<=NcomX)&&(~isnan(MM(ir+1,ic,iz)))
                    cf=cf+MM(ir+1,ic,iz);
                    avenum=avenum+1;
                end
                if (ic-1>0)&&(~isnan(MM(ir,ic-1,iz)))
                    cf=cf+MM(ir,ic-1,iz);
                    avenum=avenum+1;
                end
                if (ic+1<=NcomY)&&(~isnan(MM(ir,ic+1,iz)))
                    cf=cf+MM(ir,ic+1,iz);
                    avenum=avenum+1;
                end
                if (ir-1>0)&&(ic-1>0)&&(~isnan(MM(ir-1,ic-1,iz)))
                    cf=cf+MM(ir-1,ic-1,iz);
                    avenum=avenum+1;
                end
                if (ir+1<=NcomX)&&(ic+1<=NcomY)&&(~isnan(MM(ir+1,ic+1,iz)))
                    cf=cf+MM(ir+1,ic+1,iz);
                    avenum=avenum+1;
                end
                if (ir+1<=NcomX)&&(ic-1>0)&&(~isnan(MM(ir+1,ic-1,iz)))
                    cf=cf+MM(ir+1,ic-1,iz);
                    avenum=avenum+1;
                end
                if (ir-1>0)&&(ic+1<=NcomY)&&(~isnan(MM(ir-1,ic+1,iz)))
                    cf=cf+MM(ir-1,ic+1,iz);
                    avenum=avenum+1;
                end
                %from layer to layer+1
                if (iz<NcomZ)&&(~isnan(MM(ir,ic,iz+1)))
                    cf=cf+MM(ir,ic,iz+1);
                    avenum=avenum+1;
                end
                if (ir-1>0)&&(iz<NcomZ)&&(~isnan(MM(ir-1,ic,iz+1)))
                    cf=cf+MM(ir-1,ic,iz+1);
                    avenum=avenum+1;
                end
                if (ir+1<=NcomX)&&(iz<NcomZ)&&(~isnan(MM(ir+1,ic,iz+1)))
                    cf=cf+MM(ir+1,ic,iz+1);
                    avenum=avenum+1;
                end
                if (ic-1>0)&&(iz<NcomZ)&&(~isnan(MM(ir,ic-1,iz+1)))
                    cf=cf+MM(ir,ic-1,iz+1);
                    avenum=avenum+1;
                end
                if (ic+1<=NcomY)&&(iz<NcomZ)&&(~isnan(MM(ir,ic+1,iz+1)))
                    cf=cf+MM(ir,ic+1,iz+1);
                    avenum=avenum+1;
                end
                if (ir-1>0)&&(ic-1>0)&&(iz<NcomZ)&&(~isnan(MM(ir-1,ic-1,iz+1)))
                    cf=cf+MM(ir-1,ic-1,iz+1);
                    avenum=avenum+1;
                end
                if (ir+1<=NcomX)&&(ic+1<=NcomY)&&(iz<NcomZ)&&(~isnan(MM(ir+1,ic+1,iz+1)))
                    cf=cf+MM(ir+1,ic+1,iz+1);
                    avenum=avenum+1;
                end
                if (ir+1<=NcomX)&&(ic-1>0)&&(iz<NcomZ)&&(~isnan(MM(ir+1,ic-1,iz+1)))
                    cf=cf+MM(ir+1,ic-1,iz+1);
                    avenum=avenum+1;
                end
                if (ir-1>0)&&(ic+1<=NcomY)&&(iz<NcomZ)&&(~isnan(MM(ir-1,ic+1,iz+1)))
                    cf=cf+MM(ir-1,ic+1,iz+1);
                    avenum=avenum+1;
                end
                         
                FM(ir,ic,iz)=cf/avenum;
            else
                FM(ir,ic,iz)=MM(ir,ic,iz);
            end
        end
    end
end



%% Inpaint & Assign
FM(isnan(FM))=0; %makes pad NaN equal zero
InterpBig3D=FM; %Matrix to old term
Interp2D=sum(InterpBig3D,3);%turn 3D matrix into 2D by summing in z direction
Mn=round(((2*M)-1)/2)*2+1; %round M to next odd up
Gaussian=imgaussfilt(Interp2D,1,'FilterSize',Mn);

%% Split layers
for ilayers=1:lengthz
    name = sprintf('%sSlice00_%d_Gd158SRR.csv',SliceDirectory,ilayers);
    dlmwrite(name,FM(:, :, ilayers));
end

%% Writing Files
dlmwrite([OutDirectory 'SRRCombined3D.csv'],InterpBig3D);
dlmwrite([OutDirectory 'SRRCombined2D.csv'],Interp2D);
dlmwrite([OutDirectory 'SRRGaussian.csv'],Gaussian);
dlmwrite([OutDirectory 'SRRConventional.csv'],FirstLayer);

InterpBin=InterpBig3D;

[x y z]=size(InterpBin);
filename='3DSRR.vtk'; 
total=numel(InterpBig3D);
z=size(InterpBig3D,3); %layer of ablation
fp=fopen(sprintf('%s',filename),'w'); 
fp=fopen([OutDirectory '3DSRR.vtk'],'w');
fprintf(fp,sprintf('# vtk DataFile Version 2.0\n%s\nASCII\nDATASET STRUCTURED_POINTS\nDIMENSIONS %d %d %d\nORIGIN 0.0 0.0 0.0\nSPACING 1.0 1.0 5.0\n \nPOINT_DATA %d\nSCALARS values double\nLOOKUP_TABLE default\n',filename,x,y,z,total)');
for z1=1:z;
for y1=1:y;
for    x1=1:x;
    if x1>=x;
        fprintf(fp,'%d\n',InterpBin(x1,y1,z1));
    else
    fprintf(fp,'%d\t',InterpBin(x1,y1,z1));
    end
            
end
end
end
fclose(fp);
