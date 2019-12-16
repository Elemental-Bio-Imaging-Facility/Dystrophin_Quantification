clear
close all
clc

%% thresholding

dname='D:\';
files=dir([dname,'*.txt']);
for i=1:length(files);
    filename=files(i).name;
    csvdata=importdata([dname,files(i).name]);

    csvdata(csvdata<0)=0; %remove all negative values
    count=i;
    rawmax=max(reshape(csvdata,1,[]));
    rawmean=mean(reshape(csvdata,1,[]));
    rawmode=mode(reshape(csvdata,1,[]));
    maximum=max(reshape(csvdata,1,[])); %max value in array
    median1=nanmedian(nanmedian(csvdata,1),2);
    csvdata1=csvdata>median1; %above median (aka background) signal
    csvdataNZ=csvdata>0; %above zero (aka background) signal
    csvdatamatrix=(reshape(csvdata(csvdata1),1,[])); %check number of values above median
    average=mean(reshape(csvdata(csvdata1),1,[]));  %average of values greater than zero
    rsd=(std(reshape(csvdata(csvdata1),1,[]))/mean(reshape(csvdata(csvdata1),1,[])))*100;
    thresh = multithresh(csvdata(csvdataNZ)); %otsu's method thresholding
    thresh3 = multithresh(csvdata(csvdataNZ),3); %otsu's method thresholding
    [cluster_idx, cluster_center] = kmeans(csvdata,2,'distance','sqEuclidean'); %2 level k mean clustering
    threshmean=csvdata>thresh; %greater than ostu threshold values
    thresh2 = multithresh(csvdata(threshmean)); %otsu's method thresholding
    threshthreshmean=csvdata>thresh2; %greater than double ostu threshold values
    threshthreshmean1 = mean(reshape(csvdata(threshthreshmean),1,[])); %2nd otsu's method thresholding
    threshmatrix=(reshape(csvdata(threshmean),1,[])); %check number of values above otsu
    threshmedian=median(reshape(csvdata(threshmean),1,[])); % median of values above otsu threshold
    countall=numel(csvdata);%count all pixels
    countmedian=numel(csvdata(csvdata>=median1));%count median and above pixels
    countotsu=numel(csvdata(csvdata>=thresh));%count otsu's and above pixels
    datatable(i,1)=count;%
    datatable(i,2)=average;% median threshold average
    datatable(i,3)=median1;%median of all nonnegative values
    datatable(i,4)=rsd;% 
    datatable(i,5)=thresh;% otsu's threhold value
    datatable(i,6)=mean(reshape(csvdata(threshmean),1,[]));% otsu's threshold average
    datatable(i,7)=rawmean; % raw positive average
    datatable(i,8)=threshmedian; % otsu median
    datatable(i,9)=threshthreshmean1; % otsu 2 mean
    datatable(i,10)=thresh2; % otsu 2 threshold
    datatable(i,11)=rawmode; % mode
    datatable(i,12)=rawmax; % max
    datatable(i,13)=countall; % count all pixels
    datatable(i,14)=countmedian; % count all pixels equal or above median
    datatable(i,15)=countotsu; % count all pixels equal or above otsu
    header = {'Count',' median threshold average','median of all non-negative values','RSD','Otsu threshold','Otsu threshold Average','Raw Average','Otsu median','Otsu 2 average','Otsu 2 threshold','mode','max','count all','count median','count otsu'};
    output = [header; num2cell(datatable)];
    seg_I=imquantize(csvdata,thresh);
    seg_I=mat2gray(seg_I);
    imwrite(seg_I, sprintf('%s%sTHRESH.png',dname,filename));
end

