clear
close all
clc

%% thresholding

dname='D:\Uni\MIAK\wanagat via dave\multiplex\0Processing\multiplex\imageJ\thresholds\';
files=dir([dname,'*.csv']);
for i=1:length(files);
    csvdata=importdata([dname,files(i).name]);
    csvdata(csvdata<0)=0; %remove all negative values
    count=i;
    rawmean=mean2(csvdata);
    median=nanmedian(nanmedian(csvdata,1),2);
    csvdata1=csvdata-median; %remove median (aka background) signal
    csvdata1(csvdata1<0)=0; %remove all negative values
    average=mean2(csvdata1(csvdata1~=0));  %average of values greater than zero
    rsd=(std2(csvdata1)/mean2(csvdata1))*100;
    thresh = multithresh(csvdata); %otsu's method thresholding
    [cluster_idx, cluster_center] = kmeans(csvdata,2,'distance','sqEuclidean'); %2 level k mean clustering
    threshmean=csvdata-thresh;
    %boxplot(csvdata.data)
    datatable(i,1)=count;
    datatable(i,2)=average;% median threshold average
    datatable(i,3)=median;%median of all nonnegative values
    datatable(i,4)=rsd;% 
    datatable(i,5)=thresh;% otsu's threhold value
    datatable(i,6)=-mean2(threshmean(threshmean~=0));% otsu's threshold average
    datatable(i,7)=rawmean; % positive integer average
    %datatable(i,8)=C; %eventually kmean
end
