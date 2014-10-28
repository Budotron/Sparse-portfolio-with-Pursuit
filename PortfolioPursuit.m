% clc; clear all;
% load datapoints
clc; clear all;
warning ( 'off' , 'all' );
warning;
files=dir('/Users/akikotakeda/Dropbox/code/daily');
num_files = length(files);
datapoints=[];
indexset=[];
names=[];
startdates=[];
endDates=[];
datesSet=[]
for i=4:num_files
     data=csvread(files(i).name);
     [m,n]=size(data);
     
     if m>=3000
         dates=data(m-2999:m, 1);
         datesSet=[datesSet dates];
         string={files(i).name};
         cut1={string{1}(7:end)};
         cut2={cut1{1}(1:end-4)};
         names=[names, cut2];
         indexset=[indexset, i];
         startdate=data(m-2999, 1);
         enddate=data(m, 1);
         data=data(m-2999:m, 6);
         datapoints=[datapoints, data];
     end
end
% DATASETS = MAT2CELL(DATAPOINTS,300*ONES(10,1),440);
% 
% DATESSETS = MAT2CELL(DATES,300*ONES(10,1),1);
