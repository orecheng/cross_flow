clc
clear all
filename='³ýÊª.xlsx';
sheet=1;
data_raw=xlsread(filename,sheet);
[n,~]=size(data_raw);
%% Data Classification
SolutionMassFraction=data_raw(:,3);
SolutionTemp=data_raw(:,4);
SolutionVol=data_raw(:,5);

AirVol=data_raw(:,6);
AirInTemp=data_raw(:,8);
AirInRH=data_raw(:,9);
AirInMoi=data_raw(:,10);
AirOutTemp=data_raw(:,11);
AirOutRH=data_raw(:,12);
AirOutMoi=data_raw(:,13);
%% Design of NTU eqution
%NTU=a*V*Fa^b*Fz^c*(1-ksi)^d
function NTU=NTU_regression_cal(V,Fa,Fz,ksi,a,b,c,d)
NTU=a*V*Fa^b*Fz^c*(1-ksi)^d;
end