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

NTU=data_raw(:,17);
hd=data_raw(:,18);

H=0.5;W=0.65;L=0.15;
V=H*L*W;
%% Pre Calculation
%Density of Solution and Air
ksi=SolutionMassFraction;
for i=1:n
SolutionRho(i,:)=cal_rho_licl(SolutionTemp(i),ksi(i));
[AirInRho(i,:),~,~]= rh2da(AirInTemp(i),AirInRH(i));
end

AirInMass=AirInRho.*AirVol;
SolutionInMass=SolutionVol.*SolutionRho;
Fa=AirInMass/(L*H*3600);
Fz=SolutionInMass/(L*W*3600);
%% fit
