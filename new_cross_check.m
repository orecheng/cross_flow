function [AveAirTempratureOut,AveAirMoistureOut]=new_cross_check(Ta_in,phi,Ts_in,Ps_in,Va_in,Vs_in,H,L,NTU)
% Ta_in=30;%
% phi=0.8;%

T=Ta_in+273.15;
[rho_air,da_in,ha_in]= rh2da(Ta_in,phi);
% Ts_in=15;%
% Ps_in=0.3;%
rho_licl=cal_rho_licl(Ts_in,Ps_in);

% Va_in=4000;%
% Vs_in=4.42;%

Ms_in=Vs_in*rho_licl/3600;
Ma_in=Va_in*rho_air/3600;

% H=1;%
% L=0.15;%

hs_in=enthalpy(Ts_in,Ps_in);

meshgrid=0.005;
M=ceil(L/meshgrid);
N=ceil(H/meshgrid);

Ms_divide=Ms_in;
Ma_divide=Ma_in;

% NTU=5;%
%% 水在空气中的扩散系数
Dair=(-0.29890+1.6253e-3*T+7.5e-7*T^2)*1e-4;
%% 空气的动力粘度和运动粘度
mu_da=[-2.448e-6,0.005072,1.713]*[Ta_in^2;Ta_in;1]*1e-5; %动力粘度
nv_da=mu_da/rho_air; %运动粘度
%% Thermal Conductivity of Air
lamda=[6.993e-8,0.007618,2.442]*[Ta_in^2;Ta_in;1]*1e-2;
%% Sc
Sc=nv_da/Dair;
%% Pr
Cpa=1.01e3;
Pr=mu_da*Cpa/lamda;
Le=Sc/Pr;
%% Initial of Matrix
SolEnthalpy=zeros(N,M);
AirMoisture=zeros(N,M);
AirEnthalpy=zeros(N,M);
AirTemprature=zeros(N,M);
SolTemprature=zeros(N,M);
SolMass=zeros(N,M);
SolConcen=zeros(N,M);
%% row_1 column_1
[AirMoisture(1,1),AirTemprature(1,1),AirEnthalpy(1,1),SolTemprature(1,1),SolMass(1,1),SolConcen(1,1),SolEnthalpy(1,1)]= cross_core(Ts_in,Ps_in,Ms_divide,ha_in,da_in,Ma_divide,M,N,NTU,Le);
%% row_1
for j=2:M
    [AirMoisture(1,j),AirTemprature(1,j),AirEnthalpy(1,j),SolTemprature(1,j),SolMass(1,j),SolConcen(1,j),SolEnthalpy(1,j)]= cross_core(Ts_in,Ps_in,Ms_divide,AirEnthalpy(1,j-1),AirMoisture(1,j-1),Ma_divide,M,N,NTU,Le);
end
%% column_1
for i=2:N
    [AirMoisture(i,1),AirTemprature(i,1),AirEnthalpy(i,1),SolTemprature(i,1),SolMass(i,1),SolConcen(i,1),SolEnthalpy(i,1)]= cross_core(SolTemprature(i-1,1),SolConcen(i-1,1),SolMass(i-1,1),ha_in,da_in,Ma_divide,M,N,NTU,Le);
end
%% row_2~row_end  column_1~column_end
for i=2:N
    for j=2:M
        [AirMoisture(i,j),AirTemprature(i,j),AirEnthalpy(i,j),SolTemprature(i,j),SolMass(i,j),SolConcen(i,j),SolEnthalpy(i,j)]= cross_core(SolTemprature(i-1,j),SolConcen(i-1,j),SolMass(i-1,j),AirEnthalpy(i,j-1),AirMoisture(i,j-1),Ma_divide,M,N,NTU,Le);
    end
end

%% Average Out
AveAirMoistureOut =mean(AirMoisture(1:N,M));
AveAirEnthalpyOut =mean(AirEnthalpy(1:N,M));
AveAirTempratureOut =mean(AirTemprature(1:N,M));

AveSolTempratureOut =mean(SolTemprature(N,1:M));
AveSolMassOut=mean(SolMass(N,1:M));
AveSolConcenOut=mean(SolConcen(N,1:M));
AveSolEnthalpyOut=mean(SolEnthalpy(N,1:M));

err=(ha_in-AveAirEnthalpyOut)*Ma_in/(sum(SolEnthalpy(N,1:M).*SolMass(N,1:M)/M)-hs_in*Ms_in)-1;
backindex=AveAirTempratureOut+AveAirMoistureOut*1000;
% backindex=AveAirTempratureOut;