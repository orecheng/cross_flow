function rho_licl=cal_rho_licl(T_LiCl,ksi)
%% Input
% T_LiCl=18 ;%Temperature of solution
% ksi=0.2472;%mass fraction of LiCl

%% Initializaion
T=T_LiCl+273.15;
Tc=647.226;
Pc=22090e3;
theta=T/Tc;
%% 3.cal of LiCl density
B0=1.9937718430;
B1=1.0985211604;
B2=-0.5094492996;
B3=-1.7619124270;
B4=-45.9005480267;
B5=-723692.2618632;

rho_c_h2o=322;
tao=1-theta;
rho_h2o=rho_c_h2o*(1+B0*tao^(1/3)+B1*tao^(2/3)+B2*tao^(5/3)+B3*tao^(16/3)+B4*tao^(43/3)+B5*tao^(110/3));

rho0=1;
rho1=0.540966;
rho2=-0.303792;
rho3=0.100791;

ksical=ksi/(1-ksi);
rho_licl=rho_h2o*(rho0+rho1*ksical+rho2*ksical^2+rho3*ksical^3);%density of aqueous LiCL