function [rho_air,da,ha]= rh2da(t,phi)
% t=41.1;
% phi=0.649;
if phi>10
phi=phi/100;
end
T=t+273.15;
Pma=101325;
Rda=287;
ps=exp(7.23e-7*t^3-2.71e-4*t^2+7.2e-2*t+6.42);
da=0.622*phi*ps/(Pma-phi*ps);
rho_air=Pma*(1+0.001*da)/(Rda*T*(1+0.001606*da));% 湿空气密度 （kg/kg干空气）
ha=1.01*t+(1.83*t+2501)*da;
end


