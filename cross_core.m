function [da_out,Ta_out,ha_out,Ts_out,Ms_out,Ps_out,hs_out]= cross_core(Ts_in,Ps_in,Ms_in,ha_in,da_in,Ma_in,M,N,NTU,Le)

A=2-(1+(Ps_in/0.28)^4.3)^0.6;
B=(1+(Ps_in/0.21)^5.1)^0.49-1;
Fpt=A+B*(((Ts_in+273.15)/228)-1);
Pi25=1-(1+(Ps_in/0.362)^(-4.75))^(-0.4)-0.03*exp(-(((Ps_in-0.1)^2)/0.005));
t=1-(Ts_in+273.15)/647.14;
E=exp(((-7.85823)*t+1.83991*(t^1.5)+(-11.7811)*(t^3)+22.6705*(t^3.5)+(-15.9393)*(t^4)+1.77516*(t^7.5))/(1-t));
Pw=2.212*10^7*E;
Psol=Pw*Pi25*Fpt;                                %溶液表面的蒸汽压
de_in=0.622*(Psol/(101325-Psol));    %与溶液状态相平衡的饱和空气含湿量
%% 计算溶液焓值
X=Ps_in*100;
a=[-66.2324,11.2711,-0.79853,2.1534e-2,-1.66352e-4];
b=[4.5751,-0.146924,6.307226e-3,-1.38054e-4,1.06690e-6];
c=[-8.09689e-4,2.18145e-4,-1.36194e-5,3.20998e-7,-2.64266e-9];
x=[1;X;X^2;X^3;X^4];
A=a*x;
B=b*x;
C=c*x;
hs_in=A+B*Ts_in+C*Ts_in^2;
%% he,i-1 换算 Tsin de,i-1 
he_in=1.01*Ts_in+(2501+1.83*Ts_in)*de_in;
%% da,i ①
da_out=da_in+NTU*(de_in-da_in)/M;
%% ha,i ②
rts_in=2500-2.35*Ts_in;
ha_out=ha_in+NTU*Le*(he_in-ha_in+rts_in*(1/Le-1)*(de_in-da_in))/M;
Ta_out=(ha_out -2501*da_out)/(1.01+1.83*da_out);
%% Ms,j ③
Ms_out=Ms_in+Ma_in*(da_in-da_out)*M/N;
%% Ps,j ④
Ps_out=Ps_in*Ms_in/Ms_out;
%% hs,j ⑤
hs_out=(hs_in*Ms_in-Ma_in*(ha_out-ha_in)*M/N)/Ms_out;
%% Ts,j
X=Ps_out*100;
a=[-66.2324,11.2711,-0.79853,2.1534e-2,-1.66352e-4];
b=[4.5751,-0.146924,6.307226e-3,-1.38054e-4,1.06690e-6];
c=[-8.09689e-4,2.18145e-4,-1.36194e-5,3.20998e-7,-2.64266e-9];
x=[1;X;X^2;X^3;X^4];
A=a*x;
B=b*x;
C=c*x;
Ts_out=(-B+sqrt(B^2-4*C*(A-hs_out)))/(2*C);
