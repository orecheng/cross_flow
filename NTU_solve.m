clc
clear all
filename='0904¼ò»¯Êý¾Ý.xlsx';
sheet=1;
data_raw=xlsread(filename,sheet);
[n,~]=size(data_raw);

raw_ps=data_raw(:,3);
raw_ts=data_raw(:,4);
raw_vs=data_raw(:,5);
raw_va=data_raw(:,6);
raw_ta_in=data_raw(:,8);
raw_phi=data_raw(:,9)/100;
raw_ta_out=data_raw(:,11);
raw_da_out=data_raw(:,13)/1000;

err=+0.05;
for i=1:n
    Ta_in=raw_ta_in(i)*(1+err);
    phi=raw_phi(i)*(1+err);
    Ts_in=raw_ts(i);
    Ps_in=raw_ps(i);
    Va_in=raw_va(i);
    Vs_in=raw_vs(i);
    ta_out=raw_ta_out(i)*(1-err);
    da_out=raw_da_out(i)*(1-err);
    
    H=0.5;
    L=0.15;
    W=0.65;
    V=H*L*W;
    aw=380;
    
    
    [rho_air,da_in,ha_in]= rh2da(Ta_in,phi);
    rho_licl=cal_rho_licl(Ts_in,Ps_in);
    Ms_in=Vs_in*rho_licl/3600;
    Ma_in=Va_in*rho_air/3600;
    
    NTU_up=10;
    NTU_down=-10;
    
    opt = optimset('display','off','TolFun',1e-12,'TolX',1e-12);
    backindex=ta_out+da_out*1000;
% backindex=ta_out;
    [NTU(i,:),resnorm,residual,exitflag,output] =lsqnonlin(@(NTU)new_cross_cal(Ta_in,phi,Ts_in,Ps_in,Va_in,Vs_in,H,L,NTU)-backindex,2,NTU_down,NTU_up,opt);
% [NTU(i,:),resnorm,residual,exitflag,output] =fsolve(@(NTU)new_cross_cal(Ta_in,phi,Ts_in,Ps_in,Va_in,Vs_in,H,L,NTU)-backindex,5,opt);
    hd(i,:)=NTU(i,:)*Ma_in/(aw*V);
    [ta_out_cal(i,:),da_out_cal(i,:)]=new_cross_check(Ta_in,phi,Ts_in,Ps_in,Va_in,Vs_in,H,L,NTU(i));
end



