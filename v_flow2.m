function [vmid T]=v_flow2(x,y,t,dv,Vp_flow,SR,v_mean)
%流速子函数
%dv波动幅度
% vmean=0.05;
% vmean=0.02;
% dv=0.01;
% w_flow=8;
% k_flow=0.05;
% SR=2e-4;
% SR=10;
f_flow=SR/60;                                                      %流速的频率f
% Vp_flow=0.001;                                                     %???
% Vp_flow=0.3;
w_flow=f_flow*(2*pi);                                          %流速的角频率 w
k_flow=w_flow/Vp_flow;                                      %wavenumber=w/v
T=1/(k_flow/2/pi); %T相当于在x轴的长度。
D=0.1*10^-3;
% tmp(:,k2) =-1+(k2-1)*0.01; 
% vmid=vmean+dv*sin(1*w_flow*t-k_flow*x+pi/2);
% Vf_M=dv*sin(w_flow*t*T-k_flow*P_R(:,1)+in_t)+v_mean; 
vmid=(dv*sin(1*w_flow*t-3000*x+pi/2)+v_mean).*(1-(4.*(y.^2)./(D^2)));
% Bound_y=10e-5;
Bound_y=0.1e-3; 
% u=3e-3;
u=vmid;
Vs_R= 2.7000e-006;
p=vmid/Bound_y^2*8*Vs_R;
% uu=1/2/u*p*y*(Bound_y-y);
if y<0
    uu=1/2/Vs_R*p*(Bound_y/2-abs(y))*(Bound_y-(Bound_y/2-abs(y)));
else
    uu=1/2/Vs_R*p*(y+Bound_y/2)*(Bound_y-(y+Bound_y/2));
end


% v=vmid-SR*abs(y);                                         %shear rate
% v=vmid;
