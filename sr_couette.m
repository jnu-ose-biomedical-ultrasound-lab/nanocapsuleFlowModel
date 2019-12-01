function [v_zhijie v_max SR_r]=sr_couette(T_Psan,tt,yy, lx,Bound_y)
% Bound_y=0.1e-3;      %dan wei (mm)  0.1mm
% Bound_x=0.001;        %1mm
% lx=1000;
% xx=(1:lx)/lx*Bound_x;
% ly=200;
% yy=(-ly/2:ly/2)/ly*Bound_y;
%%%%%%%%%%%%%%%%%%
%%
% T_P=200;                                                              %100s  time
% T_Psan=100;                                                         %period 100s
% tt_z=0:5:T_P;

%    sr=(cos(srW*tt)+1)*0.5;
t=2*pi/T_Psan*(tt);
sanjiao(tt)=(sawtooth(t,0.5)+1)/2;
SR_r=sanjiao(tt);
% SR_r=1;
v_zhijie=(SR_r*(yy+Bound_y/2))'*ones(1,lx);
v_max=SR_r*Bound_y;
