function n=n_RBC(x,y,r_RBC,con)%长度，宽度，浓度
%计算RBC个数
s_RBC=x*y*con;
% r_RBC=2.7e-6;
s_RBC1=pi*r_RBC^2;
n=s_RBC/s_RBC1;