function n=n_RBC(x,y,r_RBC,con)%���ȣ���ȣ�Ũ��
%����RBC����
s_RBC=x*y*con;
% r_RBC=2.7e-6;
s_RBC1=pi*r_RBC^2;
n=s_RBC/s_RBC1;