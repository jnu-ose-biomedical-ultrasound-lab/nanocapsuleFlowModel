%% Force f_e-f_a
function F=GetMF(d_R, AGG,R)

% m=9.8*10^(-14);                                                 %mass of RBC
% R=2.7*10^(-6); %m                                             %radius of RBC
% R=4*10^(-6); %m                                             %radius of RBC
k=3*10^(-6); %N/m                                            %elastic modulus (RBC is not rigid)
da_0=11*10^(-9); %(m)                                      %Zero-force length(nm) tethe0
% da=10^(-24);                                                      %DA,Energy(J)
da=10^(-24);
% da=10^(-18);
b=10^7;%(1/m)                                                   %B, scaling factor(m-1)
da_d=d_R-2.0*R;                                                  %����RBC�ľ����ȥ�����뾶������RBC��������

switch AGG
    case 1
        %AGG-
        DA=da/10;%AGG-
        B=b;%(1/m)
    case 2
        %AGG+
        DA=da;
        B=b;%(1/m)
    case 3
        %AGG+1/2B
        DA=2*da;
        B=b/2;%(1/m)
        case 4
        %AGG+1/2B
        DA=3*da;
        B=b/3;%(1/m)
end

log_fe=d_R<(2.*R);                                              %�������RBC�ľ���С��ֱ���ĳ��ȣ�˵������RBC�Ѿ�������
f_e=k*(2*R-d_R).^(3/2).*log_fe;                           %���ֵ����� F_e

log_fa=(da_d>=da_0);                                      %����RBC����ľ����������20����RBCѪ�ܱں�Ⱦ���Ϊû����������
    f_a=(2*DA*B*(exp(2*B*(da_0-da_d))-exp(B*(da_0-da_d)))).*log_fa;   %Aggregation Force


                                    %����RBC����ľ����������20����RBCѪ�ܱں�Ⱦ���Ϊû����������
%  f_a=2*DA*B*(exp(2*B*(da_0-da_d))-exp(B*(da_0-da_d)));   %Aggregation Force
%  f_a((da_d -da_0)<=0 )=0;

% f_a=2*DA*B*(exp(2*B*(da_0-da_d))-exp(B*(da_0-da_d))); 
% if (da_d-da_0) <= 0
%     F=f_e+f_a;                                                             %��Ϊ����RBC���������ͳ����Ƿ����
% else
    F=f_e+f_a;
% end
% F=f_a;
% F=f_e;
