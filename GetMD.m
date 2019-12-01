% Get the distance and direction among particles������RBC�ľ���ͷ���
% Rnm is the matrix of distances among particles
%Unit_V is the unit direction vector among particles
%Arrange like
%Rnm��d_R��:  ��ʾ���룺��һ�д����һ���㵽������ľ���(����������ĵ�)���ڶ��д���ڶ����㵽������ľ���
%Unit_V: ��ʾ����ĸ������֣�x+yi��x�Ǻ����꣬y�������꣬abs��Unit_V��=Rnm
%d_RF:  ��ʾ���룺��һ�д����һ���㵽������ľ���(��������ĵ�)���ڶ��д���ڶ����㵽������ľ���
              %�м�ĶԽ��߱�������㵽�㱾��ľ��루0��������������100,

function [Rnm Unit_V d_RF]=GetMD(P_R,x_boundy,R)
P_n=size(P_R,1);                                                  %���ٸ�RBC��
One_Pn=ones([1,P_n]);                                        %���ǵ�λ1��һ��
One_Pn_1=ones([1,P_n-1]);                                 %RBC������һ

%xn_xm
P_x=P_R(:,1);                                                        %��һ�У�11����ĺ�����
P_xmn=(P_x*One_Pn);                                          %ÿһ�ж���ͬ�ķ��󣬺����궼�̿���
P_xmn_1=(P_x*One_Pn_1)';                                  %ÿһ�ж���ͬ��ʮ�У�������һ�У�
b=triu(P_xmn,1);                                                  %�ӵڶ��е����һ��
b=b(1:end-1,:);                                                     %һ�����Ǵӵ�һ�������һ����ĺ����꣬�ӵ�һ���㵽���е�
c=tril(P_xmn,-1);
c=c(2:end,:);                                                         % ��һ�������Ե�ĺ����굹������ �����е㵽��һ����
P_xx=b+c;                                                            %��ֵ����
P_xn_xm=P_xmn_1-P_xx;                                      %��һ�д����һ���㵽������ľ��룬�ڶ��д���ڶ����㵽������ľ���
                                                                             %�Դ�����(��ֻ�Ǻ�����)
%yn_ym
P_y=P_R(:,2);
P_ymn=(P_y*One_Pn);
P_ymn_1=(P_y*One_Pn_1)';
b=triu(P_ymn,1);
b=b(1:end-1,:);
c=tril(P_ymn,-1);
c=c(2:end,:);
P_yy=b+c;
P_yn_ym=P_ymn_1-P_yy;
P_xn_xm(P_xn_xm>x_boundy-4*R)=P_xn_xm(P_xn_xm>x_boundy-4*R)-x_boundy;
P_xn_xm(P_xn_xm<-(x_boundy-4*R))=P_xn_xm(P_xn_xm<-(x_boundy-4*R))+x_boundy;

%Rnm
Rnm=sqrt(P_xn_xm.^2+P_yn_ym.^2);                %��ʾ���룺��һ�д����һ���㵽������ľ��룬�ڶ��д���ڶ����㵽������ľ���
%Angle
Unit_V=(P_xn_xm+P_yn_ym*i)./Rnm;

%d_RF full size D to D
P_xf=P_xmn-P_xmn';                                           %��ʾ�����꣺��һ�д����һ���㵽���е�ľ��루��������ĵ㣩
P_xf(P_xf>x_boundy-4*R)=P_xf(P_xf>x_boundy-4*R)-x_boundy;
P_xf(P_xf<-(x_boundy-4*R))=P_xf(P_xf<-(x_boundy-4*R))+x_boundy;
    %�ڶ��д���ڶ����㵽���е�ľ��루��������ĵ㣩
P_yf=P_ymn-P_ymn';
d_RF=sqrt(P_xf.^2+P_yf.^2);                               %��ʾ�ܵľ���
d_RF=d_RF+100*diag(One_Pn);                          %�Խ���Ϊʲô��100������
