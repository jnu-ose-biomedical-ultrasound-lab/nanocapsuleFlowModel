function lo=dist_RBC_tube_bound(x,y,n,R,A_elastic,B_elastic,k_e,f_e,T_num);
%% input data 
% x = Bound_x-2*R;
% y = Bound_y-2*R;
% n = PNUM;
% R = R;
% T_num = T_num(1);
%RBC无重叠分布
lo=rand(n,1);
% x=x-4*R;
lo(:,1)=lo(:,1)*x;

for i = 1 : n
vessel_elasticity_tube(i) = A_elastic*cos(2*pi*f_e*T_num-k_e*lo(i,1))+(B_elastic-A_elastic);
a = vessel_elasticity_tube(i)-2*R;
b = -vessel_elasticity_tube(i)+2*R;
lo(i,2) = a + (b-a).*rand(1,1);
end
r2=(R*2)^2*1.1;
for k = 1:n
d=(lo(k,1)-lo(:,1)).^2+(lo(k,2)-lo(:,2)).^2;
    d(k)=inf;
    while min(d)<r2 

        lo(k,1) = rand()*x;
        yy = A_elastic*cos(2*pi*f_e*T_num-k_e*lo(k,1))+(B_elastic-A_elastic);
        a = yy-2*R;
        b = -yy+2*R;
        lo(k,2) =  a+(b-a).*rand();
        d=(lo(k,1)-lo(:,1)).^2+(lo(k,2)-lo(:,2)).^2;
        d(k)=inf;
    end
end
% 
% figure;
% plot(lo(:,1),lo(:,2),'ro')
% hold on 
% hold on; plot(P_R_x,A_elastic*cos(-k_e*P_R_x)+(B_elastic-A_elastic))
% hold on; plot(P_R_x,-(A_elastic*cos(-k_e*P_R_x)+(B_elastic-A_elastic)))

end