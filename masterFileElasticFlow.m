clf; clear all; clc;
%cd pwd

%--------------------------------------------------------------------------
 % masterFileElasticFlow.m

 % Last updated: December 2019, John LaRocco
 
 % Jeju National University-Biomedical Ultrasound Lab
 
 % Details: Master file that simulates nanocapsules in Womersley flow.
 % Adapted from red blood cell model by Cheong Ah Lee: 
 %  https://github.com/jnu-ose-biomedical-ultrasound-lab/bloodflowmodeling
 
 % Additional details in papers:
 % LEE, Cheong-Ah; KONG, Qi; PAENG, Dong-Guk. 
 % Depletion-model-based numerical simulation of the kinetics of red blood cell aggregation under sinusoidal pulsatile flow. 
 % Biorheology, 2018, Preprint: 1-13.
 
 % Inputs:
 % AGG: Adhesion and aggregation mode [Scalar from 1-4, in increasing order]
%--------------------------------------------------------------------------



%% 1. Tube Design
T =1:1e2:40e3+1;
t= [0.1e-3];
T_num = (T-1)*t;
r=0.000055; % tube radius in m
grid_px=300; timestep=length(T_num); 
Bound_y = 2*r; % tube height in m, diameter
P_R_y = linspace(-Bound_y/2,Bound_y/2,grid_px)';
Bound_x = 0.0005; %tube length in m
P_R_x =  linspace(0,Bound_x,500);
ru=1060; %kg/m3 (rbc) 
ru=2000; %kg/m3 (pfc) 
freq=1.0; % 60 bpm
mu= 2.7e-6; %pascal*s 
omega=2*pi*freq;
alpha=r*sqrt(omega*ru/mu); 
%% 2. wall motion
A_elastic = Bound_y/22; 
    %A_elastic = 0; % rectangular boundary
f_e = freq;
lamda = Bound_x/4;
k_e = 2*pi/lamda;
B_elastic = Bound_y/2;
Bound_vessel = linspace(-(A_elastic+B_elastic),(A_elastic+B_elastic),grid_px)';
rr = A_elastic+B_elastic;
Bound_y = 2*rr;
P_R_y = linspace(-Bound_y/2,Bound_y/2,grid_px)';

%% 3. Particle Properties;
mass=9.8*10^(-14); %kg of rbc 
mass=1.0471975512*10^(-12); %kg of pfd nanocapsule
% volume: 5.23598775598E-16
R=4*10^(-6); % radius in m for RBC
R=0.000005; % radius of PFD nanocapsule

AGG=1; % aggregation mode
PNUM= floor(n_RBC(Bound_x,Bound_y,R,0.4));
V_R=zeros([PNUM 2]);     %particle velocity field                                      
lo=rand(PNUM,1);
lo(:,1)=lo(:,1)*Bound_x;
for m = 1 : PNUM
vessel_elasticity_tube(m) = A_elastic*cos(2*pi*f_e*T_num(1)-k_e*lo(m,1))+(B_elastic);
a = vessel_elasticity_tube(m)-2*R;
b = -vessel_elasticity_tube(m)+2*R;
lo(m,2) = a + (b-a).*rand(1,1);
end
r2=(R*2)^2*1.1;
for k = 1:PNUM
d=(lo(k,1)-lo(:,1)).^2+(lo(k,2)-lo(:,2)).^2;
    d(k)=inf;
    while min(d)<r2 

        lo(k,1) = rand()*Bound_x;
        yy = A_elastic*cos(2*pi*f_e*T_num(1)-k_e*lo(k,1))+(B_elastic);
        a = yy-2*R;
        b = -yy+2*R;
        lo(k,2) =  a+(b-a).*rand();
        d=(lo(k,1)-lo(:,1)).^2+(lo(k,2)-lo(:,2)).^2;
        d(k)=inf;
    end
end
% Bound = 2*0.000045;
% lo=dist_RBC(Bound_x-2*R,Bound_y-2*R,PNUM,R);
 P_R=[lo(:,1)+R lo(:,2)];
 
 %% 4. Simulation write video
% name output file
outputName=['nanocapsuleElasticWomersleyFlowModel' num2str(AGG) '.avi'];

aviobj=VideoWriter(outputName);
aviobj.FrameRate = 1;                                                        
open(aviobj);
%% 5. RBC Simulation
uu_f = zeros(size(P_R_y,1),size(P_R_x,2));
u = zeros(size(P_R_y,1),size(P_R_x,2),timestep ); 
shear_field = zeros(size(P_R_y,1),size(P_R_x,2),timestep ); 
acc_field = zeros(size(P_R_y,1),size(P_R_x,2),timestep ); 
once_t = 0;
pn(1) = 30;
kkk = 1;
for ii = 2: timestep
    vessel_elasticity_pre = A_elastic*cos(2*pi*f_e*T_num(ii-1)-k_e*P_R_x)+(B_elastic);
    vessel_elasticity_up = A_elastic*cos(2*pi*f_e*T_num(ii)-k_e*P_R_x)+(B_elastic);
    vessel_elasticity_down = -(A_elastic*cos(2*pi*f_e*T_num(ii)-k_e*P_R_x)+(B_elastic));
    temp_num = find(vessel_elasticity_up == max(vessel_elasticity_up));  
    Area_front = pi*(vessel_elasticity_up(temp_num(1))^2);


%% Flow field %%
    
    for j = 1: size(vessel_elasticity_up,2) 
        Area_back = pi*(vessel_elasticity_up(j)^2);   
        P_R_y_tube = find(P_R_y >=vessel_elasticity_down(j)& P_R_y <=vessel_elasticity_up(j));
        P_R_y_new = linspace(-Bound_y/2,Bound_y/2,length(P_R_y_tube))';
        u_f = zeros(size(P_R_y_new,1),1);
        shear = zeros(size(P_R_y_new,1),1);
        acc = zeros(size(P_R_y_new,1),1);
        if ii == 2
            u_ff = zeros(size(P_R_y_new,1),timestep);
           for iii = 1:1: 101   
               for l = 1
                   kapa=(alpha*(l^(0.5))*i^(1.5))/rr;
                   CJA = besselj(0,kapa*rr);
                   CBJ0=besselj(0,kapa*P_R_y_new);                
                   u_ff(:,iii)=u_ff(:,iii)+real(((i*pn(l))/(ru*omega*l))*(1- CBJ0/CJA)*cos((omega*l*T_num(iii))));                 
               end 
           end
        min_u = min(min(u_ff));   
        p0 = (min_u*4*mu/(rr^2))*13;
        end
        
               for l = 1
                   kapa=(alpha*(l^(0.5))*i^(1.5))/rr;
                   CJA = besselj(0,kapa*rr);
                   CBJ0=besselj(0,kapa*P_R_y_new);  
                   CBJ1=besselj(1,kapa*P_R_y_new);
                   u_f(:,1)=u_f(:,1)+real(((i*pn(l))/(ru*omega*l))*(1- CBJ0/CJA)*cos((omega*l*T_num(ii))));  
                   shear(:,1)=shear(:,1)+abs(((i*pn(l))/(ru*omega*l))*(-(kapa/rr)*(CBJ1/CJA))*cos((omega*l*T_num(ii)))); 
                   acc(:,1)=acc(:,1)+real(((i*pn(l))/(ru*omega*l))*(1- CBJ0/CJA)*(omega*l)*-1*sin((omega*l*T_num(ii))));
               end 

            u_f(:,1)=u_f(:,1)+(p0*((P_R_y_new).^2-rr^2)/4/mu);
            shear(:,1) = shear(:,1)+(-(p0*2*P_R_y_new)/(4*mu));
            
        u(P_R_y_tube,j,ii) = u(P_R_y_tube,j,ii)+u_f(:,1)*(Area_front/Area_back); % convert to flow rate  
        shear_field(P_R_y_tube,j,ii)  =  shear_field(P_R_y_tube,j,ii) +shear(:,1)*(Area_front/Area_back);   
        acc_field(P_R_y_tube,j,ii)  =  acc_field(P_R_y_tube,j,ii) +acc(:,1)*(Area_front/Area_back);   
    end
  
      
    for k = 1: size(P_R,1)
        num1 = find(P_R_x(:)<=P_R(k,1));
        num2 = find(P_R_y<= P_R(k,2));
        uu(k,1) = u(num2(end),num1(end),ii);
        SR(k,1) = shear_field(num2(end),num1(end),ii);
        accel(k,1) = acc_field(num2(end),num1(end),ii);
    end
  
%%         
      
    [d_R FA_un d_RF]=GetMD(P_R,Bound_x,R);   %Get the distance and direction among particles    
    FA_R = GetMF(d_R,AGG,R).*FA_un;    
    Fh_R=6*pi*mu*R*([uu,zeros(PNUM,1)]-V_R);               % Hydrodynamic force
    F_Rt=[real(sum(FA_R, 1))' imag(sum(FA_R, 1))'];          % aggregation force and elastic force
    F_R=[real(sum(FA_R, 1))' imag(sum(FA_R, 1))'] +Fh_R;     % total force of particle
    AC_R=F_R./mass;                                             % accelaration (Newton's second law)
    Step_xy=V_R.*t+0.5*AC_R.*T_num(ii)^2;                            % distance                        
    Step_D=sqrt(sum((Step_xy.^2),2));            
    dist_p = max(Step_D);
    dist_m = mean(Step_D);
    V_R=V_R+AC_R.*t;                                         % velocity function [m/s]                     
    P_R_pre = P_R;
    P_R=P_R+V_R.*t+0.5*AC_R.*t^2;                            % particle position [m]
    displace = V_R.*t+0.5*AC_R.*t^2;                         % particle distance [m]    
    once_t=once_t+t  ;
    Ovx_index1=find(P_R(:,1)>Bound_x);
    P_R(Ovx_index1,1)=P_R(Ovx_index1,1)-(Bound_x);
    Ovx_index2=find(P_R(:,1)<=0);
    P_R(Ovx_index2,1)=P_R(Ovx_index2,1)+(Bound_x);

    for jj = 1: length(vessel_elasticity_up)-1
        numx = find(P_R(:,1)>=P_R_x(jj)&P_R(:,1)<=P_R_x(jj+1));
        del_tube(jj) = vessel_elasticity_pre(jj)-vessel_elasticity_up(jj);
        slop = 1-12000*del_tube(jj); 
        P_R(numx,2) = slop*P_R(numx,2);     
    end

     if  once_t> t*4
       % ii
        AZ=d_RF<2.04*R;                                      
        AZ1=AZ;
        AZ=triu(AZ);                                                
        [AGG_x,AGG_y]=find(AZ==1);                   
        AGG_p=[AGG_x AGG_y];   
        [Agg_size Agg_no Agg aa RBC_22]=getSN_Cybersys(AGG_p);
        AGG_RBC_No=unique(AGG_p); 
        %% %% Cybersys begin
%      text(P_R(AGG_RBC_No,1)-R,P_R(AGG_RBC_No,2),num2str(AGG_RBC_No),'FontSize',15);
        window_x=[0.2 0.3]*1e-3;    
        win_grid = find(P_R_x>=window_x(1)&P_R_x<=window_x(2));
        P_R_win=P_R(AGG_RBC_No,1)>window_x(1,1) & P_R(AGG_RBC_No,1)<window_x(1,2);

        if length(AGG_RBC_No)>0
            if size(AGG_RBC_No,2)>1
                AGG_RBC_No2=AGG_RBC_No';
            else
                AGG_RBC_No2=AGG_RBC_No;
            end
            P_R_win2=P_R_win.*AGG_RBC_No2;
            P_R_win2(P_R_win2==0)=[]; 
            num_agg = numel(P_R_win2);
%             DrawCircle(P_R(P_R_win2,1),P_R(P_R_win2,2),R/1, 100, 'b.');
%             text(P_R(P_R_win2,1)-R,P_R(P_R_win2,2),num2str(P_R_win2),'FontSize',15);
%             hold off          
            Agg_s2(kkk,1)=num_agg;
            Agg_s2(kkk,2) = mean(mean(u(:,win_grid,ii )));
            Agg_s2(kkk,3) = numel(AGG_RBC_No);
            net1=[find(Agg_s2(:,1)>0)];
            TT3(kkk) = T_num(ii);
        end
                once_t = 0;  


        if length(AGG_RBC_No)>0
        stack_uu(:,ii) = uu;
        stack_shear(:,ii) = SR;
        stack_acc(:,ii) = accel;
        end 
       
       figure(2)
       set(gcf,'position', [0 169 1200 300])
     
       imagesc(P_R_x,P_R_y,u(:,:,ii))
        colormap(jet)
        hold on 
        DrawCircle(P_R(:,1),P_R(:,2),R, 100, 'r')
        DrawCircle(P_R(AGG_RBC_No,1),P_R(AGG_RBC_No,2),R, 100, 'k.')
        DrawCircle(P_R(1,1),P_R(1,2),R/1, 100, 'g')
        DrawCircle(P_R(P_R_win2,1),P_R(P_R_win2,2),R/1, 100, 'b.');
%         text(P_R(AGG_RBC_No,1),P_R(AGG_RBC_No,2),num2str(uu(AGG_RBC_No)),'FontSize',10);
        plot(P_R_x,vessel_elasticity_up,'w','linewidth',3)
%         plot(P_R_x,vessel_elasticity_pre,'g','linewidth',3)
        plot(P_R_x,vessel_elasticity_down,'w','linewidth',3)
        xlabel('Vessel Length [m]')
        ylabel('Tube Depth [m]')
        title(['Time : ' num2str(T_num(ii)) 's'])
        caxis([0 0.045])
        colorbar()
        
         M=getframe(figure(2));
         writeVideo(aviobj,M)   
         kkk = kkk+1;
        
     end
%    close;
 
end