function [AGG_RBC_No2  win_aggRBC, win_field_data] = cjddk_window_rage(particle, agg_data,win_s, win_e,P_R_x)
% example.
% particle = P_R
% agg_data = AGG_RBC_No 
% win_s = win_start
% win_e = win_end
% data = [dv w_flow T*t k_flow in_t v_mean Bound_y];
% P_R_x  = P_R_x(1,:);

if win_s<win_e
    AGG_RBC_No2 = particle(agg_data,1)>=win_s& particle(agg_data,1)<=win_e;
    AGG_RBC_No2 = AGG_RBC_No2.*agg_data; 
    AGG_RBC_No2(AGG_RBC_No2==0)=[];
    win_field = find(P_R_x(1,:)>= win_s & P_R_x(1,:)<=win_e);
    win_aggRBC = particle(AGG_RBC_No2,:); 
else
    
    AGG_RBC_No2 = particle(agg_data,1)>=win_s| particle(agg_data,1)<=win_e;
    AGG_RBC_No2 = AGG_RBC_No2.*agg_data; 
    AGG_RBC_No2(AGG_RBC_No2==0)=[];
    win_field_1 = find(P_R_x>= win_s); win_field_2 = find(P_R_x<=win_e);
    win_field = [win_field_1 win_field_2];
    
    win_agg1 = particle(AGG_RBC_No2,1)<= win_e;
    win_agg1 = win_agg1.*AGG_RBC_No2;
    win_agg1(win_agg1==0)=[];
    
    win_agg2 = particle(AGG_RBC_No2,1)>=win_s;
    win_agg2 = win_agg2.*AGG_RBC_No2;
    win_agg2(win_agg2==0)=[];
    particle(win_agg1,1) = particle(win_agg1,1)+P_R_x(1,end);  
   
    win_aggRBC = particle([win_agg1;win_agg2],:);
end
win_field_data = win_field; 

 
 %% cacluation parameters field of window area
 
end 