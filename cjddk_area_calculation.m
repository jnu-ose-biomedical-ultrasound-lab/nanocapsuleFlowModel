function [area_data ] = cjddk_area_calculation(rage_parameter,x_length,P_R_y,win_rage,value,i)
% rage_parameter = win_acc_data;
% win_rage = win_acc;
% value = num2str('a')
% %% shear rate parameter
% rage_parameter = win_shear_data;
% win_rage = win_shear;
% %% temporal rate parameter
% rage_parameter = temp_rage;
% win_rage = win_temp;


 para_full(1,:)= [x_length(1) x_length(end) x_length(end) x_length(1) x_length(1)];
 para_full(2,:) = [P_R_y(1) P_R_y(1) P_R_y(end) P_R_y(end) P_R_y(1)];
 Area_full =  polyarea(para_full(1,:),para_full(2,:));

%% acceleration
if value == num2str('a')

%    for i = 10:3000:70000
    para_data   = contour( x_length, P_R_y,win_rage,[rage_parameter(i) rage_parameter(i)]);            
    para_data(:,para_data(1,:)==rage_parameter(i)) =[]; 
    para_data(2,end) = para_data(2,1);
    area_data = polyarea(para_data(1,:),para_data(2,:)); 
%      figure(1)
%      plot(para_data(1,:),para_data(2,:),'r')
%      hold on
%    end
elseif value == num2str('s')
         
        
            para_data   = contour( x_length, P_R_y,win_rage,[rage_parameter(i) rage_parameter(i)]);            
            para_data(:,para_data(1,:)==rage_parameter(i)) =[]; 
            for ii = 1:size(para_data,2)-1
                if para_data(2,ii)*para_data(2,ii+1)<0
                M1 = ii;
                end
            end
            para_data(1,1) =   para_full(1,2);
            para_data(1,M1)  = para_full(1,1);
            para_data(1,M1+1)  = para_full(1,4);
            para_data(1,end) = para_full(1,3);
            para_data(:,end+1) = para_data(:,1); 
            
            area_data = polyarea(para_data(1,:),para_data(2,:));
            
%            
         
else
       
            para_data   = contour( x_length, P_R_y,win_rage,[rage_parameter(i-1) rage_parameter(i-1)]);    
            hold on
            para_out = contour( x_length, P_R_y,win_rage,[rage_parameter(i) rage_parameter(i)]);
            para_data(:,para_data(1,:)==rage_parameter(i-1)) =[]; para_out(:,para_out(1,:)==rage_parameter(i)) =[];
            [M1] = find(para_data(2,:)==P_R_y(end));
            [M2] = find(para_out(2,:)==P_R_y(end));
            
            if rage_parameter(i-1)<0&rage_parameter(i)<0
            para_data(:,M1(1)) = para_full(:,3);
            para_data(:,M1(1)-1) = para_full(:,2);
            para_out(:,M2(1)) = para_full(:,3);
            para_out(:,M2(1)-1) = para_full(:,2);
            added_in = [size(para_data,2)+1 size(para_data,2)+2 size(para_data,2)+3];
            added_out = [size(para_out,2)+1 size(para_out,2)+2 size(para_out,2)+3];
            para_data(:,added_in) = [para_full(:,4) para_full(:,1) para_data(:,1)]; 
            para_out(:,added_out) = [para_full(:,4) para_full(:,1) para_out(:,1)]; 
           
            %% calculate the area
            Area_in = polyarea(para_data(1,:),para_data(2,:));
            Area_out =  polyarea(para_out(1,:),para_out(2,:));
            Area = abs(Area_out-Area_in);
            
            elseif rage_parameter(i-1)>0 & rage_parameter(i)>0
            para_data(:,M1(1)) = para_full(:,4);
            para_data(:,M1(1)-1) = para_full(:,1);
            para_out(:,M2(1)) = para_full(:,4);
            para_out(:,M2(1)-1) = para_full(:,1);
            added_in = [size(para_data,2)+1 size(para_data,2)+2 size(para_data,2)+3];
            added_out = [size(para_out,2)+1 size(para_out,2)+2 size(para_out,2)+3];
            para_data(:,added_in) = [para_full(:,3) para_full(:,2) para_data(:,1)]; 
            para_out(:,added_out) = [para_full(:,3) para_full(:,2) para_out(:,1)]; 
            Area_in =  polyarea(para_data(1,:),para_data(2,:));
            Area_out =  polyarea(para_out(1,:),para_out(2,:));
            Area = abs( Area_out-Area_in);
            
            elseif rage_parameter(i-1)<0 & rage_parameter(i)>0
            para_data(:,M1(1)) = para_full(:,3);
            para_data(:,M1(1)-1) = para_full(:,2);
            para_out(:,M2(1)) = para_full(:,4);
            para_out(:,M2(1)-1) = para_full(:,1);
            added_in = [size(para_data,2)+1 size(para_data,2)+2 size(para_data,2)+3];
            added_out = [size(para_out,2)+1 size(para_out,2)+2 size(para_out,2)+3];
            para_data(:,added_in) = [para_full(:,4) para_full(:,1) para_data(:,1)]; 
            para_out(:,added_out) = [para_full(:,3) para_full(:,2) para_out(:,1)];
            Area_in =   Area_full-polyarea(para_data(1,:),para_data(2,:));
            Area_out =  Area_full-polyarea(para_out(1,:),para_out(2,:));
            Area = abs( Area_full-(Area_out+Area_in));
            end
            area_data = [Area; Area_in; Area_out;Area_full]; %mm2
end


end
