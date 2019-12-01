function [para_in para_out para_full result_data para_result area_data] = cjddk_sample_space(rage_parameter,win_aggRBC,x_length,P_R_y,win_rage,value ,NUM,i,time)
% rage_parameter = accel_rage;
% win_rage = win_acc;
% value = num2str('a')
% %% shear rate parameter
% rage_parameter = shear_rage;
% win_rage = win_shear;
% %% temporal rate parameter
% rage_parameter = temp_rage;
% win_rage = win_temp;


 para_full(1,:)= [x_length(1) x_length(end) x_length(end) x_length(1) x_length(1)];
 para_full(2,:) = [P_R_y(1) P_R_y(1) P_R_y(end) P_R_y(end) P_R_y(1)];
 Area_full =  polyarea(para_full(1,:),para_full(2,:));
 xq = win_aggRBC(:,1); yq = win_aggRBC(:,2);
%% acceleration
if value == num2str('a')

            figure(1)
            para_in   = contour( x_length, P_R_y,win_rage,[rage_parameter(i-1) rage_parameter(i-1)]);            
            hold on
            para_out = contour( x_length, P_R_y,win_rage,[rage_parameter(i) rage_parameter(i)]);
            hold off
    para_in(:,para_in(1,:)==rage_parameter(i-1)) =[]; para_out(:,para_out(1,:)==rage_parameter(i)) =[];
    Area_in = polyarea(para_in(1,:),para_in(2,:));
    Area_out =  polyarea(para_out(1,:),para_out(2,:));
    [result_in re_1] = inpolygon(xq,yq,para_in(1,:),para_in(2,:));
    [result_out,re_2] = inpolygon(xq,yq,para_out(1,:),para_out(2,:));
    [result_full re_u] = inpolygon(xq,yq,para_full(1,:),para_full(2,:));
 
            if rage_parameter(i-1)*rage_parameter(i)>0
                if  rage_parameter(i-1) <rage_parameter(i)
                    para_result = result_out - result_in;  
                    Area = abs(Area_out-Area_in);
                else
                     para_result = result_in - result_out;    
                     Area = abs(Area_in-Area_out);
                end
            else
                para_result  = result_full-(result_in+result_out);
               Area = Area_full-(Area_in+Area_out);
            end
    para_result = logical(para_result);
           if i == 3
            result_data(:,i-1) = [rage_parameter(2); numel(xq(result_in))/(Area_in); numel(xq(result_in)) ];
            result_data(:,i) = [rage_parameter(3); numel(xq(para_result))/(Area); numel(xq(para_result)) ];
         else
            result_data(:,i) = [rage_parameter(i); numel(xq(para_result))/(Area);numel(xq(para_result))];
            end
    area_data = [Area; Area_in; Area_out;Area_full]; %mm2
elseif value == num2str('s')
            figure(1)
            para_in   = contour( x_length, P_R_y,win_rage,[rage_parameter(i-1) rage_parameter(i-1)]);            
            hold on
            para_out = contour( x_length, P_R_y,win_rage,[rage_parameter(i) rage_parameter(i)]);
            hold off
            para_in(:,para_in(1,:)==rage_parameter(i-1)) =[]; para_out(:,para_out(1,:)==rage_parameter(i)) =[];
            [M1] = find(para_in(2,:)==P_R_y(end));
            [M2] = find(para_out(2,:)==P_R_y(end));
            if length(M1)==0&length(M2)==0
              para_in(:,end+1) = para_in(:,1);
              para_out(:,end+1) = para_out(:,1);
            elseif length(M1)==0|length(M2)==0
            para_out(:,M2(1)) = para_full(:,4);
            para_out(:,M2(1)-1) = para_full(:,1);  
            added_out = [size(para_out,2)+1 size(para_out,2)+2 size(para_out,2)+3];
            para_out(:,added_out) = [para_full(:,3) para_full(:,2) para_out(:,1)]; 
            para_in(:,end+1) = para_in(:,1);
            else
            para_in(:,M1(1)) = para_full(:,4);
            para_in(:,M1(1)-1) = para_full(:,1);
            para_out(:,M2(1)) = para_full(:,4);
            para_out(:,M2(1)-1) = para_full(:,1);
            added_in = [size(para_in,2)+1 size(para_in,2)+2 size(para_in,2)+3];
            added_out = [size(para_out,2)+1 size(para_out,2)+2 size(para_out,2)+3];
            para_in(:,added_in) = [para_full(:,3) para_full(:,2) para_in(:,1)]; 
            para_out(:,added_out) = [para_full(:,3) para_full(:,2) para_out(:,1)]; 
            end
            Area_in = polyarea(para_in(1,:),para_in(2,:));
            Area_out =  polyarea(para_out(1,:),para_out(2,:));
            Area = abs(Area_out-Area_in);
            [result_in re_1] = inpolygon(xq,yq,para_in(1,:),para_in(2,:));
            [result_out,re_2] = inpolygon(xq,yq,para_out(1,:),para_out(2,:));
            [result_full re_u] = inpolygon(xq,yq,para_full(1,:),para_full(2,:));
            para_result  = (result_full-result_in).*result_out;
            para_result = logical(para_result);
      if i == 3
            result_data(:,i-1) = [rage_parameter(2); numel(xq(result_in))/(Area_in); numel(xq(result_in)) ];
            result_data(:,i) = [rage_parameter(3); numel(xq(para_result))/(Area); numel(xq(para_result)) ];
         else
            result_data(:,i) = [rage_parameter(i); numel(xq(para_result))/(Area);numel(xq(para_result))];
            end
    area_data = [Area; Area_in; Area_out;Area_full]; %mm2
% figure
% plot(para_in(1,:),para_in(2,:),'r',para_out(1,:),para_out(2,:),'g',para_full(1,:),para_full(2,:),'--k')
% plot(para_in(1,:),para_in(2,:),'r',para_out(1,:),para_out(2,:),'g')
% hold on
% DrawCircle(win_aggRBC(:,1),win_aggRBC(:,2),R/1, 100, 'b.')
% DrawCircle(win_aggRBC(para_result,1),win_aggRBC(para_result,2),R/1, 100, 'r.')
else
       
            para_in   = contour( x_length, P_R_y,win_rage,[rage_parameter(i-1) rage_parameter(i-1)]);    
            hold on
            para_out = contour( x_length, P_R_y,win_rage,[rage_parameter(i) rage_parameter(i)]);
            para_in(:,para_in(1,:)==rage_parameter(i-1)) =[]; para_out(:,para_out(1,:)==rage_parameter(i)) =[];
            [M1] = find(para_in(2,:)==P_R_y(end));
            [M2] = find(para_out(2,:)==P_R_y(end));
            
            if rage_parameter(i-1)<0&rage_parameter(i)<0
            para_in(:,M1(1)) = para_full(:,3);
            para_in(:,M1(1)-1) = para_full(:,2);
            para_out(:,M2(1)) = para_full(:,3);
            para_out(:,M2(1)-1) = para_full(:,2);
            added_in = [size(para_in,2)+1 size(para_in,2)+2 size(para_in,2)+3];
            added_out = [size(para_out,2)+1 size(para_out,2)+2 size(para_out,2)+3];
            para_in(:,added_in) = [para_full(:,4) para_full(:,1) para_in(:,1)]; 
            para_out(:,added_out) = [para_full(:,4) para_full(:,1) para_out(:,1)]; 
            [result_in re_1] = inpolygon(xq,yq,para_in(1,:),para_in(2,:));
            [result_out,re_2] = inpolygon(xq,yq,para_out(1,:),para_out(2,:));
            [result_full re_u] = inpolygon(xq,yq,para_full(1,:),para_full(2,:));
            para_result  = (result_full-result_out).*result_in;
            %% calculate the area
            Area_in = polyarea(para_in(1,:),para_in(2,:));
            Area_out =  polyarea(para_out(1,:),para_out(2,:));
            Area = abs(Area_out-Area_in);
            
            elseif rage_parameter(i-1)>0 & rage_parameter(i)>0
            para_in(:,M1(1)) = para_full(:,4);
            para_in(:,M1(1)-1) = para_full(:,1);
            para_out(:,M2(1)) = para_full(:,4);
            para_out(:,M2(1)-1) = para_full(:,1);
            added_in = [size(para_in,2)+1 size(para_in,2)+2 size(para_in,2)+3];
            added_out = [size(para_out,2)+1 size(para_out,2)+2 size(para_out,2)+3];
            para_in(:,added_in) = [para_full(:,3) para_full(:,2) para_in(:,1)]; 
            para_out(:,added_out) = [para_full(:,3) para_full(:,2) para_out(:,1)]; 
            [result_in re_1] = inpolygon(xq,yq,para_in(1,:),para_in(2,:));
            [result_out,re_2] = inpolygon(xq,yq,para_out(1,:),para_out(2,:));
            [result_full re_u] = inpolygon(xq,yq,para_full(1,:),para_full(2,:));
            para_result  = result_full.*(result_out-result_in);
            
            Area_in =  polyarea(para_in(1,:),para_in(2,:));
            Area_out =  polyarea(para_out(1,:),para_out(2,:));
            Area = abs( Area_out-Area_in);
            elseif rage_parameter(i-1)<0 & rage_parameter(i)>0
            para_in(:,M1(1)) = para_full(:,3);
            para_in(:,M1(1)-1) = para_full(:,2);
            para_out(:,M2(1)) = para_full(:,4);
            para_out(:,M2(1)-1) = para_full(:,1);
            added_in = [size(para_in,2)+1 size(para_in,2)+2 size(para_in,2)+3];
            added_out = [size(para_out,2)+1 size(para_out,2)+2 size(para_out,2)+3];
            para_in(:,added_in) = [para_full(:,4) para_full(:,1) para_in(:,1)]; 
            para_out(:,added_out) = [para_full(:,3) para_full(:,2) para_out(:,1)]; 
            [result_in re_1] = inpolygon(xq,yq,para_in(1,:),para_in(2,:));
            [result_out,re_2] = inpolygon(xq,yq,para_out(1,:),para_out(2,:));
            [result_full re_u] = inpolygon(xq,yq,para_full(1,:),para_full(2,:));
            Area_in =   Area_full-polyarea(para_in(1,:),para_in(2,:));
            Area_out =  Area_full-polyarea(para_out(1,:),para_out(2,:));
            Area = abs( Area_full-(Area_out+Area_in));
            
            para_result  = result_full.*result_in.*result_out;
            end
            para_result = logical(para_result);
        if i == 3
            result_data(:,i-1) = [rage_parameter(2); numel(xq(~result_in))/(Area_in); numel(xq(~result_in)) ];
            result_data(:,i) = [rage_parameter(3); numel(xq(para_result))/(Area); numel(xq(para_result)) ];
        else
            result_data(:,i) = [rage_parameter(i); numel(xq(para_result))/(Area);numel(xq(para_result))];
        end
         area_data = [Area; Area_in; Area_out;Area_full]; %mm2
end


