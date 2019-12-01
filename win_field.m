%% the function calculation position of the particles in the ROI
% field : one period velocity field in the tube
% wave_velocity : the sinusoidal pulsatile flow velocity
% frequency : the sinusoidal pulsatile flow frequency at 60 bpm
% field_x : particle x position
% ** output data**
% win_start : window of ROI start point
% win_end :window of ROI end point
function  [win_start win_end] = win_field(field, wave_velocity, frequency,field_x);
%% example 
%  field = velocity_field;
%  wave_velocity = Vp_flow; 
%  frequency = f_flow;
% field_x = P_R_x(1,:);
 wave_length  = wave_velocity/frequency;
 Bound_x = 1.0*10^-3;
 [y x]=size(field);
 max_velue = max(max(field));
 [center_rate] = find(field(y/2,:) == max_velue);
 P_R_cent = field_x(1,center_rate(1));
 win_start = P_R_cent-wave_length/2;
 win_end = P_R_cent+wave_length/2;
 if win_start<0;
     win_start = win_start+Bound_x;
 elseif win_end> Bound_x;
     win_end = win_end-Bound_x;
 end
 
end
 