function [ stereo_TD ] = stereo_matching_bp( varargin )
%   Matching stereo data captured by two event cameras
%   TD contain the two events data of the two camera
%
% TAKES IN:
%   'TD'
%       A struct of TD events with format:
%           TD(:,1) =  time stamps
%           TD(:,2) =  pixel X locations
%           TD(:,3) =  pixel Y locations
%           TD(:,4) =  event polarity
%           TD(:,5) =  left or right label left is 0 right is 1
%           TD(:,6) =  disparity
%
%
% RETURNS:f
%    stereo_TD.ts  timestamps
%    stereo_TD.x   x
%    stereo_TD.y   y
%    stereo_TD.p   disparity
%    stereo_TD.disparity_gt    depth_ground truth
%
%
% written by Xie zhen - may 2017

TD = varargin{1};
groundtruth_bool = varargin{2};   % default is 1, set to 0 if the input data has no groundtruth disparity


%  parameters for EMP 
width = 240;
height = 180;
temporal_var_Thre = 3e3;       % temporal criterion Threshold
special_var_Thre = 3;          % spatial  criterion Threshold
wradius = 1;                   % widow radius for the epipolar line the event
LABEL = 50;                    % maximum disparity
THETA = 1; %-0.1               % activation function threshold
DISC_K = 1.7;                  % truncation of discontinuity cost
DATE_TRUNC = 5;                % truncation of matching cost 
BP_ITERATIONS = 1;             % iteration times 
LAMBDA = 1;                    % weighting of data costd
WRADIUS_TS = 2;                % messages update in this radius.

old_pixel_threshold  = 20e3;   % old event Threshold
old_mrf_threshold  = 10e3;     % old event Threshold for update MRF related with the speed of the movement

% last spike map  store the latest ts and p 
last_spike_time_right = nan(height,width,2);

%define the struct of the MRF field  and set all messages to zero
msg_up = zeros( height, width, LABEL);
msg_down = zeros( height, width, LABEL);
msg_left = zeros( height, width, LABEL);
msg_right = zeros( height, width, LABEL);
msg_data = zeros( height, width, LABEL);
msg_ts = zeros( height, width, LABEL)+TD(1,1);
disparity = 0;

%Init para of matching
stereo_TD.x = [];
stereo_TD.y = [];
stereo_TD.p = [];
stereo_TD.ts = [];
stereo_TD.rx = [];
stereo_TD.ry = [];
if groundtruth_bool
    stereo_TD.disparity_gt = [];
end

%for each event
for event_index =1:length(TD)
    event_index
    current_event_index = event_index;
    current_ts = TD(current_event_index,1);
    current_x = TD(current_event_index,2);
    current_y = TD(current_event_index,3);
    current_p = TD(current_event_index,4);
    if groundtruth_bool
         current_disparity = TD(current_event_index,6);
    end

    
    if (TD(current_event_index,5) ==0 )
        min_error_row = [];
        msg_data_bool = 0;
        for l_i = 0:LABEL-1            
            if (current_x - l_i > 1 && current_y-wradius > 1 && current_y+wradius < height)
                sub_range = current_ts - last_spike_time_right(current_y-wradius:current_y+wradius,current_x-l_i,1);
                sub_range_p = current_p - last_spike_time_right(current_y-wradius:current_y+wradius,current_x-l_i,2);
                sub_range(sub_range_p~=0) = nan;                           % considert the polarity 
                sub_range(abs(sub_range)> old_pixel_threshold) = nan;      % mark old events as invalid  
                %calculate the time difference
                error_lr_t = abs(sub_range)/temporal_var_Thre;
                %calculate the distance to the epipolar line
                error_lr_d = abs(-wradius:wradius)/special_var_Thre;
                error_lr_row = error_lr_t(:) + error_lr_d(:);
                index = ~isnan(error_lr_row);
                [min_error_row,~] = min(error_lr_row(index ~= 0));
            end
            if (~isempty(min_error_row))
                msg_data(current_y,current_x,l_i+1) = min_error_row;
                msg_data_bool = msg_data_bool+1;
            else
                msg_data(current_y,current_x,l_i+1) = DATE_TRUNC;
            end
        end
        if (msg_data_bool>1)  
            if (current_y-WRADIUS_TS >1 && current_y+WRADIUS_TS < height && current_x-WRADIUS_TS>1 && current_x+WRADIUS_TS< width)
                for i = 1:LABEL                    
                    bool_mrf_update  = abs(current_ts - msg_ts(current_y-WRADIUS_TS:current_y+WRADIUS_TS,current_x-WRADIUS_TS:current_x+WRADIUS_TS,i))< old_mrf_threshold;
                    msg_up( current_y-WRADIUS_TS:current_y+WRADIUS_TS,current_x-WRADIUS_TS:current_x+WRADIUS_TS, i) = msg_up( current_y-WRADIUS_TS:current_y+WRADIUS_TS,current_x-WRADIUS_TS:current_x+WRADIUS_TS, i).*bool_mrf_update;
                    msg_down( current_y-WRADIUS_TS:current_y+WRADIUS_TS,current_x-WRADIUS_TS:current_x+WRADIUS_TS, i) =  msg_down( current_y-WRADIUS_TS:current_y+WRADIUS_TS,current_x-WRADIUS_TS:current_x+WRADIUS_TS, i).*bool_mrf_update;
                    msg_left( current_y-WRADIUS_TS:current_y+WRADIUS_TS,current_x-WRADIUS_TS:current_x+WRADIUS_TS, i) = msg_left( current_y-WRADIUS_TS:current_y+WRADIUS_TS,current_x-WRADIUS_TS:current_x+WRADIUS_TS, i).*bool_mrf_update;
                    msg_right( current_y-WRADIUS_TS:current_y+WRADIUS_TS,current_x-WRADIUS_TS:current_x+WRADIUS_TS, i) = msg_right( current_y-WRADIUS_TS:current_y+WRADIUS_TS,current_x-WRADIUS_TS:current_x+WRADIUS_TS, i).*bool_mrf_update;
                    msg_data( current_y-WRADIUS_TS:current_y+WRADIUS_TS,current_x-WRADIUS_TS:current_x+WRADIUS_TS, i) = msg_data( current_y-WRADIUS_TS:current_y+WRADIUS_TS,current_x-WRADIUS_TS:current_x+WRADIUS_TS, i).*bool_mrf_update;
                end
            end
            %% bp iteration
            for iter = 1: BP_ITERATIONS
                % compute message  
                % RIGHT
                for r = 0:WRADIUS_TS
                    for r_y = -r:r
                        for r_x = -r:r
                            y = current_y+ r_y;
                            x = current_x+ r_x;
                            if (x > 1 && x <width  && y > 1 && y < height)  
                                % aggregate and find min
                                min_val = Inf;
                                for value = 1:LABEL
                                    msg_left(y,x+1,value) = LAMBDA*msg_data(y,x,value)+ msg_left(y,x,value)+ msg_up(y,x,value)+msg_down(y,x,value);
                                    if (msg_left(y,x+1,value) < min_val)
                                        min_val = msg_left(y,x+1,value);
                                    end
                                end
                                %dt
                                for q = 2:LABEL-1
                                    prev = msg_left(y,x+1,q-1) + 1.0;
                                    if (prev < msg_left(y,x+1,q))
                                        msg_left(y,x+1,q) = prev;
                                    end
                                end
                                for q = LABEL-1:-1:1   
                                    prev =  msg_left(y,x+1,q+1) + 1.0;
                                    if (prev <  msg_left(y,x+1,q))
                                        msg_left(y,x+1,q) = prev;
                                    end
                                end
                                % truncate
                                min_val = min_val+DISC_K;  %% tune DISC_K for smooth 
                                for value = 1:LABEL
                                    if (min_val < msg_left(y,x+1,value))
                                        msg_left(y,x+1,value) = min_val;
                                    end
                                end       
                                % normalize                              
                                val = 0;
                                for  value = 1:LABEL
                                    val = val + msg_left(y,x+1,value);
                                end
                                val = val/LABEL;
                                for value = 1:LABEL
                                    msg_left(y,x+1,value) = msg_left(y,x+1,value)- val;                                    
                                end                                                            
                            end
                        end
                    end
                end  
                %LEFT
                for r = 0:WRADIUS_TS
                    for r_y = -r:r
                        for r_x = -r:r
                            y = current_y+ r_y;
                            x = current_x+ r_x;
                            if (x > 1 && x <width  && y > 1 && y < height)
                                % aggregate and find min
                                min_val = Inf;
                                for value = 1:LABEL
                                    msg_right(y,x-1,value) = LAMBDA* msg_data(y,x,value)+ msg_right(y,x,value)+ msg_up(y,x,value)+msg_down(y,x,value);
                                    if (msg_right(y,x-1,value) < min_val)
                                        min_val = msg_right(y,x-1,value);
                                    end
                                end
                                %dt
                                for q = 2:LABEL-1
                                    prev = msg_right(y,x-1,q-1) + 1.0;
                                    if (prev < msg_right(y,x-1,q))
                                        msg_right(y,x-1,q) = prev;
                                    end
                                end
                                for q = LABEL-1:-1:1
                                    prev =  msg_right(y,x-1,q+1) + 1.0;
                                    if (prev <  msg_right(y,x-1,q))
                                        msg_right(y,x-1,q) = prev;
                                    end
                                end
                                % truncate
                                min_val = min_val+DISC_K;
                                for value = 1:LABEL
                                    if (min_val < msg_right(y,x-1,value))
                                        msg_right(y,x-1,value) = min_val;
                                    end
                                end
                                % normalize
                                val = 0;
                                for  value = 1:LABEL
                                    val = val + msg_right(y,x-1,value);
                                end                                
                                val = val/LABEL;
                                for value = 1:LABEL
                                    msg_right(y,x-1,value) = msg_right(y,x-1,value)- val;
                                end
                            end
                        end
                    end
                end            
                % DOWN
                for r = 0:WRADIUS_TS
                    for r_y = -r:r
                        for r_x = -r:r
                            y = current_y+ r_y;
                            x = current_x+ r_x;
                            if (x > 1 && x <width  && y > 1 && y < height) 
                                % aggregate and find min
                                min_val = Inf;
                                for value = 1:LABEL
                                    msg_up(y+1,x,value) = LAMBDA*msg_data(y,x,value)+ msg_left(y,x,value)+ msg_up(y,x,value)+msg_right(y,x,value);
                                    if (msg_up(y+1,x,value) < min_val)
                                        min_val = msg_up(y+1,x,value);
                                    end
                                end
                                %dt
                                for q = 2:LABEL-1
                                    prev = msg_up(y+1,x,q-1) + 1.0;
                                    if (prev < msg_up(y+1,x,q))
                                        msg_up(y+1,x,q) = prev;
                                    end
                                end
                                for q = LABEL-1:-1:1
                                    prev =  msg_up(y+1,x,q+1) + 1.0;
                                    if (prev <  msg_up(y+1,x,q))
                                        msg_up(y+1,x,q) = prev;
                                    end
                                end
                                % truncate
                                min_val = min_val+DISC_K;
                                for value = 1:LABEL
                                    if (min_val < msg_up(y+1,x,value))
                                        msg_up(y+1,x,value) = min_val;
                                    end
                                end
                                % normalize
                                val = 0;
                                for  value = 1:LABEL
                                    val = val + msg_up(y+1,x,value);
                                end
                                
                                val = val/LABEL;
                                for value = 1:LABEL
                                    msg_up(y+1,x,value) = msg_up(y+1,x,value)- val;
                                end
                            end
                        end
                    end
                end
                %UP
                for r = 0:WRADIUS_TS
                    for r_y = -r:r
                        for r_x = -r:r
                            new_msg = zeros(LABEL,1);
                            y = current_y+ r_y;
                            x = current_x+ r_x;
                            if (x > 1 && x <width  && y > 1 && y < height)
                                % aggregate and find min
                                min_val = Inf;
                                for value = 1:LABEL
                                    msg_down(y-1,x,value) = LAMBDA*msg_data(y,x,value)+ msg_left(y,x,value)+ msg_right(y,x,value)+msg_down(y,x,value);
                                    if (msg_down(y-1,x,value) < min_val)
                                        min_val = msg_down(y-1,x,value);
                                    end
                                end
                                %dt
                                for q = 2:LABEL-1
                                    prev = msg_down(y-1,x,q-1) + 1.0;
                                    if (prev < msg_down(y-1,x,q))
                                        msg_down(y-1,x,q) = prev;
                                    end
                                end
                                for q = LABEL-1:-1:1
                                    prev =  msg_down(y-1,x,q+1) + 1.0;
                                    if (prev <  msg_down(y-1,x,q))
                                        msg_down(y+1,x,q) = prev;
                                    end
                                end
                                % truncate
                                min_val = min_val+DISC_K;
                                for value = 1:LABEL
                                    if (min_val < msg_down(y-1,x,value))
                                        msg_down(y-1,x,value) = min_val;
                                    end
                                end
                                % normalize
                                val = 0;
                                for  value = 1:LABEL
                                    val = val + msg_down(y-1,x,value);
                                end
                                val = val/LABEL;
                                for value = 1:LABEL
                                     msg_down(y-1,x,value) = msg_down(y-1,x,value)- val;
                                end
                            end
                        end
                    end
                end
            end
                % updatetime
                if (current_y-WRADIUS_TS >=1 && current_y+WRADIUS_TS <= height && current_x-WRADIUS_TS>=1 && current_x+WRADIUS_TS<= width)
                msg_ts(current_y-WRADIUS_TS:current_y+WRADIUS_TS,current_x-WRADIUS_TS:current_x+WRADIUS_TS,:)= current_ts;
                end
                %Finds the MAP assignment as well as calculating the energy
                %MAP assignment
                y = current_y;
                x = current_x;
                best = Inf;
                for j =1: LABEL
                    cost_event =  msg_left(y,x,j) + msg_right(y,x,j) + msg_up(y,x,j) + msg_down(y,x,j) + LAMBDA*msg_data(y,x,j);
                    if (cost_event < best)
                        best = cost_event;
                        disparity = j-1;
                    end
                end
                if best < THETA  && disparity ~=0
                    stereo_TD.x = [stereo_TD.x ;current_x];
                    stereo_TD.y = [stereo_TD.y ;current_y];
                    stereo_TD.p  = [stereo_TD.p ; disparity];
                    stereo_TD.ts = [stereo_TD.ts ; current_ts];
                    if groundtruth_bool
                        stereo_TD.disparity_gt = [stereo_TD.disparity_gt ; current_disparity];
                    end
                end
        end
    else
         last_spike_time_right(TD(event_index,3),TD(event_index,2),1) = TD(event_index,1);
         last_spike_time_right(TD(event_index,3),TD(event_index,2),2) = TD(event_index,4);
    end
    
end
end
