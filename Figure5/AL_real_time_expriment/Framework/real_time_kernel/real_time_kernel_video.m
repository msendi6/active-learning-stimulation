function real_time_kernel_video(TD)
close all;  
clc;  
dbstop if error;
ANIMAL_ID       = 'PDR008';
DEBUG           = 0;

% video_filename                              = 'Z:\Sang_video\test_NOR_nostim_1.avi';

temp = datestr(now,2);
temp(strfind(temp,'/')) =[]; % Today date
% [~, tank_name]          = fileparts(TD.GetTankName);
% block_name              = get_block_name(tank_name);


EXPERIMENT_TYPE = 'PD_grid';
% EXPERIMENT_TYPE = 'open_loop';
% EXPERIMENT_TYPE = 'pid_controller';
% EXPERIMENT_TYPE = 'cross_entropy';
% EXPERIMENT_TYPE = 'bayesian_optimization';
% EXPERIMENT_TYPE = 'bayesian_optimization_controller';

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Configure TDT Settings %
%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('TD', 'var')
    TD          = connect_to_TDT();
end

device_name = TD.GetDeviceName(0);
TD_FS       = TD.GetDeviceSF(device_name);
finishup    = onCleanup(@() clean_up(TD));

TD_BUFFER_TIME          = 1; % Seconds (Control policy checked at TD_BUFFER_TIME/2)
TD_READ_BUFFER_SIZE     = 16*floor(TD_FS*TD_BUFFER_TIME);
TD.SetTargetVal([device_name '.read_durr'], TD_READ_BUFFER_SIZE); % DOES NOT WORK! Need to modify circuit directly

% Get acquisition circuit variables
n_read_pts      = TD_READ_BUFFER_SIZE;
n_buff_pts      = n_read_pts/2; 
curindex        = TD.ReadTargetVEX([device_name '.read_index'], 0, 1, 'F32', 'F32');
buffer_offset   = n_buff_pts;

if DEBUG
    button = questdlg('You are starting in DEBUG mode','Are you sure you want to continue?');
    switch button 
        case 'Cancel'
            return;
        case 'No'
            return;
    end
    
    TD.SetSysMode(2);
else
    TD.SetSysMode(3);
end
pause(0.5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Configure Experiment Objects %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
optimization_object = [];

switch EXPERIMENT_TYPE 
    case 'grid_search'
        [optimization_object, ~, metric_objects] = ...
            configure_grid_search_video(TD, DEBUG, ANIMAL_ID);      
    case 'cross_entropy'
        [optimization_object, ~, metric_objects] = ...
            configure_cross_entropy(TD, DEBUG, ANIMAL_ID);
    case 'pid_controller'
        [optimization_object, ~, metric_objects] = ...
            configure_pid_control(TD, DEBUG);
%         display_objects = metric_objects;
    case 'bayesian_optimization'
        [optimization_object, ~, metric_objects] = ...
            configure_bayesian_optimization(TD, DEBUG);
%         display_objects = metric_objects;
    case 'bayesian_optimization_controller'
        [optimization_object, ~, metric_objects] = ...
            configure_bayesian_optimization_controller(TD, DEBUG);
    case 'PD_grid'
        [optimization_object, stim_manager, metric_objects] = ...
            configure_grid_search_PD(TD, DEBUG, ANIMAL_ID);
end


video_filename                                  = optimization_object.video_filename; %strcat('C:\Users\TDT\Desktop\PD Grid Search Videos\',animal_id,'_',temp,'_',experiment_name,block_name,'.avi');

video        = opto_video_recording(video_filename); 
video_objects = {video};




while ~optimization_object.optimization_done
    
    if buffer_offset == 0 % Second half of buffer

        % Check if second half of buffer is full
        while(curindex >= n_buff_pts)
            curindex = TD.ReadTargetVEX([device_name '.read_index'], 0, 1, 'F32', 'F32');
        end
        buffer_offset = n_buff_pts;
 
    else % First half of buffer

        % Check if first half of buffer is full
        while(curindex < n_buff_pts)
            curindex = TD.ReadTargetVEX([device_name '.read_index'], 0, 1, 'F32', 'F32');
        end
        buffer_offset = 0; 
        
        
        
    end
    
    % Read data from buffer and reshape

    
    
    new_data    = TD.ReadTargetVEX([device_name '.read_buff'], ...
        buffer_offset, n_buff_pts, 'F32', 'F32');   
    new_data    = reshape(new_data, 16, n_buff_pts/16)';
    
    % Update metrics
    for c1 = 1:size(metric_objects,2)
        metric_objects{c1}.update_buffer(new_data, get_current_time(TD, device_name, TD_FS));
    end
    
    % Iterate optimization step
    optimization_object.optimize();

    % Display metrics
    if exist('display_objects', 'var')
        realtime_metric_display(display_objects); 
    end
    % Video
    video_objects{1}.write_recording();
end  % End While

video_objects{1}.close_recording();
TD.SetSysMode(0)
end

 function t = get_current_time(TD,device_name, sampling_frequency)
     t = TD.ReadTargetVEX([device_name '.current_time'], 0, 1, 'I32', 'F64')/sampling_frequency;
 end
        
function clean_up(TD)
    TD.SetSysMode(0);
    fclose('all');
    close all;
end
