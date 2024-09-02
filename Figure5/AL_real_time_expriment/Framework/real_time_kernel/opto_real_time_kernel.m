function opto_real_time_kernel(TD)
close all;  
clc;  
dbstop if error;


DEBUG           = 0;
% EXPERIMENT_TYPE = 'grid_search1'; old version do not use (2019_12_09 MJC)
EXPERIMENT_TYPE = 'grid_search3'; 
% EXPERIMENT_TYPE = 'electrical_stimulation'; 
% EXPERIMENT_TYPE = 'open_loop';
% EXPERIMENT_TYPE = 'pid_controller';
% EXPERIMENT_TYPE = 'cross_entropy';
% EXPERIMENT_TYPE = 'bayesian_optimization';
% EXPERIMENT_TYPE = 'theta_optimization';

%%%%%%%%%
% Created on 12/09/19 known to work
%%%%%%%%%
% EXPERIMENT_TYPE = 'gamma_maximization'; % for gamma optimization
% EXPERIMENT_TYPE = 'random_nested_pulse_train'; % SANG'S CODE
% EXPERIMENT_TYPE = 'gp_gamma_uncertainty_sampling'; %active learning of gamma response vs. standard train param
% EXPERIMENT_TYPE = 'gp_gamma_uncertainty_sampling_2';
% EXPERIMENT_TYPE = 'nested_grid_search'; %Code by eric to do grid NPT parameters
% EXPERIMENT_TYPE = 'nested_random_search'; %Random NPT parameters
% EXPERIMENT_TYPE = 'random_grid_search';                                         %search of NPT stimulation
% EXPERIMENT_TYPE = 'random_sine_search';
% EXPERIMENT_TYPE = 'grid_sine_search';
% EXPERIMENT_TYPE = 'grid_poisson_search';
% EXPERIMENT_TYPE = 'random_poisson_search';
% EXPERIMENT_TYPE = 'grid_search_seizure_2';
% EXPERIMENT_TYPE = 'grid_search_seizure_3';
% EXPERIMENT_TYPE = 'recording';
% EXPERIMENT_TYPE = 'active_learning';



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
tic
switch EXPERIMENT_TYPE 
    case 'grid_search1'
        [optimization_object, ~, metric_objects, display_objects] = ...
            opto_configure_grid_search_1(TD, DEBUG); %1 is for NOR, conventional grid search.. %3 is for train pulse
    case 'grid_search3'
        [optimization_object, ~, metric_objects, display_objects] = ...
            opto_configure_grid_search_3(TD, DEBUG); %1 is for NOR, conventional grid search.. %3 is for train pulse
    case 'electrical_stimulation'
        [optimization_object, ~, metric_objects] = ...
            opto_configure_electrical_stimulation(TD, DEBUG); %1 is for NOR, conventional grid search.. %3 is for train pulse
    case 'cross_entropy'
        [optimization_object, ~, metric_objects] = ...
            opto_configure_cross_entropy(TD, DEBUG);
    case 'pid_controller'
        [optimization_object, ~, metric_objects] = ...
            configure_pid_control(TD, DEBUG);
        %         display_objects = metric_objects;
    case 'bayesian_optimization'
        [optimization_object, ~, metric_objects] = ...
            opto_configure_bayesian_optimization_controller(TD, DEBUG);
        
    case 'bayesian_optimization2'
        [optimization_object, ~, metric_objects] = ...
            opto_configure_bayesian_optimization_controller2(TD, DEBUG);
%         display_objects = metric_objects;
    case 'theta_optimization'
        [optimization_object, ~, metric_objects] = ...
            opto_configure_bayesian_optimization_controller_110918_archive(TD, DEBUG);
    case 'gamma_maximization'
        [optimization_object, ~, metric_objects] = ...
            opto_configure_maximize_gamma_single_frequency(TD, DEBUG);
    case 'random_nested_pulse_train'
        [optimization_object, ~, metric_objects] = ...
            opto_configure_random_nested_pulse_train(TD, DEBUG);
    case 'gp_gamma_uncertainty_sampling'
        [optimization_object, ~, metric_objects] = ...
            opto_configure_gamma_gp_uncertainty_sampling(TD, DEBUG);
    case 'gp_gamma_uncertainty_sampling_2'
        [optimization_object, ~, metric_objects] = ...
            opto_configure_gamma_gp_uncertainty_sampling_2(TD, DEBUG);    
    case 'nested_random_search'
         [optimization_object, ~, metric_objects, display_objects] = ...
            opto_configure_nested_search(TD, DEBUG); 
    case 'nested_grid_search'
         [optimization_object, ~, metric_objects, display_objects] = ...
            opto_configure_nested_grid_search(TD, DEBUG);     
    case 'random_grid_search'
         [optimization_object, ~, metric_objects, display_objects] = ...
             opto_configure_random_grid_search(TD, DEBUG);
    case 'random_sine_search'
        [optimization_object, ~, metric_objects, display_objects] = ...
            opto_configure_sine_search(TD, DEBUG);
    case 'grid_sine_search'
        [optimization_object, ~, metric_objects, display_objects] = ...
            opto_configure_grid_sine_search(TD, DEBUG);
    case 'grid_poisson_search'
        [optimization_object, ~, metric_objects, display_objects] = ...
            opto_configure_grid_poisson(TD, DEBUG);
    case 'random_poisson_search'
        [optimization_object, ~, metric_objects, display_objects] = ...
            opto_configure_random_poisson_search(TD, DEBUG);
    case 'grid_search_seizure_2'
        [optimization_object, ~, metric_objects, display_objects] = ...
            opto_configure_grid_search_sz_monitor_2(TD, DEBUG);
    case 'grid_search_seizure_3'
        [optimization_object, ~, metric_objects, display_objects] = ...
            opto_configure_grid_search_sz_monitor_3(TD, DEBUG);
    case 'recording'
        [optimization_object, ~, metric_objects, display_objects] = ...
            opto_configure_recording(TD, DEBUG);
    case 'active_learning'
        [optimization_object, ~, metric_objects] = ...
            opto_configure_mo_al_gamma(TD, DEBUG);
end


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
    if exist('display_objects', 'var') && ~isempty(display_objects)
        display_objects{c1}.update_buffer(new_data, 0);
%         realtime_metric_display(display_objects, new_data); 
        display_objects{c1}.display_metric();
    end
    % Video
%     video_objects{1}.write_recording();
    
end  % End While
% video_objects{1}.close_recording();
TD.SetSysMode(0)
end

 function t = get_current_time(TD,device_name, sampling_frequency)
     t = TD.ReadTargetVEX([device_name '.current_time'], 0, 1, 'I32', 'F64')/sampling_frequency;
 end
        
function clean_up(TD)
    TD.SetSysMode(0);
    
    fclose('all');
%     close all;
end
