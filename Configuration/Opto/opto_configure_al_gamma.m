function [optimization_object, stimulation_manager, metric_objects] = opto_configure_mo_al_gamma( TD, DEBUG, TURBO)

%Adapted version of gamma_gp_uncertainty sampling: expanded for flexibility
%with multiple parameter spaces; incorporates hill climbing for comp. speed

temp = datestr(now,2);
temp(strfind(temp,'/')) =[]; % Today date
animal_id                                   = 'STV003';%'STV005';%'ARN088';
behavior                                    = 'Anesthesia';
result_dir                                  = ['results/' animal_id '/mo_active_learning'];

param_space                                 = 'Standard';
model_type                                  = '';
query_type                                  = '';
% reg_type: 1-ridge, 2-lasso, 3-Decision tree regression 
%query_type:1-EMCM, 2-QBC, 3- GS, 4-RD, 5-RD+EMCM, 6-RD+QBC, 7-RD+GS

experiment_name                             = sprintf('%s_GP_Uncertainty_Sampling_%s_%s',animal_id,param_space,behavior);

if strcmp(animal_id, '') ||  strcmp(experiment_name, '')
    [animal_id, experiment_name]                = get_experiment_info(DEBUG);    
end

% if DEBUG
%     experiment_name = [experiment_name 'DEBUG'];
% end

time_str                                    = datestr(now, 30);
log_pattern                                 = [result_dir '/' experiment_name '-%s_%s'];
exp_directory                               = sprintf(log_pattern, animal_id, time_str);
mkdir(exp_directory);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Define stim parameters and stimulation_object 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %parameters common across different spaces
load('Opto_param.mat');

V_bounds = [0.4 1.45];%[0.975 2.4375];%[0.4 1.45];%[0.975 2.4375];%[0.4 1.4];
fprintf('Voltage Boundary: %.2f - %.2f V\n',V_bounds(1),V_bounds(2));

amp_temp = opto_param(:,1);
amp_temp = (amp_temp - min(amp_temp))/range(amp_temp);
opto_param(:,1) = amp_temp*(V_bounds(2)-V_bounds(1))+V_bounds(1);

optimization_object             = opto_mo_al_controller();
optimization_object.param_space  = opto_param;

switch param_space
    case 'Standard'
        stimulation_manager                         = opto_stimulator_3();
        stimulation_manager.stimulation_type        = 'standard';
        
        stimulation_manager.stimulation_pulse_width_pulse   = 0;
        stimulation_manager.stimulation_pulse_width_train   = 0;
        
        optim_param = {'amplitude','pulse_frequency','pulse_width'};
        
    case 'Poisson'
        stimulation_manager                         = opto_stimulator_poisson();
        stimulation_manager.stimulation_type        = 'poisson';
        stimulation_manager.stimulation_refractory  = 1*10^-3;
        
        lower_bound         = [V_bounds(1) 5  2*10^-3];
        upper_bound         = [V_bounds(2) 42 10*10^-3];
        
        [input_cell,input_space] = create_input_arrays(lower_bound, upper_bound, 50);
        
        optim_param = {'amplitude','pulse_frequency','pulse_width'};
        
    case 'Sine'
        stimulation_manager                         = opto_stimulator_3_sine();
        stimulation_manager.stimulation_type        = 'sinusoid';
        stimulation_manager.stimulation_N_sinusoid  = 1;
        
        lower_bound         = [V_bounds(1) 2];
        upper_bound         = [V_bounds(2) 50];
        
        [input_cell,input_space] = create_input_arrays(lower_bound, upper_bound, 50);

        optim_param = {'amplitude','sine_frequency'};
        
    case 'Sine_2N'
        stimulation_manager                         = opto_stimulator_3_sine();
        stimulation_manager.stimulation_type        = 'sinusoid';
        stimulation_manager.stimulation_N_sinusoid  = 2;
        
        lower_bound         = [V_bounds(1)  2   2   0];
        upper_bound         = [V_bounds(2)  50  50  1];
        
        [input_cell,input_space] = create_input_arrays(lower_bound, upper_bound, 50);

        optim_param = {'amplitude','sine_frequency','sine_frequency_2','amplitude_ratio'};
        
    case 'NPT'
        stimulation_manager                         = opto_stimulator_3();
        stimulation_manager.stimulation_type        = 'train';
        
        lower_bound         = [V_bounds(1)   2    25    30   2*10^-3];
        upper_bound         = [V_bounds(2)   12   100   90   10*10^-3];
        
        [input_cell,input_space] = create_input_arrays(lower_bound, upper_bound, 50);

        optim_param = {'amplitude','train_frequency','pulse_frequency','train_width','pulse_width'};
        
    otherwise
        error('Error: Did not select a valid parameter space')
        
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

stimulation_manager.TD                      = TD;
stimulation_manager.device_name             = TD.GetDeviceName(0);
stimulation_manager.sampling_frequency      = TD.GetDeviceSF(stimulation_manager.device_name);
stimulation_manager.stimulation_channels    = [10];
            
stimulation_manager.headstage_type          = 'ZC16-OB1-16CH';
stimulation_manager.electrode_location      = 'R_HPC';
stimulation_manager.logging_directory       = exp_directory;
stimulation_manager.tank_name               = strcat('C:\TDT\OpenEx\MyProjects\CustomStimActiveX_Opto\DataTanks\CustomStimActiveX_Opto_DT1_',temp);
%stimulation_manager.block_name              = get_block_name(stimulation_manager.tank_name);

stimulation_manager.animal_id               = animal_id;
stimulation_manager.experiment_name         = experiment_name;
stimulation_manager.experiment_start_time   = posixtime(datetime('now'));
stimulation_manager.display_log_output      = 0;
stimulation_manager.initialize();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameter definition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
model_directory     = ''; %Only for PriSM optimization
model_name          = ''; %Only for PriSM optimization

recording_channels  = [1:16]; % Channel 13 is broken only on zif-clip headstage
sampling_frequency                      = stimulation_manager.sampling_frequency;
metric_type                             = 'PSD';
metric_def                              = [33 50];
metric                                  = opto_model_state(sampling_frequency, recording_channels, metric_def, metric_type, [model_directory model_name]);
metric_objects                          = {metric};
delta_or_raw                            = 'raw';
max_or_min                              = 'maximize';  %CHANGE

target_metric       = NaN;%%


%%%%%%%% use ferrule info to find V from chosen intensity values %%%%%%%%%
% ferrule_info                        = [4.468 -16.838];  %[slope, intercept] of calibration, with V and mW units
% ferrule_r                          = 0.20/2; %inner radius of ferrule, mm 
% I_bounds                            = [10 50];
% V_bounds = (I_bounds*pi*(ferrule_r)^2 - ferrule_info(2))./ferrule_info(1);
% V_bounds = [3.8 3.98];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
% Configure Bayesian optimization object
%

optimization_object.objective_function      = metric; 
optimization_object.objective_type          = delta_or_raw;
optimization_object.target_metric           = target_metric;
optimization_object.optimization_direction  = max_or_min;

optimization_object.model_type              = model_type;
optimization_object.query_type              = query_type;
optimization_object.TD                      = TD;
optimization_object.device_name             = TD.GetDeviceName(0);
optimization_object.sampling_frequency      = TD.GetDeviceSF(optimization_object.device_name);
optimization_object.stimulation_parameter   = optim_param;


optimization_object.n_parameters            = numel(optimization_object.stimulation_parameter);
optimization_object.logging_directory       = exp_directory;


%%%%%%%%%%%%%%%%%%%%%%%%%%%

optimization_object.stimulation_time_s      = 5;
stimulation_manager.stimulation_duration    = 5;
optimization_object.objective_window_s      = stimulation_manager.stimulation_duration;
optimization_object.evaluate_delay_s        = 5;

if exist('TURBO', 'var')
    optimization_object.stimulation_time_s      = 0.01;
    stimulation_manager.stimulation_duration    = .3;
    optimization_object.objective_window_s      = stimulation_manager.stimulation_duration;
    optimization_object.evaluate_delay_s        = 0.01;
end
optimization_object.stimulator              = stimulation_manager;

optimization_object.n_burn_in               = 3; 
optimization_object.n_trials                = 75; 

optimization_object.initialize();
end

function [param_cell,param_space] = create_input_arrays(lower_bounds, upper_bounds, resolution)
    
    param_cell = cell(length(lower_bounds),1);
    for k = 1:length(lower_bounds)
        param_cell{k} = linspace(lower_bounds(k),upper_bounds(k),resolution);
    end
    param_space = transpose(combvec(param_cell{:}));
end

