classdef opto_mo_al_controller < handle
    %Adapted version of opto_bayesian_optimization_controller to run
    %Gaussian process uncertainty sampling data collection; version 2
    %extends code for multiple parameter spaces and adds hill climbing for
    %computational efficiency
    
    properties
        TD
        device_name
        sampling_frequency
        logging_directory
        
        % Timing parameters
        run_start_time_s
        cycle_start_time_s
        this_time_s
        last_time_s
        
        stimulation_time_s
        evaluate_delay_s
        objective_function
        objective_window_s
        objective_type
        optimization_direction
        sample_skip
        sample_results
        
        % Optimization parameters
        control_time
        n_samples
        n_trials
        next_sample_set
        n_burn_in
        parameter_vector
        pre_stimulation_metric
        state_threshold
        
        %Mo AL parameters
        init_stim_param
        sample_idx
        model_type
        query_type
        exp_type
        
        param_space
        stimulation_parameter
        n_parameters
        target_metric
        
        % Acquisition parameters
        lower_bound
        upper_bound
        
        % Logging parameters
        objective_fid   
        parameter_fid
        state_fid
        
        optimization_done
        stimulator
    end
    
    methods
      
        %%%%%%%%%%%%%%%
        %
        %%%%%%%%%%%%%%%
        function initialize(obj, gp_model)
            
            
            obj.run_start_time_s        = posixtime(datetime);

            obj.optimization_done       = 0;

            obj.objective_fid           = fopen(sprintf('%s/objective_function.csv',...
                                            obj.logging_directory),'a');
            obj.parameter_fid           = fopen(sprintf('%s/parameter.csv',...
                                            obj.logging_directory),'a');
            obj.state_fid               = fopen(sprintf('%s/state_function.csv',...
                obj.logging_directory),'a');
           
            obj.n_samples               = 1;
            
            if nargin == 2
                 obj.n_samples = size(gp_model.y_data,1);
                 obj.parameter_vector = gp_model.x_data;
                 obj.sample_results = gp_model.y_data;
                 obj.gp_model = gp_model;
            end
            
            if strcmp(obj.exp_type,'AL') || strcmp(obj.exp_type,'Random')
                [obj.init_stim_param,obj.sample_idx] = active_learning_init(obj.param_space);
                obj.n_burn_in = size(obj.init_stim_param,1);
                obj.sample_results  = nan(obj.n_trials,1);
            else
                obj.parameter_vector = obj.param_space(randperm(size(obj.param_space,1),obj.n_trials),:);
            end
            
            obj.reset_time();
        end
        
        %%%%%%%%%%%%%%%
        %
        %%%%%%%%%%%%%%%
        function optimize(obj)
            
%             obj.stimulator.check_stimulation()
            
            this_time  = obj.get_current_time();

            if this_time > ...
                obj.stimulation_time_s +  obj.cycle_start_time_s && ...
                obj.stimulator.stimulation_armed
            
                obj.select_next_sample();

                %%%%%%%%%%%%%%%%%%%%
                %                 obj.TD.SetTargetVal([obj.stimulator.device_name '.nPeriod'],obj.stimulator.sampling_frequency/obj.stimulator.stimulation_frequency);
                %                 obj.TD.SetTargetVal([obj.stimulator.device_name '.nPulses'],obj.stimulator.stimulation_frequency * obj.stimulator.stimulation_duration);
                %                 % optimizer.width                     = 2/stimulation_manager.sampling_frequency;
                %                 obj.TD.SetTargetVal([obj.stimulator.device_name '.nDur-A'],obj.stimulator.sampling_frequency/1000*obj.stimulator.stimulation_pulse_width);
                %                 obj.TD.SetTargetVal([obj.stimulator.device_name '.nDur-B'],5); % Useless
                %                 obj.TD.SetTargetVal([obj.stimulator.device_name '.Amp-A'],obj.stimulator.stimulation_amplitude);
                obj.stimulator.generate_stimulation_signal();
                obj.stimulator.write_stimulation_signals();
                %%%%%%%%%%%%%%%%%%%%
                
                obj.stimulator.stimulate()  
                
            elseif this_time  > ...
                    obj.evaluate_delay_s + obj.stimulator.stimulation_duration + obj.stimulator.stimulation_time && ...
                    ~obj.stimulator.stimulation_armed
                
                obj.evaluate_objective_function(); 
                
                obj.reset_time();
                obj.stimulator.stimulation_armed = 1;
            end  
            if obj.n_samples == obj.n_trials+1
                obj.optimization_done = 1;
                al_results.param_space = obj.param_space;
                al_results.reg_type = obj.model_type;
                al_results.query_type = obj.query_type;
                al_results.params = obj.parameter_vector;
                al_results.obj_samples = obj.sample_results;
                save(sprintf('%s/al_results.mat',obj.logging_directory),'al_results');
            end
        end
        
        function select_next_sample(obj)
            if strcmp(obj.exp_type,'AL') || strcmp(obj.exp_type,'Random')
                if obj.n_samples <= obj.n_burn_in
                    obj.parameter_vector(obj.n_samples,:) = obj.init_stim_param(obj.n_samples,:);
                    
                elseif strcmp(obj.exp_type,'AL')
                    param_temp = active_learning_function(obj.param_space,...
                        obj.parameter_vector(1:obj.n_samples-1,:),obj.sample_results(1:obj.n_samples-1),...
                        obj.model_type,obj.query_type);
                    obj.parameter_vector = [obj.parameter_vector; param_temp];
                    
                elseif strcmp(obj.exp_type,'Random')
                    param_temp = obj.param_space(randperm(size(obj.param_space,1),1),:);
                    obj.parameter_vector = [obj.parameter_vector; param_temp];
                end
                
                
            end
            
            fprintf('Sample: %d\n', obj.n_samples);
            for c1 = 1:obj.n_parameters
                % Update the stimulation parameters
                
                fprintf('\tStimulation %s = %e\n', obj.stimulation_parameter{c1}, obj.parameter_vector(obj.n_samples, c1));
                

                switch obj.stimulation_parameter{c1}
                    case 'train_frequency'
                        obj.stimulator.stimulation_frequency_train      = obj.parameter_vector(obj.n_samples, c1);
                    case 'pulse_frequency'
                        obj.stimulator.stimulation_frequency_pulse      = obj.parameter_vector(obj.n_samples, c1);
                    case 'duration'
                        obj.stimulator.stimulation_duration             = obj.parameter_vector(obj.n_samples, c1);
                    case 'amplitude'
                        obj.stimulator.stimulation_amplitude            = obj.parameter_vector(obj.n_samples, c1);
                    case 'state_threshold'
                        obj.state_threshold                             = obj.parameter_vector(obj.n_samples, c1);
                    case 'train_width'
                        obj.stimulator.stimulation_pulse_width_train    = obj.parameter_vector(obj.n_samples, c1);
                    case 'pulse_width'
                        obj.stimulator.stimulation_pulse_width_pulse    = obj.parameter_vector(obj.n_samples, c1);
                    case 'sine_frequency'
                        obj.stimulator.stimulation_frequency            = obj.parameter_vector(obj.n_samples, c1);
                    case 'sine_frequency_2'
                        obj.stimulator.stimulation_frequency            = [obj.stimulator.stimulation_frequency obj.parameter_vector(obj.n_samples, c1)];
                    case 'amplitude_ratio'
                        obj.stimulator.stimulation_amplitude            = [obj.stimulator.stimulation_amplitude obj.parameter_vector(obj.n_samples, c1)];
                end
                
            end
            
            s = sprintf('%e,',obj.parameter_vector(obj.n_samples,:));
            fprintf(obj.parameter_fid,  [s(1:end-1) '\n']);
        end
        
        %%%%%%%%%%%%%%%
        % Resets the clock to start another stimulation trial
        %%%%%%%%%%%%%%%
        function reset_time(obj)
            obj.cycle_start_time_s   = obj.get_current_time();
            obj.last_time_s          = obj.get_current_time();
            obj.this_time_s          = obj.get_current_time();
        end

        %%%%%%%%%%%%%%%
        % 
        %  
        %%%%%%%%%%%%%%%
        function t = get_current_time(obj)
            t = obj.TD.ReadTargetVEX([obj.device_name '.current_time'], 0, 1, 'I32', 'F64')/obj.sampling_frequency;
        end
        
         %%%%%%%%%%%%%%%
        % Evaluate objective function - this needs to be split into it's
        % own object 
        %%%%%%%%%%%%%%%
        function evaluate_objective_function(obj)
            objective_window_prestart   = obj.stimulator.stimulation_time - obj.objective_window_s;
            objective_window_start      = obj.stimulator.stimulation_time;
            objective_window_end        = objective_window_start + obj.objective_window_s;
            if(objective_window_prestart < 0)
                a = 1;
            end
                
            obj.pre_stimulation_metric  = obj.objective_function.get_metric(objective_window_prestart, objective_window_start,0);
            
            fprintf('Pre-stimulation state = %e\n', obj.pre_stimulation_metric);
            m = obj.objective_function.get_metric(objective_window_start, objective_window_end,1);
            
            if ~strcmp(obj.objective_type,'threshold') || obj.pre_stimulation_metric<0.1
                if ~isnan(m) && ~isnan(obj.pre_stimulation_metric)
                    switch(obj.objective_type)
                        case 'raw'
                            if strcmp(obj.optimization_direction,'maximize')
                                sign = 1;
                            elseif strcmp(obj.optimization_direction,'minimize')
                                sign = -1;
                            end
                            obj.sample_results(obj.n_samples,1) = sign*m;
                            
                        case 'delta'
                            if strcmp(obj.optimization_direction,'maximize')
                                sign = 1;
                            elseif strcmp(obj.optimization_direction,'minimize')
                                sign = -1;
                            end
                            obj.sample_results(obj.n_samples,1) = sign*(m-obj.pre_stimulation_metric);
                            
                       case 'complex'
                            if strcmp(obj.optimization_direction,'maximize')
                                sign = 1;
                            elseif strcmp(obj.optimization_direction,'minimize')
                                sign = -1;
                            end
                            obj.sample_results(obj.n_samples,1) = sign*((m-obj.pre_stimulation_metric+1)*(1+obj.pre_stimulation_metric));
                            
                        case 'target'
                            if strcmp(obj.optimization_direction,'maximize')
                                sign = 1;
                            elseif strcmp(obj.optimization_direction,'minimize')
                                sign = -1;
                            end
                            obj.sample_results(obj.n_samples,1) = sign*(abs(m-obj.target_metric));

                    end
                    
                    fprintf(obj.objective_fid,'%f, %f, %e\n', ...
                        obj.stimulator.stimulation_time, posixtime(datetime('now')), obj.sample_results(obj.n_samples,1));
                    
                    fprintf(obj.state_fid,'%f, %f, %e\n', ...
                        obj.stimulator.stimulation_time, posixtime(datetime('now')), obj.pre_stimulation_metric);
                    fprintf('Objective: %.5e\n', obj.sample_results(obj.n_samples,1));
                    
                    obj.n_samples = obj.n_samples+1;
                else
                    %                 obj.parameter_vector(obj.n_samples,:) = [];
                end
            end
        end
        
    end
    

end

