classdef opto_stimulator_3 < handle
    %STIMULATOR Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        TD
        device_name
        sampling_frequency
        headstage_type
        electrode_location
        animal_id
        experiment_name
        tank_name
        block_name
        experiment_start_time
                
        stimulation_type
        stimulation_frequency
        stimulation_frequency_train
        stimulation_frequency_pulse
        stimulation_duration
        stimulation_pulse_width_pulse
        stimulation_pulse_width_train
        stimulation_amplitude
        stimulation_channels
        
        
        stimulation_signal
        stimulation_time
        stimulation_armed
        channels_closed
        
        logging_directory
        display_log_output
        stimulation_fid
        stimulation_uid
    end
    
    methods
        
        %%%%%%%%%%%%%%%
        %
        %%%%%%%%%%%%%%%
        function initialize(obj)
            obj.stimulation_armed       = 1;
            obj.channels_closed         = 0;
            obj.stimulation_fid         = fopen(sprintf('%s/stimulation_table.csv',...
                                            obj.logging_directory), 'a'); 
            obj.stimulation_time        = 1;        
            obj.write_log_header();
            obj.open_channels();
        end
        
        %%%%%%%%%%%%%%%
        %
        %%%%%%%%%%%%%%%
        function generate_stimulation_signal(obj)
            % Determine the type of stimulation signal
            switch obj.stimulation_type, 
                case 'synchronous'
                
                    signal = generate_biphasic(obj.sampling_frequency,  ...
                        obj.stimulation_frequency,                      ...
                        obj.stimulation_duration,                       ...
                        obj.stimulation_amplitude,                      ...
                        obj.stimulation_pulse_width,                    ...
                        size(obj.stimulation_channels,2));
                
                
                case 'bipolar'
                
                    signal = generate_bipolar(obj.sampling_frequency,       ...
                        obj.stimulation_frequency,                          ...
                        obj.stimulation_duration,                           ...
                        obj.stimulation_amplitude,                          ...
                        obj.stimulation_pulse_width,                        ...
                        size(obj.stimulation_channels,1));
                    
                case 'train'
                    signal = opto_generate_train(obj.sampling_frequency,       ...
                        obj.stimulation_frequency_train,                       ...
                        obj.stimulation_frequency_pulse,                    ...
                        obj.stimulation_duration,                           ...
                        obj.stimulation_amplitude,                          ...
                        obj.stimulation_pulse_width_train,                        ...
                        obj.stimulation_pulse_width_pulse,                        ...
                        size(obj.stimulation_channels,1));
                case 'calibration'
                    signal = opto_generate_calibration(obj.sampling_frequency,       ...
                        obj.stimulation_duration,                           ...
                        obj.stimulation_amplitude,                          ...
                        size(obj.stimulation_channels,1));
                case 'standard'
                    signal = opto_generate_standard(obj.sampling_frequency,                     ...
                        obj.stimulation_frequency_pulse,                    ...
                        obj.stimulation_duration,                           ...
                        obj.stimulation_amplitude,                        ...
                        obj.stimulation_pulse_width_pulse,                        ...
                        size(obj.stimulation_channels,1));
                case 'arbitrary'
                    signal = opto_generate_arbitrary(obj.sampling_frequency,                     ...
                        obj.stimulation_frequency_pulse,                    ...
                        obj.stimulation_duration,                           ...
                        obj.stimulation_amplitude,                        ...
                        obj.stimulation_pulse_width_pulse,                        ...
                        size(obj.stimulation_channels,1));
             end
            
            obj.stimulation_signal = signal;
        end
        
        %%%%%%%%%%%%%%%
        %
        %%%%%%%%%%%%%%%
        function open_channels(obj)
            for c1 = 1:size(obj.stimulation_channels,2);
                obj.TD.WriteTargetVEX( ...
                    [obj.device_name '.stim_buff~' num2str(obj.stimulation_channels(c1))], ...
                    0, 'F32',  10000*ones(1,10));
                
%                 obj.TD.SetTargetVal([obj.device_name '.dur~' num2str(obj.stimulation_channels(c1))], 10);
                obj.TD.SetTargetVal([obj.device_name '.nPeriod'],24000);
                obj.TD.SetTargetVal([obj.device_name '.nPulses'],10);
                obj.TD.SetTargetVal([obj.device_name '.nDur-A'],244);
                obj.TD.SetTargetVal([obj.device_name '.nDur-B'],5);
                obj.TD.SetTargetVal([obj.device_name '.Amp-A'],3);
            end
        end
        
        %%%%%%%%%%%%%%%
        %
        %%%%%%%%%%%%%%%
        function v = is_stimulating(obj)
            v = 0;
            for c1 = 1:numel(obj.stimulation_channels)
                v = v + obj.TD.ReadTargetVEX(...
                    [obj.device_name '.not_done'],...
                    0, 1, 'F32', 'F32');
            end
            
        end
        
        %%%%%%%%%%%%%%%
        %
        %%%%%%%%%%%%%%%
        function write_stimulation_signals(obj)
            obj.stimulation_uid     = posixtime(datetime('now'));
%             for c1 = 1:length(obj.stimulation_channels)

%                 obj.TD.WriteTargetVEX(...
%                     [obj.device_name '.stim_buff~' num2str(obj.stimulation_channels(c1))],...
%                     0, 'F32',   [0 obj.stimulation_signal(c1,:) ]);
%                 
%                 obj.TD.SetTargetVal(...
%                     [obj.device_name '.dur~' num2str(obj.stimulation_channels(c1))], ...
%                     1+numel(obj.stimulation_signal(c1,:)));

            obj.TD.WriteTargetVEX(...
                [obj.device_name '.stim_buff'],...
                0, 'F32',   [0 obj.stimulation_signal ]);

            obj.TD.SetTargetVal(...
                [obj.device_name '.dur'], ...
                1+numel(obj.stimulation_signal));
            %             end
        end
        
        %%%%%%%%%%%%%%%
        %
        %%%%%%%%%%%%%%%
        function stimulate(obj)       
            obj.TD.SetTargetVal([obj.device_name '.trig_stim'], 1);

            pause(.1)
            obj.TD.SetTargetVal([obj.device_name '.trig_stim'], 0);
            
            obj.channels_closed     = 1;
            obj.stimulation_armed   = 0;
            obj.stimulation_time    = obj.get_last_stimulation_time();
            obj.log_stimulation();
        end

        %%%%%%%%%%%%%%%
        %
        %%%%%%%%%%%%%%%
        function check_stimulation(obj)
            if  obj.channels_closed && ~obj.is_stimulating()
                obj.open_channels();
                obj.channels_closed     = 0;
            end
        end

        %%%%%%%%%%%%%%%
        %
        %%%%%%%%%%%%%%%
        function log_stimulation(obj)
            log_string = '';
            log_string = [log_string sprintf('%f,', obj.stimulation_uid)];
            log_string = [log_string sprintf('%s,', obj.headstage_type)];
            log_string = [log_string sprintf('%s,', obj.electrode_location)];
            log_string = [log_string sprintf('%f,', obj.experiment_start_time)];
            log_string = [log_string sprintf('%s,', obj.experiment_name)];
            log_string = [log_string sprintf('%s,', obj.animal_id)];
            log_string = [log_string sprintf('%s,', obj.block_name)];
            log_string = [log_string sprintf('%s,', obj.sampling_frequency)];
            log_string = [log_string sprintf('%f,', obj.stimulation_time)];
            log_string = [log_string sprintf('%f,', obj.stimulation_frequency_train)];
%             if strcmp(obj.stimulation_type,'train')
                log_string = [log_string sprintf('%f,', obj.stimulation_frequency_pulse)];
%             else
%                 log_string = [log_string sprintf('%f,', 0.0)];
%             end
            log_string = [log_string sprintf('%f,', obj.stimulation_duration)];
            log_string = [log_string sprintf('%f,', obj.stimulation_amplitude/2)];
            log_string = [log_string sprintf('%f,', obj.stimulation_amplitude/2)];
            log_string = [log_string sprintf('%f,', obj.stimulation_pulse_width_train/2)];
            log_string = [log_string sprintf('%f,', obj.stimulation_pulse_width_train/2)];
%             if strcmp(obj.stimulation_type,'train')
                log_string = [log_string sprintf('%f,', obj.stimulation_pulse_width_pulse/2)];
                log_string = [log_string sprintf('%f,', obj.stimulation_pulse_width_pulse/2)];
%             else
%                 log_string = [log_string sprintf('%f,', 0.0)];clea
%                 log_string = [log_string sprintf('%f,', 0.0)];
%             end
%             if strcmp(obj.stimulation_type,'train')
                log_string = [log_string sprintf('%f,', 1.0)];
%             else
%                 log_string = [log_string sprintf('%f,', 0.0)];
%             end
            
            log_string = [log_string sprintf('%s,', obj.stimulation_type)];
            log_string = [log_string sprintf('10')];
%             log_string = [log_string sprintf('%d', obj.stimulation_channels)];
            log_string = [log_string sprintf('\n')];
            
            switch obj.display_log_output
                case 'verbose'
                    fprintf(log_string);
                case 'simple'
                    fprintf('Frequency_pulse: %.1f, Frequency_train: %.1f Duration: %.2f, Amplitude: %.2f\n', ...
                        obj.stimulation_frequency_pulse,obj.stimulation_frequency_train, obj.stimulation_duration,...
                        obj.stimulation_amplitude);
            end
            fprintf(obj.stimulation_fid, log_string);
        end
        
        %%%%%%%%%%%%%%%
        %
        %%%%%%%%%%%%%%%
        function write_log_header(obj)
            header_string = ['stimulation_uid,' ...
                'headstage_type,' ...
                'electrode_location,' ...
                'experiment_start_time,' ...
                'experiment_name,' ...
                'animal_id,' ...
                'block_name,' ...
                'sampling_frequency,' ...
                'stimulation_time,' ...
                'pulse_frequency_train,' ...
                'pulse_frequency_pulse,' ...
                'stimulation_duration,' ...
                'stimulation_amplitude_a,' ...
                'stimulation_amplitude_b,' ...
                'stimulation_pulse_width_train_a,' ...
                'stimulation_pulse_width_train_b,' ...
                'stimulation_pulse_width_pulse_a,' ...
                'stimulation_pulse_width_pulse_b,' ...
                'stimulation_train_TorF,' ...
                'stimulation_type,' ...
                'stimulation_channel_order\n' ...
                ];
            fprintf(obj.stimulation_fid, header_string);

        end
        
        %%%%%%%%%%%%%%%
        %
        %%%%%%%%%%%%%%%
        function t = get_last_stimulation_time(obj)
            t = obj.TD.ReadTargetVEX([obj.device_name '.last_stim_time'], 0, 1, 'I32', 'F64')/obj.sampling_frequency;
        end
    end
    
end

