function stim_signal = opto_generate_nested_train(f_sample, f_stim_T,f_stim_P , duration, amplitude, PWT, PWP, n_channels)
%Updated version of opto_generate_train.m: Given a set of stimulation
%parameters for nested pulse train, generates a time series for the stimulation pattern
%to be delivered. Updated version handles irregular edge cases based on
%stimulation parameters and returns a flag of -1 when they are selected (uncomment 
%portions where flags are returned to default to non-train stimulus instead):
%
% 1) If train duration is too long, stimulation is no longer a nested train
% - just continuous high frequency pulsing
% 2) Uneven pulse width+freq./train duration ratio: truncate to whole number
% of pulses
% 3) High pulse width -> train just becomes a continuous big pulse
%
% Inputs
%   f_sample    - sampling rate
%   f_stim_T    - train-frequency
%   f_stim_P    - pulse-frequency
%   duration    - duration of nested-pulse-train stimulation
%   amplitude   - pulse amplitude
%   PWT         - pulse-width for the train
%   PWP         - pulse-width for each pulse
%   n_channels  - replicates signal


if duration < 10.1
    t   = 1/f_sample:1/f_sample:duration ;
    dP  = PWP/2:1/f_stim_P:duration;
    
    if PWP == 1
        yP = ones(1,length(t));
    else
        yP = pulstran(t,dP,'rectpuls',PWP);
    end
    
    if PWT == 1
        yT = amplitude*ones(1,length(yP));
    else
        dT = PWT/2:1/f_stim_T:duration+ PWT/2;
        yT = amplitude*pulstran(t,dT,'rectpuls',PWT);
    end
    
    if (PWP >= 1/f_stim_P) && (PWT > (1/f_stim_T))
        stim_signal = -1;
        return
        %disp('Warning: both pulse and train width too high. Now relentlessly cooking brain with continuous stimulus.');
        %stim_signal = ones(size(yT))*amplitude;
    elseif PWP >= 1/f_stim_P
        stim_signal = -1;
        return
        %disp('Warning: pulse width too high. Now using train of long pulses instead of nested pulse train.');
        %stim_signal = yT;
    elseif PWT > (1/f_stim_T)
        stim_signal = -1;
        return
        %disp('Warning: train width too high. Now using train of short pulses instead of nested pulse train.');
        %stim_signal     = amplitude*yP;
    else
        %Edge Case: predict need to truncate/extend cut-off pulse.
        %Detect whether last pulse in a train gets cut off, and round the PWT for a
        %whole number of pulses if so.
        num_pulses = PWT*f_stim_P;
        t_pstart = floor(num_pulses)/f_stim_P;
        t_pend = t_pstart + PWP;
        
        if (t_pend>PWT) && (t_pstart<PWT) && (PWP<1/f_stim_P)
            %disp('Warning: Truncating for cut-off pulse')
            stim_signal = -1;
            return
            
%             PWT = round(num_pulses)/f_stim_P;
%             dT = PWT/2:1/f_stim_T:duration+ PWT/2;
%             yT = amplitude*pulstran(t,dT,'rectpuls',PWT);
        end
        stim_signal = yP.*yT;
    end
    
    stim_signal     = repmat(stim_signal,n_channels,1);
    
else % if duration>10
    temp_duration   = 5;
    N_rep           = floor(duration/temp_duration);
    Sec_remain      = round(duration - N_rep*temp_duration);
    
    t               = 0:1/f_sample:temp_duration ;
    dP              = PWP/2:1/f_stim_P:temp_duration;
    
    if PWP == 1
        yP  = ones(1,length(t));
    else
        yP  = pulstran(t,dP,'rectpuls',PWP);
    end
    
    if PWT == 1
        yT = amplitude*ones(1,length(yP));
    else
        dT  = PWT/2:1/f_stim_T:temp_duration+ PWT/2;
        yT  = amplitude*pulstran(t,dT,'rectpuls',PWT);
    end
    
    if (PWP >= 1/f_stim_P) && (PWT > (1/f_stim_T))
        %         disp('Warning: both pulse and train width too high. Now relentlessly cooking brain with continuous stimulus.');
        %         stim_signal_temp = ones(size(yT))*amplitude;
        stim_signal = -1;
        return
    elseif PWP >= 1/f_stim_P
        %         disp('Warning: pulse width too high. Now using train of long pulses instead of nested pulse train.');
        %         stim_signal_temp = yT;
        stim_signal = -1;
        return
    elseif PWT > (1/f_stim_T)
        %         disp('Warning: train width too high. Now using train of short pulses instead of nested pulse train.');
        %         stim_signal_temp     = amplitude*yP;
        stim_signal = -1;
        return
    else
        %Edge Case: predict need to truncate/extend cut-off pulse.
        %Detect whether last pulse in a train gets cut off, and round the PWT for a
        %whole number of pulses if so.
        num_pulses = PWT*f_stim_P;
        t_pstart = floor(num_pulses)/f_stim_P;
        t_pend = t_pstart + PWP;
        
        if (t_pend>PWT) && (t_pstart<PWT) && (PWP<1/f_stim_P)
            stim_signal = -1;
            return
            %             disp('Warning: Truncating for cut-off pulse')
%             PWT = round(num_pulses)/f_stim_P;
%             dT = PWT/2:1/f_stim_T:duration+ PWT/2;
%             yT = amplitude*pulstran(t,dT,'rectpuls',PWT);
        end
        stim_signal_temp = yP.*yT;
    end
    
    stim_signal_temp    = repmat(stim_signal_temp,n_channels,N_rep);
    
    temp_duration       = Sec_remain;
    t                   = 0:1/f_sample:temp_duration ;
    dP                  = 0:1/f_stim_P:temp_duration;
    
    if PWP == 1
        yP = ones(1,length(t));
    else
        yP = pulstran(t,dP,'rectpuls',PWP);
    end
    
    if PWT == 1
        yT = amplitude*ones(1,length(yP));
    else
        dT = PWT/2:1/f_stim_T:temp_duration+ PWT/2;
        yT = amplitude*pulstran(t,dT,'rectpuls',PWT);
    end
    
    stim_signal = [stim_signal_temp repmat(yT.*yP,n_channels,1)];
    
    
end
end
