function stim_signal = opto_generate_sinusoid_stim(f_sample, duration, N_sin, amp, fr)
%opto_generate_sinusoid.m: Generates sinusoidal stimulation signal given
%amplitude and frequency parameters
%
%f_sample: sampling rate
%duration: length of signal in seconds
%N_sin: number of summed sinusoids of signal (1 = default sine wave, 2 =
%sum of two sinusoids)
%amp/fr: vector of N_sin x 1 of amplitude or frequency parameters for each
%summed sinusoid

tt = 1/f_sample:1/f_sample:duration;
stim_signal = zeros(1,length(tt));
if N_sin==1
%     stim_signal = stim_signal + max(amp);
%     stim_signal = stim_signal + amp(1)*sin(2*pi*fr(1)*tt);
stim_signal = amp(1)*sin(2*pi*fr(1)*tt);
stim_signal = (stim_signal + max(stim_signal))/2;

    
elseif N_sin == 2
    %fr = 2D vector
    %amp = main amplitude, amp. ratio
    stim_signal = stim_signal + (1-amp(2))*sin(2*pi*fr(1)*tt);
    stim_signal = stim_signal + amp(2)*sin(2*pi*fr(2)*tt);
    %ratio = amp(1)/max(stim_signal+min(stim_signal));
    stim_signal = (stim_signal-min(stim_signal));
    stim_signal = stim_signal*(amp(1)/max(stim_signal));
end



end



















