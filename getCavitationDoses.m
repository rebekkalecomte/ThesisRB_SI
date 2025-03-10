%% Define functions
% Fourier on 1 time domain
function [freqaxis, absolutefft, amplitude_ss] = fourier(timedata, voltdata, frequency_Hz, cycles, plt)
    %define relevant part of the signal
    start_idx = find(timedata>= -2e-6, 1);                                                  % Signal starts 2 microseconds earlier than expected
    skipcycles = 5;                                                                 
    start_time_real = timedata(start_idx);
    start_time = start_time_real + (skipcycles / frequency_Hz);                             % Eleminate the 5-cylce ramping up of the gaussian modulated pulse
    end_time = start_time_real + (1/frequency_Hz) * cycles - (skipcycles / frequency_Hz);   % Eleminate the 5-cylce ramping down of the gaussian modulated pulse
    end_idx = find(timedata >= end_time, 1);
    if plt == 1
        figure  
        plot(timedata*1000000, voltdata*1000); %to us and milivolt
        xlabel('Time (us)'), ylabel("Amplitude (mV)"), title('Time domain')
        xline(start_time*1000000, 'r--');
        xline(end_time*1000000, 'r--');
    end
    
    %cut relevant time domain and define parameters
    t = timedata(start_idx:end_idx);
    v = voltdata(start_idx:end_idx);
    signallength = max(t) - min(t);                                                         % Acquisition time
    N = length(t);                                                                          % Number of samples per acquisition
    Fs = N/ signallength;                                                                   % Sampling frequency
    T = 1/ Fs;                                                                              % Time between 2 succesive samples
    Nzero = 10*N;                                                                           % Zero padding FFT
    dF = Fs / Nzero;                                                                        % Size of a frequency bin is inversely proportional to the total duration of signal 
    
    % Generate double sided frequency axes, centered around 0. Max determined by nyquist
    freqaxis = (Fs/Nzero)*linspace(-Nzero/2,Nzero/2-1,Nzero); 
    % Add windowing
    v_win = v.*hanning(length(v))';                          
    % Simple fast fourier transform, plotting ignores the imaginary parts of the complex result so amplitude and phase information is not visibile
    plainfft = fft(v_win, Nzero);   
    % FFT SHIFT: transfers the frequency information from the sides to the center.
    shiftedfft= fftshift(plainfft);
    % ABS AND SCALING: retrieve amplitude information from the complex nature of the signal: the length
    absolutefft = abs(shiftedfft)/(N/2);                                                    % Use N to keep amplitudes correct, not Nzero!
    % Save amplitude 
    amplitude_ss = 2*max(absolutefft); %in volts
    
    % Plot double sided fd
    if plt == 1
        figure
        plot(freqaxis/1000000, absolutefft*1000);                                           % Convert to MHz and mV
        xlabel('Frequency (MHz)'), ylabel("Amplitude (mV)"), title('Fourier domain')
    
        % Display amplitude_ss in red in the upper-right corner
        xlim_range = xlim;                                                                  % Get current x-axis limits
        ylim_range = ylim;                                                                  % Get current y-axis limits
        text_position_x = xlim_range(2) * 0.95;                                             % Position slightly left of the maximum x-limit
        text_position_y = ylim_range(2) * 0.95;                                             % Position slightly below the maximum y-limit
        text(text_position_x, text_position_y, sprintf('Amplitude_{ss}: %.3f V', amplitude_ss), ...
            'Color', 'r', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'right');
    end

end 


% Mean Fourier of the received therapeutic burst
function [freqaxis, mean_absolutefft, amplitudelist] = meanfourierburst(inputstruct, frequency_Hz, cycles, plt)
    
    % Preallocate storage for spectra and amplitudes
    numFields = length(fieldnames(inputstruct));                                            % Input struct the struct received from the GUI
    numMeasurements = numFields/6;
    all_absolutefft = [];
    amplitudelist = [zeros(numMeasurements, 1)];

    % Fourier transform all 30 pulses within the therapeutic burst
    for i= 1:numMeasurements
        timedata = inputstruct.(sprintf('Time_%d', i));
        voltdata = inputstruct.(sprintf('Volt_%d', i));
        [freqaxis, absolutefft, amplitude_ss] =  fourier(timedata, voltdata, frequency_Hz, cycles, 0);

        % Store max amplitude 
        amplitudelist(i, 1) = amplitude_ss;       

        % Initialize storage with the correct size
        if isempty(all_absolutefft)
            all_absolutefft = zeros(length(absolutefft), numMeasurements);
        end

        % Store spectrum
        all_absolutefft(:, i) = absolutefft; 
    end 

    % Compute mean absolute fft and mean maximal amplitude
    mean_absolutefft = mean(all_absolutefft, 2);                                            % Mean returns a column vector containing the mean of each row
    amplitude_ss = 2*max(mean_absolutefft);                                                 % In volts

   

    % Find the frequency of the maximum amplitude in the mean spectrum
    [~, max_idx] = max(mean_absolutefft);                                                   % Index of the max amplitude
    freq_max = freqaxis(max_idx);                                                           % Frequency corresponding to the max amplitude

    % Find the frequency that is twice the frequency of max amplitude
    freq_double = 2 * freq_max;
    [~, double_idx] = min(abs(freqaxis - freq_double));                                     % Closest index to 2*freq_max

    % Compute the SD at the max amplitude frequency and twice that frequency
    sd_at_max = std(all_absolutefft(max_idx, :)); 
    sd_at_double = std(all_absolutefft(double_idx, :)); 

    % Plot the mean Fourier spectrum
    if plt == 1
        figure
        plot(freqaxis / 1e6, mean_absolutefft * 1e3)
        xlabel('Frequency (MHz)'), ylabel("Amplitude (mV)"), title('Mean Fourier domain Rx')
        xlim([-10, 10]); % Set the x-axis limits to -10 MHz to 10 MHz

        % Display information in the upper-right corner of the plot
        xlim_range = xlim;                                                                  % Get current x-axis limits
        ylim_range = ylim;                                                                  % Get current y-axis limits
        text_position_x = xlim_range(2) * 0.95;                                             % Slightly left of the maximum x-limit
        text_position_y = ylim_range(2) * 0.90;                                             % Slightly below the maximum y-limit

        % Display amplitude_ss and SDs
        text(text_position_x, text_position_y, ...
            sprintf(['Amplitude_{ss}: %.3f mV\n', ...
                     'SD_{fund}: %.3f mV\n', ...
                     'SD_{harm}: %.3f mV'], ...
                amplitude_ss * 1e3, ...                                                     % Amplitude in mV
                sd_at_max * 1e3, ...                                                        % SD at fundamental frequency (mV)
                sd_at_double * 1e3), ...                                                    % SD at harmonic frequency (mV)
            'Color', 'r', 'FontSize', 9, 'FontWeight', 'bold', 'HorizontalAlignment', 'right');
    end

    
end



% Get list of amplitudes at harmonic frequency for all 30 pulses within therapeutic burst
function  harmamplitudelist = getHarmAmpl(inputstruct, frequency_Hz, cycles)
   
    % Preallocate storage for spectra and amplitudes
    numFields = length(fieldnames(inputstruct));
    numMeasurements = numFields/6;
    harmamplitudelist = [zeros(numMeasurements, 1)];

    % Fourier transform all 30 pulses within the therapeutic burst
    for i= 1:numMeasurements
        timedata = inputstruct.(sprintf('Time_%d', i));
        voltdata = inputstruct.(sprintf('Volt_%d', i));
        [freqaxis, absolutefft, ~] =  fourier(timedata, voltdata, frequency_Hz, cycles, 0);

       % Find the frequency of the maximum amplitude in the spectrum
       [~, max_idx] = max(absolutefft);                                                     % Index of the max amplitude
       freq_max = freqaxis(max_idx);                                                        % Frequency corresponding to the max amplitude
    
        % Find the 1st harmonic frequency
        freq_double = 2 * freq_max;
        [~, double_idx] = min(abs(freqaxis - freq_double));                                 % Closest index to 2*freq_max
        
       
        
        % Define the frequency search range (±0.15 MHz) for the highest peak, because the peak
        % is never exactly at the expected location
        bandwidth = 0.3e6; % 0.2 MHz in Hz
        freq_lower = freq_double - 0.5*bandwidth;
        freq_upper = freq_double + 0.5*bandwidth;


        % Find indices and extract the aboslutefft values within the search range
        valid_indices = find(freqaxis >= freq_lower & freqaxis <= freq_upper);
        absolutefft_band = absolutefft(valid_indices);
        
        % Find the highest peak within this range and get corresponding
        % index in the original absolute fft array
        [ss_amp, band_peak_idx] = max(absolutefft_band);
        harmonic_idx = valid_indices(band_peak_idx);

        % Retrieve the corresponding frequency
        %freq_harmonic = freqaxis(harmonic_idx);

        % Find the corresponding amplitude
        harmamplitudelist(i, 1) = 2*ss_amp;
       
    end

end


% Changing data unit; V --> Pa. 
% Go from Rx amplifier output (mV) to Rx output (mV) to emitted pressures by microbubbles (Pa)

function Pressureamplitude = getPressureamplitude(boosteramplifier_output_mV)
    
    %from booster amplifier output (mV) to Rx output (mV)
    slope_calibline_BA = 26.738933426280940; % mV output / mV input
    Rx_output_mV = boosteramplifier_output_mV/slope_calibline_BA;

    %from Rx output (mV) to cavitation pressure amplitude (Pa)
    slope_calibline_Rx = 0.0222e3; %Pa / mV
    Pressureamplitude = Rx_output_mV * slope_calibline_Rx;

end


function [cavitation_levels, cavitation_dose] = getCavitationDose(inputstruct, frequency_Hz, cycles, type)


    switch type
        
        % 1 = harmonics
        case 1
            harmonics = [2,3,4,5,6];
            numMeasurements = 30; 
    
             % Prepardatacontainer
            cavitation_levels = zeros(numMeasurements, 1); 
            bandwidth = 20000; % 20 Khz
    
            %loop over all measurments within the burst of the inputstruct
            for i= 1:numMeasurements
        
                % Compute fourier for each of the 30 measurements
                timedata = inputstruct.(sprintf('Time_%d', i));
                voltdata = inputstruct.(sprintf('Volt_%d', i));
                [freqaxis, absolutefft, ~] =  fourier(timedata, voltdata, frequency_Hz, cycles, 0);

                % Make it single sided, and select the sensitive bandwith
                % of transduces (>10% sensitivity)
                f_min = 1e6; % 1 MHz
                f_max = 10e6; % 10 MHz
                
                % Retain only the values within the specified frequency range
                bw_indices = (freqaxis >= f_min) & (freqaxis <= f_max);
                freqaxis_filtered = freqaxis(bw_indices);
                absolutefft_filtered = getPressureamplitude(absolutefft(bw_indices)*1000); %filter, and V --> mV --> Pa
                
                % Initialize a row vector to store squared amplitudes for
                % all harmonics
                squared_amplitudes = [];
                
                % Compute cavitation level for each of the 30 measurements
                % bandwidth of 20 KHz
                for j = 1:length(harmonics) 
                    harmonic = frequency_Hz * harmonics(j);
                    lower_bound = harmonic - bandwidth / 2;
                    higher_bound = harmonic + bandwidth / 2;
    
                    % Find the indices in freqaxis that fall within the lower and higher bounds
                    indices = find(freqaxis_filtered >= lower_bound & freqaxis_filtered <= higher_bound);
    
                    % Square the corresponding values of absolutefft and add to squared_amplitudes
                    squared_amplitudes = [squared_amplitudes, absolutefft_filtered(indices).^2]; 
    
                end 
    
                RMS_ampl  = sqrt(mean(squared_amplitudes));
                cavitation_levels(i,1) = RMS_ampl;
    
            end 
    
            % Compute cavitation dose for the burst
            cavitation_dose = sum(cavitation_levels);
    
        % 2 = ultraharmonics
        case 2
            %m/2*f0,
            ultraharmonics = [3, 5, 7, 9, 11];
            numMeasurements = 30; 
    
             % Prepardatacontainer
            cavitation_levels = zeros(numMeasurements, 1); 
            bandwidth = 20000; % 20 Khz
    
            %loop over all measurments within the burst of the inputstruct
            for i= 1:numMeasurements
        
                % Compute fourier for each of the 30 measurements
                timedata = inputstruct.(sprintf('Time_%d', i));
                voltdata = inputstruct.(sprintf('Volt_%d', i));
                [freqaxis, absolutefft, ~] =  fourier(timedata, voltdata, frequency_Hz, cycles, 0);
                
                % Make it single sided, and select the sensitive bandwith
                % of transduces (>10% sensitivity)
                f_min = 1e6; % 1 MHz
                f_max = 10e6; % 10 MHz
                
                % Retain only the values within the specified frequency range
                bw_indices = (freqaxis >= f_min) & (freqaxis <= f_max);
                freqaxis_filtered = freqaxis(bw_indices);
                absolutefft_filtered = getPressureamplitude(absolutefft(bw_indices)*1000); %filter, and V --> mV --> Pa

                % Initialize a row vector to store squared amplitudes for each harmonic
                squared_amplitudes = [];
                
                % Compute cavitation level for each of the 30 measurements
                % bandwidth of 20 KHz
                for j = 1:length(ultraharmonics) 
                    ultraharmonic = frequency_Hz * (ultraharmonics(j)/2);
                    lower_bound = ultraharmonic - bandwidth / 2;
                    higher_bound = ultraharmonic + bandwidth / 2;
    
                    % Find the indices in freqaxis that fall within the lower and higher bounds
                    indices = find(freqaxis_filtered >= lower_bound & freqaxis_filtered <= higher_bound);
    
                    % Square the corresponding values of absolutefft and add to squared_amplitudes
                    squared_amplitudes = [squared_amplitudes, absolutefft_filtered(indices).^2];
    
                end 
    
                RMS_ampl  = sqrt(mean(squared_amplitudes));
                cavitation_levels(i,1) = RMS_ampl;
    
            end 
    
            % Compute cavitation dose for the burst
            cavitation_dose = sum(cavitation_levels);
    
        % 3 = interial
        case 3

            %exclude 360 KHz around harmonics
            %exclude 100 kHz around ultraharmonics

            numMeasurements = 30; 
            ultraharmonics = [3, 5, 7, 9, 11];
            harmonics = [2,3,4,5,6];

    
            % Prepardatacontainer
            cavitation_levels = zeros(numMeasurements, 1); 
            harmonic_bandwidth = 360e3; 
            ultraharmonic_bandwidth = 100e3;
            lower_frequency_threshold = frequency_Hz + harmonic_bandwidth/2; % 1.68 MHz

    
            %loop over all measurments within the burst of the inputstruct
            for i= 1:numMeasurements
        
                % Compute fourier for each of the 30 measurements
                timedata = inputstruct.(sprintf('Time_%d', i));
                voltdata = inputstruct.(sprintf('Volt_%d', i));
                [freqaxis, absolutefft, ~] =  fourier(timedata, voltdata, frequency_Hz, cycles, 0);
                
                % Make it single sided, and select the sensitive bandwith
                % of transduces (>10% sensitivity)
                f_min = 1e6; % 1 MHz
                f_max = 10e6; % 10 MHz
                
                % Retain only the values within the specified frequency range
                bw_indices = (freqaxis >= f_min) & (freqaxis <= f_max);
                freqaxis_filtered = freqaxis(bw_indices);
                absolutefft_filtered = getPressureamplitude(absolutefft(bw_indices)*1000); %filter, and V --> mV --> Pa

                % Find indices of frequencies to exclude (harmonics and ultraharmonics)
                excluded_indices = false(size(freqaxis_filtered)); % Initialize mask for excluded indices

                % Exclude harmonics
                for j = 1:length(harmonics)
                    harmonic = frequency_Hz * harmonics(j);
                    lower_bound = harmonic - harmonic_bandwidth / 2;
                    higher_bound = harmonic + harmonic_bandwidth / 2;
                    excluded_indices = excluded_indices | (freqaxis_filtered >= lower_bound & freqaxis_filtered <= higher_bound);
                end

                % Exclude ultraharmonics
                for j = 1:length(ultraharmonics)
                    ultraharmonic = frequency_Hz * (ultraharmonics(j) / 2);
                    lower_bound = ultraharmonic - ultraharmonic_bandwidth / 2;
                    higher_bound = ultraharmonic + ultraharmonic_bandwidth / 2;
                    excluded_indices = excluded_indices | (freqaxis_filtered >= lower_bound & freqaxis_filtered <= higher_bound);
                end

                % Exclude frequencies below the threshold (and negative spectrum)
                excluded_indices = excluded_indices | (freqaxis_filtered < lower_frequency_threshold);
    
                % Find the remaining indices (not excluded)
                included_indices = ~excluded_indices;
    
                % Square the corresponding values of absolutefft and calculate RMS
                squared_amplitudes = absolutefft_filtered(included_indices).^2;
                RMS_ampl = sqrt(mean(squared_amplitudes));
                cavitation_levels(i, 1) = RMS_ampl;
            end
    
            % Compute cavitation dose for the burst
            cavitation_dose = sum(cavitation_levels);


end

end

%% Define data
MB_1 = PPO2_TherapeuticBurst_23_4V;
noMB_1 = PPO2_Background_1_23_4V;
noMB_2 = PPO2_Background_2_23_4V;
noMB_3 = PPO2_Background_3_23_4V;

%%
figure
fourier(MB_1.Time_1, getPressureamplitude(MB_1.Volt_1), 1.5e6, 100, 1 );
figure
fourier(noMB_1.Time_1, getPressureamplitude(noMB_1.Volt_1), 1.5e6, 100, 1 );


%% Explore data (Pa)
% Define shades of red and blue (from light to dark)
red = [1 0 0];
blues = [0.6 0.6 1; 0.3 0.3 1; 0 0 1];

pulse_no = (1:30);

% Loop over voltages, and make a figure with 1 subplot with 6 lines (for all MB and noMB measurements) per voltage
figure
[~, ~, amplitudelist_1] = meanfourierburst(MB_1, 1.5e6, 100, 0);                  % Get a list of the amplitudes of the fundamental of each pulse withing the therapeutic burst
[~, ~, amplitudelist_4] = meanfourierburst(noMB_1, 1.5e6, 100, 0);
[~, ~, amplitudelist_5] = meanfourierburst(noMB_2, 1.5e6, 100, 0);
[~, ~, amplitudelist_6] = meanfourierburst(noMB_3, 1.5e6, 100, 0);
hold on
% Plot the blue lines (noMB)
plot(pulse_no, getPressureamplitude(amplitudelist_4 * 1e3), 'Color', blues(1,:), 'LineWidth', 1);
plot(pulse_no, getPressureamplitude(amplitudelist_5 * 1e3), 'Color', blues(2,:), 'LineWidth', 1);
plot(pulse_no, getPressureamplitude(amplitudelist_6 * 1e3), 'Color', blues(3,:), 'LineWidth', 1);
% Plot the red linee (MB)
hold on
plot(pulse_no, getPressureamplitude(amplitudelist_1 * 1e3), 'Color', red(1,:), 'LineWidth', 1);
grid on

% Set axis limits and labels and add legend
ylim([0 1000]);
legend({'NoMB rep1', 'NoMB rep2', 'NoMB rep3', 'MB'}, ...
       'Location', 'northeast', 'FontSize', 6);
ylabel('Pressure amplitude (Pa)');
xlabel('Measurement');
title('Evolution of amplitude of the fundamental throughout burst');
hold off; 


%% Explore data (Pa)

pulse_no = (1:30);

% Loop over voltages, and make a figure with 1 subplot with 6 lines (for all MB and noMB measurements) per voltage
figure
harm_amplitudelist_1 = getHarmAmpl(MB_1, 1.5e6, 100);                         % Get a list of the amplitudes of the 1st harmonic of each pulse withing the therapeutic burst
harm_amplitudelist_4 = getHarmAmpl(noMB_1, 1.5e6, 100);
harm_amplitudelist_5 = getHarmAmpl(noMB_2, 1.5e6, 100);
harm_amplitudelist_6 = getHarmAmpl(noMB_3, 1.5e6, 100);
hold on
% Plot the blue lines (noMB)
plot(pulse_no, getPressureamplitude(harm_amplitudelist_4) * 1e3, 'Color', blues(1,:), 'LineWidth', 1);
plot(pulse_no, getPressureamplitude(harm_amplitudelist_5) * 1e3, 'Color', blues(2,:), 'LineWidth', 1);
plot(pulse_no, getPressureamplitude(harm_amplitudelist_6) * 1e3, 'Color', blues(3,:), 'LineWidth', 1);
hold on
% Plot the red line (MB)
plot(pulse_no, getPressureamplitude(harm_amplitudelist_1) * 1e3, 'Color', red(1,:), 'LineWidth', 1);
grid on


% Set axis limits and labels and add legend
ylim([0 100]);
legend({'noMB rep1', 'noMB rep2', 'noMB rep3', 'MB'}, ...
       'Location', 'northeast', 'FontSize', 6);
ylabel('Pressure amplitude (Pa)');
xlabel('Measurement');
title('Evolution of first harmonic amplitude  throughout burst');



hold off

%% Compute mean fourier domain of the burst
% mean FD of  the therapeutic burst on MB
% mean FD of the 3 therapeutic burst on noMB

% MB
[freqaxis_MB, meanfftlist_MB(:, 1), meanamplitudelist_MB(:, 1)] = meanfourierburst(MB_1, 1.5e6, 100, 0);

% noMB
[freqaxis_noMB, meanfftlist_noMB(:, 1), meanamplitudelist_noMB(:, 1)] = meanfourierburst(noMB_1, 1.5e6, 100, 0);
[~, meanfftlist_noMB(:, 2), meanamplitudelist_noMB(:, 2)] = meanfourierburst(noMB_2, 1.5e6, 100, 0);
[~, meanfftlist_noMB(:, 3), meanamplitudelist_noMB(:, 3)] = meanfourierburst(noMB_3, 1.5e6, 100, 0);

% Compute the mean absolute fft over the three repeats
meanfft_MB = mean(meanfftlist_MB, 2); %mean of 1 repeat is the repeat itself
meanfft_noMB = mean(meanfftlist_noMB, 2);

% Plot
figure
subplot(1,2,1)
plot(freqaxis_MB/1e6, getPressureamplitude(meanfft_MB*1e3), 'r', freqaxis_noMB/1e6, getPressureamplitude(meanfft_noMB*1e3), 'b') %Pa and MHz
legend('MB', 'noMB')
title('Mean Fourier Spectrum of therapeutic burst at 0.4 MPa');
xlim([-15 15]); xlabel('Frequency (MHz)'), ylabel(' Pressure amplitude (Pa)')
grid on
xticks(-15:1.5:15)   % Change x-axis grid spacing (adjust as needed)

        
% From here on, data is displayed and saved as pressure amplitudes
% Going to dB for a double sided spectrum
meanfft_MB_dB_ds = 20*log10(getPressureamplitude(meanfft_MB) / 1);
meanfft_noMB_dB_ds = 20*log10(getPressureamplitude(meanfft_noMB) / 1);

% Going to dB 
subplot(1,2,2)
plot(freqaxis_MB/1e6, meanfft_MB_dB_ds, 'r', freqaxis_noMB/1e6, meanfft_noMB_dB_ds, 'b') %mV and MHz
legend('MB', 'noMB')
title('Mean Decibels Fourier Spectrum of therapeutic burst at 0.4 MPa');
xlim([-15 15]); ylim([-75 -10]); xlabel('MHz'), ylabel('dB')
grid on
xticks(-15:1.5:15)   % Change x-axis grid spacing (adjust as needed)

    
  
% plot signle sided fourier spectrum
meanfft_MB_dB_ss = 20*log10(getPressureamplitude(meanfft_MB*2) / 1); % just double the ampltidues
meanfft_noMB_dB_ss = 20*log10(getPressureamplitude(meanfft_noMB*2) / 1);

figure
plot(freqaxis_MB/1e6, meanfft_MB_dB_ss, 'r', freqaxis_noMB/1e6, meanfft_noMB_dB_ss, 'b') %mV and MHz
legend('MB', 'noMB')
title('Mean Decibels Fourier Spectrum of therapeutic burst at 0.4 MPa');
xlim([0 15]); ylim([-75 -10]); xlabel('Frequency (MHz)'), ylabel('Pressure amplitude (dB)') %only plot positive halve
grid on
xticks(0:1.5:15)   % Change x-axis grid spacing (adjust as needed)

    

%% Compute cavitaiton dose: Harmonic
[levels_MB_1_h, dose_MB_1_h] = getCavitationDose(MB_1, 1.5e6, 100, 1);
[levels_noMB_1_h, dose_noMB_1_h] = getCavitationDose(noMB_1, 1.5e6, 100, 1);
[levels_noMB_2_h, dose_noMB_2_h] = getCavitationDose(noMB_2, 1.5e6, 100, 1);
[levels_noMB_3_h, dose_noMB_3_h] = getCavitationDose(noMB_3, 1.5e6, 100, 1);

% now the mean
meandose_MB_h = dose_MB_1_h;
meandose_noMB_h = (dose_noMB_1_h + dose_noMB_2_h + dose_noMB_3_h)/3;
meandifference_h = meandose_MB_h-meandose_noMB_h; % ACTUAL CAVITATION DOSE

% Figure
measurements = (1:30);
figure
plot(measurements, levels_MB_1_h, 'r*-')
hold on
plot(measurements, levels_noMB_1_h, 'g*-')
hold on
plot(measurements, levels_noMB_2_h, 'b*-')
hold on
plot(measurements, levels_noMB_3_h, 'y*-')
grid on
title('Harmonic cavitation levels throughout burst')
legend('MB', 'noMB_1', 'noMB_2', 'noMB_3')
hold off

%% Compute cavitaiton doses: Ultraharmonic
[levels_MB_1_uh, dose_MB_1_uh] = getCavitationDose(MB_1, 1.5e6, 100, 2);
[levels_noMB_1_uh, dose_noMB_1_uh] = getCavitationDose(noMB_1, 1.5e6, 100, 2);
[levels_noMB_2_uh, dose_noMB_2_uh] = getCavitationDose(noMB_2, 1.5e6, 100, 2);
[levels_noMB_3_uh, dose_noMB_3_uh] = getCavitationDose(noMB_3, 1.5e6, 100, 2);

% now the mean
meandose_MB_uh = dose_MB_1_uh;
meandose_noMB_uh = (dose_noMB_1_uh + dose_noMB_2_uh + dose_noMB_3_uh)/3;
meandifference_uh = meandose_MB_uh-meandose_noMB_uh; % ACTUAL CAVITATION DOSE

% Figure
measurements = (1:30);
figure
plot(measurements, levels_MB_1_uh, 'r*-')
hold on
plot(measurements, levels_noMB_1_uh, 'g*-')
hold on
plot(measurements, levels_noMB_2_uh, 'b*-')
hold on
plot(measurements, levels_noMB_3_uh, 'y*-')
grid on
title('Ultraharmonic cavitation levels throughout burst')
legend('MB', 'noMB_1', 'noMB_2', 'noMB_3')
hold off


%% Compute cavitaiton doses: Inertial
[levels_MB_1_in, dose_MB_1_in] = getCavitationDose(MB_1, 1.5e6, 100, 3);
[levels_noMB_1_in, dose_noMB_1_in] = getCavitationDose(noMB_1, 1.5e6, 100, 3);
[levels_noMB_2_in, dose_noMB_2_in] = getCavitationDose(noMB_2, 1.5e6, 100, 3);
[levels_noMB_3_in, dose_noMB_3_in] = getCavitationDose(noMB_3, 1.5e6, 100, 3);

% now the mean
meandose_MB_in = dose_MB_1_in;
meandose_noMB_in = (dose_noMB_1_in + dose_noMB_2_in + dose_noMB_3_in)/3;
meandifference_in = meandose_MB_in-meandose_noMB_in; % ACTUAL CAVITATION DOSE

% Figure
measurements = (1:30);
figure
plot(measurements, levels_MB_1_in, 'r*-')
hold on
plot(measurements, levels_noMB_1_in, 'g*-')
hold on
plot(measurements, levels_noMB_2_in, 'b*-')
hold on
plot(measurements, levels_noMB_3_in, 'y*-')
hold on
title('Inertial cavitation levels throughout burst')
legend('MB', 'noMB_1', 'noMB_2', 'noMB_3')
grid on
hold off


