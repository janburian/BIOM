clc
close all
clear all

%% Loading data
xml_file = './data/ekg.XML'; 
ekg_data = xml2struct(xml_file);

sequences = ekg_data.Children(16).Children(2).Children(16).Children(2).Children;

% Data extraction
MDC_ECG_LEAD_I = sequences(4).Children(2).Children(4).Children(6).Children.Data; 
MDC_ECG_LEAD_II = sequences(6).Children(2).Children(4).Children(6).Children.Data; 
MDC_ECG_LEAD_III = sequences(8).Children(2).Children(4).Children(6).Children.Data; 
MDC_ECG_LEAD_AVR = sequences(10).Children(2).Children(4).Children(6).Children.Data;
MDC_ECG_LEAD_AVL = sequences(12).Children(2).Children(4).Children(6).Children.Data;
MDC_ECG_LEAD_AVF = sequences(14).Children(2).Children(4).Children(6).Children.Data;
MDC_ECG_LEAD_V1 = sequences(16).Children(2).Children(4).Children(6).Children.Data;
MDC_ECG_LEAD_V2 = sequences(18).Children(2).Children(4).Children(6).Children.Data;
MDC_ECG_LEAD_V3 = sequences(20).Children(2).Children(4).Children(6).Children.Data;
MDC_ECG_LEAD_V4 = sequences(22).Children(2).Children(4).Children(6).Children.Data;
MDC_ECG_LEAD_V5 = sequences(24).Children(2).Children(4).Children(6).Children.Data;
MDC_ECG_LEAD_V6 = sequences(26).Children(2).Children(4).Children(6).Children.Data;

% Converting to numeric vectors
MDC_ECG_LEAD_I_num_vec = str2double(strsplit(MDC_ECG_LEAD_I(1:end-1), ' ')); 
MDC_ECG_LEAD_II_num_vec = str2double(strsplit(MDC_ECG_LEAD_II(1:end-1), ' '));
MDC_ECG_LEAD_III_num_vec = str2double(strsplit(MDC_ECG_LEAD_III(1:end-1), ' '));
MDC_ECG_LEAD_AVR_num_vec = str2double(strsplit(MDC_ECG_LEAD_AVR(1:end-1), ' '));
MDC_ECG_LEAD_AVL_num_vec = str2double(strsplit(MDC_ECG_LEAD_AVL(1:end-1), ' '));
MDC_ECG_LEAD_AVF_num_vec = str2double(strsplit(MDC_ECG_LEAD_AVF(1:end-1), ' '));

%% Data visualisation
% Generate sample data for three vectors
x = linspace(0, 10, 10000);  % Assuming 10,000 steps

% Create subplots
figure;

subplot(6, 1, 1);
plot(x, MDC_ECG_LEAD_I_num_vec);
set(gca,'ytick',[])
title('I');
xlabel('t [s]')
ylabel('U [\muV]')

subplot(6, 1, 2);
plot(x, MDC_ECG_LEAD_II_num_vec);
set(gca,'ytick',[])
title('II');
xlabel('t [s]')
ylabel('U [\muV]')

subplot(6, 1, 3);
plot(x, MDC_ECG_LEAD_III_num_vec);
set(gca,'ytick',[])
title('III');
xlabel('t [s]')
ylabel('U [\muV]')

subplot(6, 1, 4);
plot(x, MDC_ECG_LEAD_AVR_num_vec);
set(gca,'ytick',[])
title('aVR');
xlabel('t [s]')
ylabel('U [\muV]')

subplot(6, 1, 5);
plot(x, MDC_ECG_LEAD_AVL_num_vec);
set(gca,'ytick',[])
title('aVL');
xlabel('t [s]')
ylabel('U [\muV]')

subplot(6, 1, 6);
plot(x, MDC_ECG_LEAD_AVF_num_vec);
set(gca,'ytick',[])
title('aVF');
xlabel('t [s]')
ylabel('U [\muV]')

%sgtitle('EKG') 

%% Filter design
fs = 1000; % sampling frequency
f_cutoff = 50; % cutoff frequency
n = 2; % filter order
Wn = f_cutoff/(fs/2); % normalized cutoff frequency

% Lowpass filter
[b, a] = butter(2, Wn, 'low');

freqz(b, a, fs);
title('Frequency response')

% Filter application
filtered_signal = filtfilt(b, a, MDC_ECG_LEAD_II_num_vec);

% Visualisation of the original and filtered signal
figure;
subplot(2,1,1);
plot(MDC_ECG_LEAD_II_num_vec); 
title('Original ECG signal');
subplot(2,1,2); 
plot(filtered_signal); 
title('Filtered ECG signal');

%% 2 periods
figure;
subplot(2,1,1);
plot(x, MDC_ECG_LEAD_II_num_vec); 
xlabel('t [s]')
ylabel('U [\muV]')
xlim([3.8 5.5]-3.8)
title('Original ECG signal');
subplot(2,1,2); 
plot(x, filtered_signal); 
xlabel('t [s]')
ylabel('U [\muV]')
xlim([3.8 5.5]-3.8)
title('Filtered ECG signal');
sgtitle('Comparison of the original and filtered signals')

%%
figure;
subplot(2,1,1);
spectrogram(MDC_ECG_LEAD_II_num_vec,'yaxis', fs) 
ylabel('Normalized frequency')
title('Original ECG signal');
subplot(2,1,2); 
spectrogram(filtered_signal,'yaxis', fs)
ylabel('Normalized frequency')
title('Filtered ECG signal');
sgtitle('Comparison of the spectrograms')

%% 1 period
signal_1_period = filtered_signal(1, 3980:4700); 
x = linspace(0, length(signal_1_period) * 0.001, length(signal_1_period));

figure;
hold on
plot(x, signal_1_period); 
xlabel('t [s]')
ylabel('U [\muV]')
title('Filtered ECG signal (1 period)');

%% Segments
% P_vec
% PQ_vec
% QRS_vec
% ST_vec
% T_vec

% Plotting the data
figure
hold on
plot(x(1:165), signal_1_period(1:165), 'Color', '#0072BD', 'DisplayName', 'P')
plot(x(165:215), signal_1_period(165:215), 'Color', '#D95319', 'DisplayName', 'PQ')
plot(x(215:230), signal_1_period(215:230), 'Color', '#EDB120', 'DisplayName', 'Q')
plot(x(230:248), signal_1_period(230:248), 'Color', '#7E2F8E', 'DisplayName', 'QR')
scatter(x(248), signal_1_period(248), 'filled', 'DisplayName', 'R')
plot(x(248:282), signal_1_period(248:282), 'Color', '#77AC30', 'DisplayName', 'S')
plot(x(282:330), signal_1_period(282:330), 'Color', '#4DBEEE', 'DisplayName', 'ST')
plot(x(330:end), signal_1_period(330:end), 'Color', '#A2142F', 'DisplayName', 'T')
xlabel('t [s]')
ylabel('U [\muV]')

% Displaying the legend
legend('show')
title('ECG signal (PQRST)')

%% FFT
N = length(signal_1_period);  % Length of the signal
Fs = 1000;  % Sampling frequency
frequencies = Fs * (0:(N/2))/N;  % Frequency axis for positive frequencies

% Perform FFT
fft_result = fft(signal_1_period);

% Take only the positive half (single-sided spectrum)
single_sided_spectrum = 2/N * abs(fft_result(1:N/2 + 1));

% Compute the power spectrum (square of magnitude)
power_spectrum = single_sided_spectrum.^2;

figure;
subplot(2, 1, 1);
plot(x, signal_1_period);
title('ECG signal (1 period)');
grid on;
xlabel('t [s]');
ylabel('U [\muV]');
subplot(2, 1, 2);
plot(frequencies, power_spectrum)
xlim([0 30])
xlabel('Frequency [Hz]')
ylabel('Power')
title('Power Spectrum of the Signal')
sgtitle('ECG signal and power spectrum');

%% Find R-peaks in the ECG signal
[peaks, locs] = findpeaks(filtered_signal, 'MinPeakHeight', 0.6, 'MinPeakDistance', 300);
t = linspace(0, 10, 10000);  % Time vector

% Visualize the ECG signal and the detected R-peaks
figure;
plot(t, filtered_signal);
hold on;
plot(t(locs), peaks, 'rx', 'MarkerSize', 10);
title('ECG Signal with detected R-peaks');
xlabel('t [s]');
ylabel('U [\muV]');
legend('ECG Signal', 'R-peaks');

%% Counting BPM
% Calculate Inter-Beat Intervals (IBIs) in seconds
ibis = diff(t(locs));

% Convert IBIs to Beats Per Minute (BPM)
bpm = 60 ./ mean(ibis);

% Display the BPM
disp(['Heart rate (BPM): ', num2str(bpm)]);




