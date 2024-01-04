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
f_cutoff = 25; 

% Lowpass filter
[b, a] = butter(2, f_cutoff/(fs/2), 'low');

freqz(b, a, 1024, fs);

% Filter application
filtered_signal = filtfilt(b, a, MDC_ECG_LEAD_II_num_vec);

% Visualisation of the original and filtered signal
figure;
subplot(2,1,1);
plot(MDC_ECG_LEAD_II_num_vec); 
title('Original EKG signal');
subplot(2,1,2); 
plot(filtered_signal); 
title('Filtered EKG signal');

%% 1 period
figure
plot(x, filtered_signal);
xlabel('t [s]')
ylabel('U [\muV]')
yline(-490)
xlim([3.8 4.8]-3.8)
title('1 period');

%%
% Urèení adaptivních prahù
mean_threshold = mean(filtered_signal);
std_threshold = std(filtered_signal);

% Aplikace prahù na detekci vrcholù P, QRS, T
threshold_P = mean_threshold * 30;
threshold_Q = mean_threshold * 1.2;
threshold_R = mean_threshold * 5;
threshold_S = mean_threshold * 0.04;
threshold_T = mean_threshold * 0.8;

min_distance_P = 700;
min_distance_Q = 0.5;
min_distance_R = 0.5;
min_distance_S = 0.5;
min_distance_T = 0.5;

% Detekce vrcholù P, Q, R, S, T
[pks_P, locs_P] = findpeaks(filtered_signal, 'MinPeakHeight', threshold_P, 'MinPeakDistance', min_distance_P);
[pks_Q, locs_Q] = findpeaks(-filtered_signal, 'MinPeakHeight', threshold_Q, 'MinPeakDistance', min_distance_Q);
[pks_R, locs_R] = findpeaks(filtered_signal, 'MinPeakHeight', threshold_R, 'MinPeakDistance', min_distance_R);
[pks_S, locs_S] = findpeaks(-filtered_signal, 'MinPeakHeight', threshold_S, 'MinPeakDistance', min_distance_S);
[pks_T, locs_T] = findpeaks(filtered_signal, 'MinPeakHeight', threshold_T, 'MinPeakDistance', min_distance_T);

% Zobrazení výsledkù
figure;
plot(filtered_signal);
hold on;
plot(locs_P, pks_P, 'rv', 'MarkerFaceColor', 'r', 'DisplayName', 'Vrcholy P');
plot(locs_Q, -pks_Q, 'bv', 'MarkerFaceColor', 'b', 'DisplayName', 'Vrcholy Q');
plot(locs_R, pks_R, 'gx', 'MarkerFaceColor', 'g', 'DisplayName', 'Vrcholy R');
plot(locs_S, -pks_S, 'ms', 'MarkerFaceColor', 'm', 'DisplayName', 'Vrcholy S');
plot(locs_T, pks_T, 'co', 'MarkerFaceColor', 'c', 'DisplayName', 'Vrcholy T');
title('Detekce vrcholù P, Q, R, S, T');
xlabel('Èas');
ylabel('Amplituda');
legend('EKG signál', 'Location', 'Best');
hold off;

