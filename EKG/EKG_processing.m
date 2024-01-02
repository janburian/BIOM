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
title('I');
xlabel('t [s]')
ylabel('U [\muV]')

subplot(6, 1, 2);
plot(x, MDC_ECG_LEAD_II_num_vec);
title('II');
xlabel('t [s]')
ylabel('U [\muV]')

subplot(6, 1, 3);
plot(x, MDC_ECG_LEAD_III_num_vec);
title('III');
xlabel('t [s]')
ylabel('U [\muV]')

subplot(6, 1, 4);
plot(x, MDC_ECG_LEAD_AVR_num_vec);
title('aVR');
xlabel('t [s]')
ylabel('U [\muV]')

subplot(6, 1, 5);
plot(x, MDC_ECG_LEAD_AVL_num_vec);
title('aVL');
xlabel('t [s]')
ylabel('U [\muV]')

subplot(6, 1, 6);
plot(x, MDC_ECG_LEAD_AVF_num_vec);
title('aVF');
xlabel('t [s]')
ylabel('U [\muV]')

sgtitle('EKG') 

%% Filter design
fs = 1000;
f_cutoff = 50;

% Lowpass filter
[b, a] = butter(2, f_cutoff/(fs/2), 'low');

freqz(b, a, 1024, fs);

%% 
% Aplikace filtru
filtered_signal = filtfilt(b, a, MDC_ECG_LEAD_II_num_vec);

% Zobrazeni puvodniho a vyfiltrovaneho signalu
figure;
subplot(2,1,1);
plot(MDC_ECG_LEAD_II_num_vec); 
title('Original EKG signal');
subplot(2,1,2); 
plot(filtered_signal); 
title('Filtered EKG signal');

%% P, QRS, T 
% Urceni adaptivnich prahu
mean_threshold = mean(filtered_signal);
std_threshold = std(filtered_signal);

% Aplikace prahu na detekci vrcholu P, QRS, T
threshold_P = mean_threshold * 0.6;
threshold_QRS = mean_threshold * 1.2;
threshold_T = mean_threshold * 0.8;

% Detekce vrcholu P, QRS, T
[pks_P, locs_P] = findpeaks(filtered_signal, 'MinPeakHeight', threshold_P);
[pks_QRS, locs_QRS] = findpeaks(-filtered_signal, 'MinPeakHeight', threshold_QRS);
[pks_T, locs_T] = findpeaks(filtered_signal, 'MinPeakHeight', threshold_T);

% Zobrazeni v�sledku
figure;
plot(filtered_signal);
hold on;
plot(locs_P, pks_P, 'rv', 'MarkerFaceColor', 'r');
plot(locs_QRS, -pks_QRS, 'bx', 'MarkerFaceColor', 'b');
plot(locs_T, pks_T, 'g^', 'MarkerFaceColor', 'g');
title('Detekce vrchol� P, QRS, T');
xlabel('�as');
ylabel('Amplituda');
legend('EKG sign�l', 'Vrcholy P', 'Vrcholy QRS', 'Vrcholy T');
hold off;

% % Vypocet delek intervalu
% PQ_interval = locs_QRS(1) - locs_P;
% QRS_interval = locs_QRS(2) - locs_QRS(1);
% QT_interval = locs_T(1) - locs_QRS(2);

