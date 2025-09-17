% =========================================================================
% Demonstration of key limitations of short-time Fourier transform (STFT):
% 1. Non-concentrated time-frequency energy distribution
% 2. Amplitude loss
% This demonstration was written by Yabo Wang at Xidian University
% =========================================================================
clc;clear;close all
%% Basic parameters
fs=1024; % Sampling frequency (Hz)
N=1; % Signal duration set to 1s
time=0:1/fs:N-1/fs; % Time vector
Nt=length(time); % Number of sampling points
fre=(0:(Nt-1)/2)/Nt*2*fs/2; % Frequency vector

%% Signal model: multi-component signal
a=1; % Constant amplitude of 1 for all signal components
s1= a*cos(2*pi*(35*time)); % Component 1: Fixed frequency (35 Hz) signal
s2= a*cos(2*pi*( 80*time+140*time.^2)); % Component 2: Linear frequency modulated (LFM) signal
s3= a*cos(2*pi*( 195*time+50*time.^2+70*time.^3)); % Component 3: Quadratic frequency modulated (QFM) signal
s=s1'+s2'+s3';

%% Ideal instantaneous frequencies (IF) for each signal component
IF1(1,1:Nt)=35; % IF of signal component 1
IF2=80+280*time; % IF of signal component 2
IF3=195+100*time+210*time.^2; % IF of signal component 3

%% Signal waveform and instantaneous frequencies
figure;
subplot(2,1,1)
plot(time,s,'linewidth',1);
xlabel('Time / s');ylabel('Amplitude');
xlim([0 max(time)]);
title('(a) Signal wave');
set(gca, 'FontSize', 10, 'Fontname', 'Times New Roman');

subplot(2,1,2)
plot(time,IF1,'linewidth',1.5);
hold on;
plot(time,IF2,'linewidth',1.5);
plot(time,IF3,'linewidth',1.5);
hold off;
xlabel('Time / s');ylabel('Frequency / Hz');
ylim([0 fs/2]);
title('(b) Ideal instantaneous frequency trajectories');
legend('Component 1 (35 Hz)','Component 2 (LFM)','Component 3 (QFM)',"Location","northwest");
set(gca, 'FontSize', 10, 'Fontname', 'Times New Roman');

%% STFT
w_length=170; % Window length
[tfr] = STFT_Y(s,w_length); % STFT

% The time-frequency representation (TFR) result of STFT
figure;
imagesc(time,fre,abs(tfr));
xlabel('Time / s');ylabel('Frequency / Hz');
title('The TFR result of STFT');
axis xy;
caxis([0 1]);
set(gca, 'FontSize', 10, 'Fontname', 'Times New Roman');

% 3D visualization of STFT
figure
mesh(time,fre,abs(tfr),'LineWidth', 1);
xlabel('Time / s');ylabel('Frequency / Hz');
zlabel('Amplitude');
xlim([0 max(time)]);
ylim([0 max(fre)]);
zlim([0 1]);
axis xy
shading interp;
set(gcf, 'Color', 'white');
view(-19,72)
caxis([0 1]);
title('The 3D view of TFR');
set(gca, 'FontSize', 10, 'Fontname', 'Times New Roman');

%% Frequency slice analysis (TFR slices at specific time points)
% Key amplitude characteristics note
% All signal components exhibit energy diffusion
% Component 1 (35 Hz): Almost no amplitude loss, only energy diffusion
% Component 2 (LFM): Constant amplitude loss
% Component 3 (QFM): Amplitude loss that increases over time
figure
% Slice 1: First time point
f_slice1=100;
t1=time(f_slice1); % 0.0967s
subplot(2,2,1)
plot(fre,abs(tfr(:,f_slice1)),"LineWidth",1.5);
hold on
xlim([0 max(fre)]);
ylim([0 1]);
xlabel('Frequency / Hz');ylabel('Amplitude');
title(['(a) The frequency slice at ', num2str(t1, '%.4f'), 's']);
set(gca, 'FontSize', 10, 'Fontname', 'Times New Roman');

% Slice 2: Second time point
f_slice2=350;
t2=time(f_slice2); % 0.3408s
subplot(2,2,2)
plot(fre,abs(tfr(:,f_slice2)),"LineWidth",1.5);
xlim([0 max(fre)]);
ylim([0 1]);
xlabel('Frequency / Hz');ylabel('Amplitude');
title(['(b) The frequency slice at ', num2str(t2, '%.4f'), 's']);
set(gca, 'FontSize', 10, 'Fontname', 'Times New Roman');

% Slice 3: Third time point
f_slice3=700;
t3=time(f_slice3); % 0.6826s
subplot(2,2,3)
plot(fre,abs(tfr(:,f_slice3)),"LineWidth",1.5);
xlim([0 max(fre)]);
ylim([0 1]);
xlabel('Frequency / Hz');ylabel('Amplitude');
title(['(c) The frequency slice at ', num2str(t3, '%.4f'), 's']);
set(gca, 'FontSize', 10, 'Fontname', 'Times New Roman');

% Slice 4: Fourth time point
f_slice4=950;
t4=time(f_slice4); % 0.9268s
subplot(2,2,4)
plot(fre,abs(tfr(:,f_slice4)),"LineWidth",1.5);
hold off
xlim([0 max(fre)]);
ylim([0 1]);
xlabel('Frequency / Hz');ylabel('Amplitude');
title(['(d) The frequency slice at ', num2str(t4, '%.4f'), 's']);
set(gca, 'FontSize', 10, 'Fontname', 'Times New Roman');
set(gca, 'FontSize', 10, 'Fontname', 'Times New Roman');