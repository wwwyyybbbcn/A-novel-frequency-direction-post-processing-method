% =========================================================================
% Demonstration of time-frequency analysis for the simulated multi-component signal:
% Time-frequency analysis methods (STFT, SST, SET, MSST, SRT, and the proposed method)
% In terms of TFR, energy concentration, ridge extraction, and signal reconstruction accuracy.
% written by Yabo Wang, School of Electronic Engineering, Xidian University
% Repository: https://github.com/wwwyyybbbcn/A-novel-frequency-direction-post-processing-method.git
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
s_noise=awgn(s,5,"measured");

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

%% TFA and Renyi entropy
% TFA
w_length=170;
tic;[tfr] = STFT_Y(s,w_length);toc;
tic;[Ts] = SST_Y(s,w_length);  toc;
tic;[~,Te] = SET_Y(s,w_length); toc;
tic;[MTs] = MSST_Y_new(s,w_length,20);toc;
tic;[srt] = SRT(s,w_length,fs);toc;
tic;[tfr_pm] = PM_W(s,w_length);toc;

% Renyi entropy
renyi_stft=renyi(abs(tfr));
renyi_sst=renyi(abs(Ts));
renyi_set=renyi(abs(Te));
renyi_MSST=renyi(abs(MTs));
renyi_srt=renyi(abs(srt));
renyi_pm=renyi(abs(tfr_pm));

%% plot time-frequency representation (TFR)
% Fig. 1
figureUnits = 'inches';
figureWidth = 3.5;
figureHeight = 1.35;
fontSize = 7.5;
figure('Unit', figureUnits, 'Position', [1, 1, figureWidth, figureHeight]);
leftMargin = 0.45;
rightMargin = 0.04;
bottomMargin = 0.23; 
topMargin = 0.13;
hSpace = 0.55;
plotWidth = (figureWidth - leftMargin - rightMargin - hSpace)/2;
plotHeight = figureHeight - topMargin - bottomMargin;

subplot(1,2,1);
mesh(time,fre,abs(tfr),'LineWidth', 0.5);
xlabel('Time (s)', 'FontSize', fontSize);
ylabel('Frequency (Hz)', 'FontSize', fontSize,"Rotation",-70);
zlabel('Amplitude', 'FontSize', fontSize);
title('(a)', 'FontSize', fontSize);
xlim([0 max(time)]);
ylim([0 max(fre)]);
zlim([0 1]);
axis xy
set(gca, 'FontSize', 10, 'Fontname', 'Times New Roman');
view(-19,72)
shading interp;
set(gcf, 'Color', 'white');
xticks([0,0.4,0.8]);
yticks([0,200,400]);
caxis([0 1]);
pos1 = [leftMargin, bottomMargin, plotWidth, plotHeight];
set(gca,'Unit', figureUnits,'FontSize', fontSize, 'Fontname', 'Times New Roman', 'Position', pos1);

subplot(1,2,2);
mesh(time,fre,abs(tfr_pm),'LineWidth', 0.5);
xlabel('Time (s)', 'FontSize', fontSize);
ylabel('Frequency (Hz)', 'FontSize', fontSize,"Rotation",-70);
zlabel('Amplitude', 'FontSize', fontSize);
title('(b)', 'FontSize', fontSize);
xlim([0 max(time)]);
ylim([0 max(fre)]);
zlim([0 1]);
axis xy
set(gca, 'FontSize', 10, 'Fontname', 'Times New Roman');
view(-19,72)
shading interp;
set(gcf, 'Color', 'white');
xticks([0,0.4,0.8]);
yticks([0,200,400]);
caxis([0 1]);
pos2 = [leftMargin + plotWidth + hSpace, bottomMargin, plotWidth, plotHeight];
set(gca,'Unit', figureUnits,'FontSize', fontSize, 'Fontname', 'Times New Roman', 'Position', pos2);

% The TFR result of STFT
figure;
imagesc(time,fre,abs(tfr));
xlabel('Time (s)');
ylabel('Frequency (Hz)');
set(gca, 'FontSize', 10, 'Fontname', 'Times New Roman');
axis xy;
caxis([0 1]);
% 3D visualization of STFT
figure
mesh(time,fre,abs(tfr),'LineWidth', 1);
xlabel('Time (s)');
ylabel('Frequency (Hz)');
zlabel('Amplitude');
xlim([0 max(time)]);
ylim([0 max(fre)]);
zlim([0 1]);
axis xy
shading interp;
set(gcf, 'Color', 'white');
view(-19,72)
caxis([0 1]);
set(gca, 'FontSize', 10, 'Fontname', 'Times New Roman');

% The TFR result of SST
figure;
imagesc(time,fre,abs(Ts));
xlabel('Time (s)');
ylabel('Frequency (Hz)');
set(gca, 'FontSize', 10, 'Fontname', 'Times New Roman');
axis xy;
caxis([0 1]);
% 3D visualization of SST
figure
mesh(time,fre,abs(Ts),'LineWidth', 1);
xlabel('Time (s)');
ylabel('Frequency (Hz)');
zlabel('Amplitude');
xlim([0 max(time)]);
ylim([0 max(fre)]);
zlim([0 1]);
axis xy
shading interp;
set(gcf, 'Color', 'white');
view(-19,72)
caxis([0 1]);
set(gca, 'FontSize', 10, 'Fontname', 'Times New Roman');

% The TFR result of SET
figure;
imagesc(time,fre,abs(Te));
xlabel('Time (s)');
ylabel('Frequency (Hz)');
set(gca, 'FontSize', 10, 'Fontname', 'Times New Roman');
axis xy;
caxis([0 1]);
% 3D visualization of SET
figure
mesh(time,fre,abs(Te),'LineWidth', 1);
xlabel('Time (s)');
ylabel('Frequency (Hz)');
zlabel('Amplitude');
xlim([0 max(time)]);
ylim([0 max(fre)]);
zlim([0 1]);
axis xy
shading interp;
set(gcf, 'Color', 'white');
view(-19,72)
caxis([0 1]);
set(gca, 'FontSize', 10, 'Fontname', 'Times New Roman');

% The TFR result of MSST
figure;
imagesc(time,fre,abs(MTs));
xlabel('Time (s)');
ylabel('Frequency (Hz)');
set(gca, 'FontSize', 10, 'Fontname', 'Times New Roman');
axis xy;
caxis([0 1]);
% 3D visualization of MSST
figure
mesh(time,fre,abs(MTs),'LineWidth', 1);
xlabel('Time (s)');
ylabel('Frequency (Hz)');
zlabel('Amplitude');
xlim([0 max(time)]);
ylim([0 max(fre)]);
zlim([0 1]);
axis xy
shading interp;
set(gcf, 'Color', 'white');
view(-19,72)
caxis([0 1]);
set(gca, 'FontSize', 10, 'Fontname', 'Times New Roman');

% The TFR result of SRT
figure;
imagesc(time,fre,abs(srt));
xlabel('Time (s)');
ylabel('Frequency (Hz)');
set(gca, 'FontSize', 10, 'Fontname', 'Times New Roman');
axis xy;
caxis([0 1]);
% 3D visualization of SRT
figure
mesh(time,fre,abs(srt),'LineWidth', 1);
xlabel('Time (s)');
ylabel('Frequency (Hz)');
zlabel('Amplitude');
xlim([0 max(time)]);
ylim([0 max(fre)]);
zlim([0 1]);
axis xy
shading interp;
set(gcf, 'Color', 'white');
view(-19,72)
caxis([0 1]);
set(gca, 'FontSize', 10, 'Fontname', 'Times New Roman');

% The TFR result of the proposed method
figure;
imagesc(time,fre,abs(tfr_pm));
xlabel('Time (s)');
ylabel('Frequency (Hz)');
set(gca, 'FontSize', 10, 'Fontname', 'Times New Roman');
axis xy;
caxis([0 1]);
% 3D visualization of the proposed method
figure
mesh(time,fre,abs(tfr_pm),'LineWidth', 1);
xlabel('Time (s)');
ylabel('Frequency (Hz)');
zlabel('Amplitude');
xlim([0 max(time)]);
ylim([0 max(fre)]);
zlim([0 1]);
axis xy
shading interp;
set(gcf, 'Color', 'white');
view(-19,72)
caxis([0 1]);
set(gca, 'FontSize', 10, 'Fontname', 'Times New Roman');

%% Multi-ridge extraction
% The ridge extraction result of STFT
[Nf,~]=size(tfr);
Cs_STFT = brevridge_mult(abs(tfr), 1:Nf, 3, 0.02, 20);
Cs_STFT1 = Cs_STFT/Nf*fs/2; % Normalized ridges

% The ridge extraction result of SST
[Nf,~]=size(Ts);
Cs_SST = brevridge_mult(abs(Ts), 1:Nf, 3, 0.02, 20);
Cs_SST1 = Cs_SST/Nf*fs/2; % Normalized ridges

% The ridge extraction result of SET
[Nf,~]=size(Te);
Cs_SET = brevridge_mult(abs(Te), 1:Nf, 3, 0.02, 20);
Cs_SET1 = Cs_SET/Nf*fs/2; % Normalized ridges

% The ridge extraction result of MSST
[Nf,~]=size(MTs);
Cs_MSST = brevridge_mult(abs(MTs), 1:Nf, 3, 0.02, 20);
Cs_MSST1 = Cs_MSST/Nf*fs/2; % Normalized ridges

% The ridge extraction result of SRT
[Nf,~]=size(srt); 
Cs_SRT = brevridge_mult(abs(srt), 1:Nf, 3, 0.02, 20); 
Cs_SRT1 = Cs_SRT/Nf*fs/2; % Normalized ridges

% The ridge extraction result of the proposed method
[Nf,~]=size(tfr_pm);
Cs_PM = brevridge_mult(abs(tfr_pm), 1:Nf, 3, 0.02, 20);
Cs_PM1 = Cs_PM/Nf*fs/2; % Normalized ridges

%% Ridge-only signal reconstruction
% The reconstruction result of STFT
k=1;
j=1;
RS=zeros(3,Nt);
for k=1:3
    for j=1:Nt
        RS(k,j)=sum(real(tfr(Cs_STFT(k,j),j)));
        tfr(Cs_STFT(k,j),j)=0;
    end
end
STFT_RS=sum(real(RS)); % Reconstructed signal
% The residual and mean absolute error
signal_stft=s(:,1)'-STFT_RS(1,:); % residual error
MAE_STFT = mean(abs(signal_stft)); % mean absolute error (MAE)

% The reconstruction result of SST
k=1;
j=1;
RS=zeros(3,Nt);
for k=1:3
    for j=1:Nt
        RS(k,j)=sum(real(Ts(Cs_SST(k,j),j)));
        Ts(Cs_SST(k,j),j)=0;
    end
end
SST_RS=sum(real(RS)); % Reconstructed signal
% The residual and mean absolute error
signal_sst=s(:,1)'-SST_RS(1,:); % residual error
MAE_SST = mean(abs(signal_sst)); % mean absolute error (MAE)

% The reconstruction result of SET
k=1;
j=1;
RS=zeros(3,Nt);
for k=1:3
    for j=1:Nt
        RS(k,j)=sum(real(Te(Cs_SET(k,j),j)));
        Te(Cs_SET(k,j),j)=0;
    end
end
SET_RS=sum(real(RS)); % Reconstructed signal
% The residual and mean absolute error
signal_set=s(:,1)'-SET_RS(1,:); % residual error
MAE_SET = mean(abs(signal_set)); % mean absolute error (MAE)

% The reconstruction result of MSST
k=1;
j=1;
RS=zeros(3,Nt);
for k=1:3
    for j=1:Nt
        RS(k,j)=sum(real(MTs(Cs_MSST(k,j),j)));
        MTs(Cs_MSST(k,j),j)=0;
    end
end
MSST_RS=sum(real(RS)); % Reconstructed signal
% The residual and mean absolute error
signal_msst=s(:,1)'-MSST_RS(1,:); % residual error
MAE_MSST = mean(abs(signal_msst)); % mean absolute error (MAE)

% The reconstruction result of SRT
k=1;
j=1;
RS=zeros(3,Nt);
for k=1:3
    for j=1:Nt
        RS(k,j)=sum(real(srt(Cs_SRT(k,j),j)));
        srt(Cs_SRT(k,j),j)=0;
    end
end
SRT_RS=sum(real(RS)); % Reconstructed signal
% The residual and mean absolute error
signal_srt=s(:,1)'-SRT_RS(1,:); % residual error
MAE_srt = mean(abs(signal_srt)); % mean absolute error (MAE)

% The reconstruction result of the proposed method
k=1;
j=1;
RS=zeros(3,Nt);
for k=1:3
    for j=1:Nt
        RS(k,j)=sum(real(tfr_pm(Cs_PM(k,j),j)));
        tfr_pm(Cs_PM(k,j),j)=0;
    end
end
PM_RS=sum(real(RS)); % Reconstructed signal
% The residual and mean absolute error
signal_pm=s(:,1)'-PM_RS(1,:); % residual error
MAE_pm = mean(abs(signal_pm)); % mean absolute error (MAE)

% Residual error curve
figure
plot(time, s, 'LineWidth', 1);
hold on
plot(time, signal_stft, 'LineWidth', 1);
plot(time, signal_sst, 'LineWidth', 1);
plot(time, signal_set, 'LineWidth', 1);
plot(time, signal_msst, 'LineWidth', 1);
plot(time, signal_srt, 'LineWidth', 1);
plot(time, signal_pm, 'LineWidth', 1);
hold off
xlabel('Time (s)', 'FontSize', 10);
ylabel('Amplitude', 'FontSize', 10);
xlim([0 max(time)]);
ylim([-3 3]);
legend("Original signal",'STFT','SST','SET','MSST','SRT','PM','NumColumns', 4,'orientation', 'horizontal','FontSize', 10,"Location",'north');
set(gca, 'FontSize', 10, 'Fontname', 'Times New Roman');