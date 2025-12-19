% =========================================================================
% Bat Echolocation Signal Analysis: A Comprehensive Comparison of 
% Time-Frequency Analysis Methods (STFT, SST, SET, MSST, SRT, and Proposed Method)
% This demonstration reproduces the bat echolocation signal experiment from the paper.
% Fig. 2, Fig. 3, Fig. 4, Fig. 5, Fig. 6, Fig. 7.
% Reference:
%   Yabo Wang et al., "A Frequency-Direction Post-Processing Method for
%   Improving Energy Concentration and Signal Reconstruction"
% Author: Yabo Wang, School of Electronic Engineering, Xidian University
% GitHub: https://github.com/wwwyyybbbcn/A-novel-frequency-direction-post-processing-method.git
% =========================================================================

%% Load and preprocess the bat echolocation signal
clear;
load('batdata1.mat'); % Load bat echolocation signal data
SampFreq = 1000000/7; % Sampling frequency (Hz)
n=length(bat(401:1460,1)); % Select a segment of the signal (1060 points)
time=(1:n)/SampFreq; % Time vector
fre=(SampFreq/2)/(n/2):(SampFreq/2)/(n/2):(SampFreq/2); % Frequency vector
time=time*1000;
fre=fre/1000;
signal=detrend(bat(401:1460,1)/max(abs(bat(:)))); % Normalize the signal

%% TFA
w_length=150;
tic;[tfr] = STFT_Y(signal,w_length);toc;
tic;[Ts] = SST_Y(signal,w_length); toc;
tic;[~,Te] = SET_Y(signal,w_length); toc;
tic;[MTs] = MSST_Y_new(signal,w_length,20);toc;
tic;[srt] = SRT(signal,w_length,SampFreq);toc;
tic;[tfr_pm] = PM_W(signal,w_length);toc;

%% Renyi entropy
renyi_stft=renyi(abs(tfr));
renyi_sst=renyi(abs(Ts));
renyi_set=renyi(abs(Te));
renyi_MSST=renyi(abs(MTs));
renyi_srt=renyi(abs(srt));
renyi_pm=renyi(abs(tfr_pm));

%% Signal waveform and STFT result Fig. 2
figureUnits = 'inches';
figureWidth = 3.5;
figureHeight = 1.35;
fontSize = 7.5;
figure('Unit', figureUnits, 'Position', [1, 1, figureWidth, figureHeight]);
leftMargin = 0.35;
rightMargin = 0.08;
bottomMargin = 0.26; 
topMargin = 0.12;
hSpace = 0.4;
plotWidth = (figureWidth - leftMargin - rightMargin - hSpace-0.16)/2;
plotHeight = figureHeight - topMargin - bottomMargin;

subplot(1,2,1);
plot(time, signal, 'LineWidth', 1);
title('(a)', 'FontSize', fontSize);
xlabel('Time (ms)', 'FontSize', fontSize);
ylabel('Amplitude', 'FontSize', fontSize);
xlim([0 max(time)]);
ylim([-1.05 1.05]);
pos1 = [leftMargin, bottomMargin, plotWidth, plotHeight];
set(gca,'Unit', figureUnits,'FontSize', fontSize, 'Fontname', 'Times New Roman', 'Position', pos1);

subplot(1,2,2);
imagesc(time,fre,abs(tfr));
xlabel('Time (ms)', 'FontSize', fontSize);
ylabel('Frequency (kHz)', 'FontSize', fontSize);
title('(b)', 'FontSize', fontSize);
axis xy;
load('MyColormap2','MyColormap2');
colormap(MyColormap2);
xlim([0 max(time)]);
ylim([0 max(fre)]);
colorbar;
colorbar('Ticks',[])
text(0.8,0.01,sprintf('R = %.2f\n', renyi_stft),'Unit', figureUnits,'FontSize', fontSize, 'Fontname', 'Times New Roman','Color', 'red','FontWeight','bold');
pos2 = [leftMargin + plotWidth + hSpace, bottomMargin, plotWidth, plotHeight];
set(gca,'Unit', figureUnits,'FontSize', fontSize, 'Fontname', 'Times New Roman', 'Position', pos2);

%% plot time-frequency representation (TFR) Fig. 3
figureUnits = 'inches';
figureWidth = 3.5;
figureHeight = 4.5;
fontSize = 7.5;
figure('Unit', figureUnits, 'Position', [1, 1, figureWidth, figureHeight]);
leftMargin = 0.27;
rightMargin = 0.04; 
bottomMargin = 0.28;
topMargin = 0.13;
hSpace = 0.35;
vSpace = 0.42;
lSpace = 0.25;
hlSpace = 0.05;
availableWidth = figureWidth - leftMargin - rightMargin;
availableHeight = figureHeight - topMargin - bottomMargin;
subWidth = (availableWidth - hSpace -2*lSpace-2*hlSpace) / 2;
subHeight = (availableHeight - 2*vSpace) / 3;
lastwidth1=1.1*(availableWidth - hSpace -0.15) / 2;
lastwidth2=0.5*(availableWidth - hSpace -0.15) / 2;
lastwidth3=0.15*(availableWidth - hSpace -0.15) / 2;
x_min = 1;
x_max = 2;
y_min = 31;
y_max = 42;

subplot(3,4,2);
imagesc(time,fre,abs(Ts));
axis xy;
load('MyColormap2','MyColormap2');
colormap(MyColormap2);
xlim([x_min x_max]);
ylim([y_min y_max]);
set(gca,'xtick',[],'xticklabel',[])
set(gca,'ytick',[],'yticklabel',[])
set(gca, 'Layer', 'top');
pos1 = [leftMargin+subWidth+hlSpace, bottomMargin+2*subHeight+2*vSpace, lSpace, subHeight];
set(gca,'Unit', figureUnits,'FontSize', fontSize, 'Fontname', 'Times New Roman', 'Position', pos1);

subplot(3,4,1);
imagesc(time,fre,abs(Ts));
xlabel('Time (ms)', 'FontSize', fontSize);
ylabel('Frequency (kHz)', 'FontSize', fontSize);
title('(a)', 'FontSize', fontSize);
axis xy;
load('MyColormap2','MyColormap2');
colormap(MyColormap2);
xlim([0 max(time)]);
ylim([0 max(fre)]);
hold on
plot([x_min, x_max, x_max, x_min, x_min], ...
      [y_min, y_min, y_max, y_max, y_min],'r-', 'LineWidth', 1);
pos1 = [leftMargin, bottomMargin+2*subHeight+2*vSpace, subWidth, subHeight];
set(gca,'Unit', figureUnits,'FontSize', fontSize, 'Fontname', 'Times New Roman', 'Position', pos1);

subplot(3,4,4);
imagesc(time,fre,abs(Ts));
axis xy;
load('MyColormap2','MyColormap2');
colormap(MyColormap2);
xlim([x_min x_max]);
ylim([y_min y_max]);
set(gca,'xtick',[],'xticklabel',[])
set(gca,'ytick',[],'yticklabel',[])
pos1 = [leftMargin + 2*subWidth + hSpace+1*lSpace+2*hlSpace, bottomMargin+2*subHeight+2*vSpace, lSpace, subHeight];
set(gca,'Unit', figureUnits,'FontSize', fontSize, 'Fontname', 'Times New Roman', 'Position', pos1);

subplot(3,4,3);
imagesc(time,fre,abs(Te));
xlabel('Time (ms)', 'FontSize', fontSize);
ylabel('Frequency (kHz)', 'FontSize', fontSize);
title('(b)', 'FontSize', fontSize);
axis xy;
load('MyColormap2','MyColormap2');
colormap(MyColormap2);
xlim([0 max(time)]);
ylim([0 max(fre)]);
hold on
plot([x_min, x_max, x_max, x_min, x_min], ...
      [y_min, y_min, y_max, y_max, y_min],'r-', 'LineWidth', 1);
pos2 = [leftMargin + subWidth + hSpace+1*lSpace+1*hlSpace, bottomMargin+2*subHeight+2*vSpace, subWidth, subHeight];
set(gca,'Unit', figureUnits,'FontSize', fontSize, 'Fontname', 'Times New Roman', 'Position', pos2);

subplot(3,4,6);
imagesc(time,fre,abs(MTs));
axis xy;
load('MyColormap2','MyColormap2');
colormap(MyColormap2);
xlim([x_min x_max]);
ylim([y_min y_max]);
set(gca,'xtick',[],'xticklabel',[])
set(gca,'ytick',[],'yticklabel',[])
pos1 = [leftMargin+1*subWidth+1*hlSpace, bottomMargin+1*subHeight+1*vSpace, lSpace, subHeight];
set(gca,'Unit', figureUnits,'FontSize', fontSize, 'Fontname', 'Times New Roman', 'Position', pos1);

subplot(3,4,5);
imagesc(time,fre,abs(MTs));
xlabel('Time (ms)', 'FontSize', fontSize);
ylabel('Frequency (kHz)', 'FontSize', fontSize);
title('(c)', 'FontSize', fontSize);
axis xy;
load('MyColormap2','MyColormap2');
colormap(MyColormap2);
xlim([0 max(time)]);
ylim([0 max(fre)]);
hold on
plot([x_min, x_max, x_max, x_min, x_min], ...
      [y_min, y_min, y_max, y_max, y_min],'r-', 'LineWidth', 1);
pos1 = [leftMargin, bottomMargin+1*subHeight+1*vSpace, subWidth, subHeight];
set(gca,'Unit', figureUnits,'FontSize', fontSize, 'Fontname', 'Times New Roman', 'Position', pos1);

subplot(3,4,8);
imagesc(time,fre,abs(srt));
axis xy;
load('MyColormap2','MyColormap2');
colormap(MyColormap2);
xlim([x_min x_max]);
ylim([y_min y_max]);
set(gca,'xtick',[],'xticklabel',[])
set(gca,'ytick',[],'yticklabel',[])
pos1 = [leftMargin + 2*subWidth + hSpace+1*lSpace+2*hlSpace,  bottomMargin+1*subHeight+1*vSpace, lSpace, subHeight];
set(gca,'Unit', figureUnits,'FontSize', fontSize, 'Fontname', 'Times New Roman', 'Position', pos1);

subplot(3,4,7);
imagesc(time,fre,abs(srt));
xlabel('Time (ms)', 'FontSize', fontSize);
ylabel('Frequency (kHz)', 'FontSize', fontSize);
title('(d)', 'FontSize', fontSize);
axis xy;
load('MyColormap2','MyColormap2');
colormap(MyColormap2);
xlim([0 max(time)]);
ylim([0 max(fre)]);
hold on
plot([x_min, x_max, x_max, x_min, x_min], ...
      [y_min, y_min, y_max, y_max, y_min],'r-', 'LineWidth', 1);
pos2 = [leftMargin + subWidth + hSpace+1*lSpace+1*hlSpace, bottomMargin+1*subHeight+1*vSpace, subWidth, subHeight];
set(gca,'Unit', figureUnits,'FontSize', fontSize, 'Fontname', 'Times New Roman', 'Position', pos2);

subplot(3,4,10);
imagesc(time,fre,abs(tfr_pm));
axis xy;
load('MyColormap2','MyColormap2');
colormap(MyColormap2);
xlim([x_min x_max]);
ylim([y_min y_max]);
c1=colorbar;
set(c1,'Unit', figureUnits, 'Position', [2.2*leftMargin+lastwidth1+0.5*hSpace+lastwidth2, bottomMargin, lastwidth3, 1.08*subHeight]);
set(gca,'xtick',[],'xticklabel',[])
set(gca,'ytick',[],'yticklabel',[])
pos2 = [2.2*leftMargin + lastwidth1 + 0.25*hSpace, 0.95*bottomMargin, lastwidth2, 1.08*subHeight];
set(gca,'Unit', figureUnits,'FontSize', fontSize, 'Fontname', 'Times New Roman', 'Position', pos2);

subplot(3,4,9);
imagesc(time,fre,abs(tfr_pm));
xlabel('Time (ms)', 'FontSize', fontSize);
ylabel('Frequency (kHz)', 'FontSize', fontSize);
title('(e)', 'FontSize', fontSize);
axis xy;
load('MyColormap2','MyColormap2');
colormap(MyColormap2);
xlim([0 max(time)]);
ylim([0 max(fre)]);
hold on
plot([x_min, x_max, x_max, x_min, x_min], ...
      [y_min, y_min, y_max, y_max, y_min],'r-', 'LineWidth', 1);
pos1 = [2.2*leftMargin, 0.95*bottomMargin, lastwidth1, 1.08*subHeight];
set(gca,'Unit', figureUnits,'FontSize', fontSize, 'Fontname', 'Times New Roman', 'Position', pos1);

%% Frequency slice analysis Fig. 4
figureUnits = 'inches';
figureWidth = 3.5; 
figureHeight = 2.45;
fontSize = 7.5;
f_slice=265;
figure('Unit', figureUnits, 'Position', [1, 1, figureWidth, figureHeight]);
leftMargin = 0.32;
rightMargin = 0.11; 
bottomMargin = 0.29;
topMargin = 0.12;
hSpace = 0.5;
vSpace = 0.4;
availableWidth = figureWidth - leftMargin - rightMargin;
availableHeight = figureHeight - topMargin - bottomMargin;
subWidth = (availableWidth - hSpace) / 2;
subHeight = (availableHeight - 1*vSpace) / 2;

subplot(2,2,1);
plot(fre,abs(tfr(:,f_slice)),"LineWidth",2,"Color","#FEB40B");
hold on;plot(fre,abs(Ts(:,f_slice)),"LineWidth",2,"Color","#443295");
plot(fre,abs(Te(:,f_slice)),"LineWidth",2,"Color","#F57F4B");
plot(fre,abs(MTs(:,f_slice)),"LineWidth",2,"Color","#518CD8");
plot(fre,abs(srt(:,f_slice)),"LineWidth",2,"Color","#994487");
plot(fre,abs(tfr_pm(:,f_slice)),"LineWidth",2,"Color","#6DC354");
hold off
xlim([0 max(fre)]);
ylim([0 0.9]);
xlabel('Frequency (kHz)', 'FontSize', fontSize);
ylabel('Amplitude', 'FontSize', fontSize);
title('(a)', 'FontSize', fontSize);
legend("STFT","SST","SET","MSST","SRT","PM");
pos1 = [leftMargin, bottomMargin+subHeight+vSpace, availableWidth, subHeight];
set(gca,'Unit', figureUnits,'FontSize', fontSize, 'Fontname', 'Times New Roman', 'Position', pos1);

subplot(2,2,3);
plot(fre,abs(tfr(:,f_slice)),"LineWidth",2,"Color","#FEB40B");
hold on;plot(fre,abs(Ts(:,f_slice)),"LineWidth",2,"Color","#443295");
plot(fre,abs(Te(:,f_slice)),"LineWidth",2,"Color","#F57F4B");
plot(fre,abs(MTs(:,f_slice)),"LineWidth",2,"Color","#518CD8");
plot(fre,abs(srt(:,f_slice)),"LineWidth",2,"Color","#994487");
plot(fre,abs(tfr_pm(:,f_slice)),"LineWidth",2,"Color","#6DC354");
hold off
xlim([32.2 33.4]);
ylim([0 0.9]);
xlabel('Frequency (kHz)', 'FontSize', fontSize);
ylabel('Amplitude', 'FontSize', fontSize);
title('(b)', 'FontSize', fontSize);
pos1 = [leftMargin, bottomMargin, subWidth, subHeight];
set(gca,'Unit', figureUnits,'FontSize', fontSize, 'Fontname', 'Times New Roman', 'Position', pos1);

subplot(2,2,4);
plot(fre,abs(tfr(:,f_slice)),"LineWidth",2,"Color","#FEB40B");
hold on;plot(fre,abs(Ts(:,f_slice)),"LineWidth",2,"Color","#443295");
plot(fre,abs(Te(:,f_slice)),"LineWidth",2,"Color","#F57F4B");
plot(fre,abs(MTs(:,f_slice)),"LineWidth",2,"Color","#518CD8");
plot(fre,abs(srt(:,f_slice)),"LineWidth",2,"Color","#994487");
plot(fre,abs(tfr_pm(:,f_slice)),"LineWidth",2,"Color","#6DC354");
hold off
xlim([62.6 64]);
ylim([0 0.035]);
xlabel('Frequency (kHz)', 'FontSize', fontSize);
ylabel('Amplitude', 'FontSize', fontSize);
title('(c)', 'FontSize', fontSize);
pos2 = [leftMargin + subWidth + hSpace, bottomMargin, subWidth, subHeight];
set(gca,'Unit', figureUnits,'FontSize', fontSize, 'Fontname', 'Times New Roman', 'Position', pos2);

%% The 3D view of TFR Fig. 5
% The ridge extraction result of STFT
[Nf,Nt]=size(tfr);
Cs_STFT = brevridge_mult(abs(tfr), 1:Nf, 2, 0.02, 80);
Cs_STFT1 = Cs_STFT/Nf*SampFreq/2/1000;

% The ridge extraction result of SST
[Nf,Nt]=size(Ts);
Cs_SST = brevridge_mult(abs(Ts), 1:Nf, 2, 0.02, 20);
Cs_SST1 = Cs_SST/Nf*SampFreq/2/1000;

% The ridge extraction result of SET
[Nf,Nt]=size(Te);
Cs_SET = brevridge_mult(abs(Te), 1:Nf, 2, 0.02, 30);
Cs_SET1 = Cs_SET/Nf*SampFreq/2/1000;

% The ridge extraction result of MSST
[Nf,Nt]=size(MTs);
Cs_MSST = brevridge_mult(abs(MTs), 1:Nf, 2, 0.02, 20);
Cs_MSST1 = Cs_MSST/Nf*SampFreq/2/1000;

% The ridge extraction result of SRT
[Nf,Nt]=size(srt);
Cs_SRT = brevridge_mult(abs(srt), 1:Nf, 2, 0.02, 80);
Cs_SRT1 = Cs_SRT/Nf*SampFreq/2/1000;

% The ridge extraction result of the proposed method
[Nf,Nt]=size(tfr_pm);
Cs_PM = brevridge_mult(abs(tfr_pm), 1:Nf, 2, 0.02, 80);
Cs_PM1 = Cs_PM/Nf*SampFreq/2/1000;

% Plot the 3D view of TFR
figureUnits = 'inches';
figureWidth = 3.5;
figureHeight = 4;
fontSize = 7.5;
figure('Unit', figureUnits, 'Position', [1, 1, figureWidth, figureHeight]);
leftMargin = 0.45;
rightMargin = 0.04;
bottomMargin = 0.23;
topMargin = 0.13;
hSpace = 0.5;
vSpace = 0.35;
availableWidth = figureWidth - leftMargin - rightMargin;
availableHeight = figureHeight - topMargin - bottomMargin;
subWidth = (availableWidth - hSpace) / 2;
subHeight = (availableHeight - 2*vSpace) / 3;

subplot(3,2,1);
mesh(time,fre,abs(tfr),'LineWidth', 0.5);
xlabel('Time (ms)', 'FontSize', fontSize);
ylabel('Frequency (kHz)', 'FontSize', fontSize,"Rotation",-65);
zlabel('Amplitude', 'FontSize', fontSize);
title('(a)', 'FontSize', fontSize);
axis xy;
hold on
plot3(time(1:1060),Cs_STFT1(1,1:1060),abs(tfr(sub2ind(size(tfr),Cs_STFT(1,1:1060),1:1060))),'color',"r",'linewidth',1,'LineStyle','-');
plot3(time(174:1060),Cs_STFT1(2,174:1060),abs(tfr(sub2ind(size(tfr),Cs_STFT(2,174:1060),174:1060))),'color',"r",'linewidth',1,'LineStyle','-');
xlim([0 max(time)]);
ylim([0 max(fre)]);
zlim([0 0.95]);
view(-19,72)
xticks([0,2,4,6]);
yticks([0,20,40,60]);
shading interp;
set(gcf, 'Color', 'white');
pos1 = [leftMargin, bottomMargin+2*subHeight+2*vSpace, subWidth, subHeight];
set(gca,'Unit', figureUnits,'FontSize', fontSize, 'Fontname', 'Times New Roman', 'Position', pos1);

subplot(3,2,2);
mesh(time,fre,abs(Ts),'LineWidth', 0.5);
xlabel('Time (ms)', 'FontSize', fontSize);
ylabel('Frequency (kHz)', 'FontSize', fontSize,"Rotation",-65);
zlabel('Amplitude', 'FontSize', fontSize);
title('(b)', 'FontSize', fontSize);
axis xy;
hold on
plot3(time(1:1060),Cs_SST1(1,1:1060),abs(Ts(sub2ind(size(Ts),Cs_SST(1,1:1060),1:1060))),'color',"r",'linewidth',1,'LineStyle','-');
plot3(time(174:1060),Cs_SST1(2,174:1060),abs(Ts(sub2ind(size(Ts),Cs_SST(2,174:1060),174:1060))),'color',"r",'linewidth',1,'LineStyle','-');
xlim([0 max(time)]);
ylim([0 max(fre)]);
zlim([0 0.95]);
view(-19,72)
xticks([0,2,4,6]);
yticks([0,20,40,60]);
shading interp;
set(gcf, 'Color', 'white');
pos2 = [leftMargin + subWidth + hSpace, bottomMargin+2*subHeight+2*vSpace, subWidth, subHeight];
set(gca,'Unit', figureUnits,'FontSize', fontSize, 'Fontname', 'Times New Roman', 'Position', pos2);

subplot(3,2,3);
mesh(time,fre,abs(Te),'LineWidth', 0.5);
xlabel('Time (ms)', 'FontSize', fontSize);
ylabel('Frequency (kHz)', 'FontSize', fontSize,"Rotation",-65);
zlabel('Amplitude', 'FontSize', fontSize);
title('(c)', 'FontSize', fontSize);
axis xy;
hold on
plot3(time(1:1060),Cs_SET1(1,1:1060),abs(Te(sub2ind(size(Te),Cs_SET(1,1:1060),1:1060))),'color',"r",'linewidth',1,'LineStyle','-');
plot3(time(174:1060),Cs_SET1(2,174:1060),abs(Te(sub2ind(size(Te),Cs_SET(2,174:1060),174:1060))),'color',"r",'linewidth',1,'LineStyle','-');
xlim([0 max(time)]);
ylim([0 max(fre)]);
zlim([0 0.95]);
view(-19,72)
xticks([0,2,4,6]);
yticks([0,20,40,60]);
shading interp;
set(gcf, 'Color', 'white');
pos1 = [leftMargin, bottomMargin+1*subHeight+1*vSpace, subWidth, subHeight];
set(gca,'Unit', figureUnits,'FontSize', fontSize, 'Fontname', 'Times New Roman', 'Position', pos1);

subplot(3,2,4);
mesh(time,fre,abs(MTs),'LineWidth', 0.5);
xlabel('Time (ms)', 'FontSize', fontSize);
ylabel('Frequency (kHz)', 'FontSize', fontSize,"Rotation",-65);
zlabel('Amplitude', 'FontSize', fontSize);
title('(d)', 'FontSize', fontSize);
axis xy;
hold on
plot3(time(1:1060),Cs_MSST1(1,1:1060),abs(MTs(sub2ind(size(MTs),Cs_MSST(1,1:1060),1:1060))),'color',"r",'linewidth',1,'LineStyle','-');
plot3(time(174:1060),Cs_MSST1(2,174:1060),abs(MTs(sub2ind(size(MTs),Cs_MSST(2,174:1060),174:1060))),'color',"r",'linewidth',1,'LineStyle','-');
xlim([0 max(time)]);
ylim([0 max(fre)]);
zlim([0 0.95]);
view(-19,72)
xticks([0,2,4,6]);
yticks([0,20,40,60]);
shading interp;
set(gcf, 'Color', 'white');
pos2 = [leftMargin + subWidth + hSpace, bottomMargin+1*subHeight+1*vSpace, subWidth, subHeight];
set(gca,'Unit', figureUnits,'FontSize', fontSize, 'Fontname', 'Times New Roman', 'Position', pos2);

subplot(3,2,5);
mesh(time,fre,abs(srt),'LineWidth', 0.5);
xlabel('Time (ms)', 'FontSize', fontSize);
ylabel('Frequency (kHz)', 'FontSize', fontSize,"Rotation",-65);
zlabel('Amplitude', 'FontSize', fontSize);
title('(e)', 'FontSize', fontSize);
axis xy;
hold on
plot3(time(1:1060),Cs_SRT1(1,1:1060),abs(srt(sub2ind(size(srt),Cs_SRT(1,1:1060),1:1060))),'color',"r",'linewidth',1,'LineStyle','-');
plot3(time(174:1060),Cs_SRT1(2,174:1060),abs(srt(sub2ind(size(srt),Cs_SRT(2,174:1060),174:1060))),'color',"r",'linewidth',1,'LineStyle','-');
xlim([0 max(time)]);
ylim([0 max(fre)]);
zlim([0 0.95]);
view(-19,72)
xticks([0,2,4,6]);
yticks([0,20,40,60]);
shading interp;
set(gcf, 'Color', 'white');
pos1 = [leftMargin, bottomMargin, subWidth, subHeight];
set(gca,'Unit', figureUnits,'FontSize', fontSize, 'Fontname', 'Times New Roman', 'Position', pos1);

subplot(3,2,6);
mesh(time,fre,abs(tfr_pm),'LineWidth', 0.5);
xlabel('Time (ms)', 'FontSize', fontSize);
ylabel('Frequency (kHz)', 'FontSize', fontSize,"Rotation",-65);
zlabel('Amplitude', 'FontSize', fontSize);
title('(f)', 'FontSize', fontSize);
axis xy;
hold on
plot3(time(1:1060),Cs_PM1(1,1:1060),abs(tfr_pm(sub2ind(size(tfr_pm),Cs_PM(1,1:1060),1:1060))),'color',"r",'linewidth',1,'LineStyle','-');
plot3(time(174:1060),Cs_PM1(2,174:1060),abs(tfr_pm(sub2ind(size(tfr_pm),Cs_PM(2,174:1060),174:1060))),'color',"r",'linewidth',1,'LineStyle','-');
xlim([0 max(time)]);
ylim([0 max(fre)]);
zlim([0 0.95]);
view(-19,72)
xticks([0,2,4,6]);
yticks([0,20,40,60]);
shading interp;
set(gcf, 'Color', 'white');
pos2 = [leftMargin + subWidth + hSpace, bottomMargin, subWidth, subHeight];
set(gca,'Unit', figureUnits,'FontSize', fontSize, 'Fontname', 'Times New Roman', 'Position', pos2);

%% Ridge-only signal reconstruction
% The reconstruction result of STFT
k=1;
j=1;
RS=zeros(2,Nt);
for k=1:2
    for j=1:Nt
        RS(k,j)=sum(real(tfr(Cs_STFT(k,j),j)));
        tfr(Cs_STFT(k,j),j)=0;
    end
end
STFT_RS=sum(real(RS)); % Reconstructed signal
% The residual and mean absolute error
signal_stft=signal(:,1)'-STFT_RS(1,:); % residual error
MAE_STFT = mean(abs(signal_stft)); % mean absolute error (MAE)

% The reconstruction result of SST
k=1;
j=1;
RS=zeros(2,Nt);
for k=1:2
    for j=1:Nt
        RS(k,j)=sum(real(Ts(Cs_SST(k,j),j)));
        Ts(Cs_SST(k,j),j)=0;
    end
end
SST_RS=sum(real(RS)); % Reconstructed signal
% The residual and mean absolute error
signal_sst=signal(:,1)'-SST_RS(1,:); % residual error
MAE_SST = mean(abs(signal_sst)); % mean absolute error (MAE)

% The reconstruction result of SET
k=1;
j=1;
RS=zeros(2,Nt);
for k=1:2
    for j=1:Nt
        RS(k,j)=sum(real(Te(Cs_SET(k,j),j)));
        Te(Cs_SET(k,j),j)=0;
    end
end
SET_RS=sum(real(RS)); % Reconstructed signal
% The residual and mean absolute error
signal_set=signal(:,1)'-SET_RS(1,:); % residual error
MAE_SET = mean(abs(signal_set)); % mean absolute error (MAE)

% The reconstruction result of MSST
k=1;
j=1;
RS=zeros(2,Nt);
for k=1:2
    for j=1:Nt
        RS(k,j)=sum(real(MTs(Cs_MSST(k,j),j)));
        MTs(Cs_MSST(k,j),j)=0;
    end
end
MSST_RS=sum(real(RS)); % Reconstructed signal
% The residual and mean absolute error
signal_msst=signal(:,1)'-MSST_RS(1,:); % residual error
MAE_MSST = mean(abs(signal_msst)); % mean absolute error (MAE)

% The reconstruction result of SRT
k=1;
j=1;
RS=zeros(2,Nt);
for k=1:2
    for j=1:Nt
        RS(k,j)=sum(real(srt(Cs_SRT(k,j),j)));
        srt(Cs_SRT(k,j),j)=0;
    end
end
SRT_RS=sum(real(RS)); % Reconstructed signal
% The residual and mean absolute error
signal_srt=signal(:,1)'-SRT_RS(1,:); % residual error
MAE_srt = mean(abs(signal_srt)); % mean absolute error (MAE)

% The reconstruction result of the proposed method
k=1;
j=1;
RS=zeros(2,Nt);
for k=1:2
    for j=1:Nt
        RS(k,j)=sum(real(tfr_pm(Cs_PM(k,j),j)));
        tfr_pm(Cs_PM(k,j),j)=0;
    end
end
PM_RS=sum(real(RS)); % Reconstructed signal
% The residual and mean absolute error
signal_pm=signal(:,1)'-PM_RS(1,:); % residual error
MAE_pm = mean(abs(signal_pm)); % mean absolute error (MAE)

% Plot signal reconstruction result Fig. 6
figureUnits = 'inches';
figureWidth = 3.5;
figureHeight = 2.5;
fontSize = 7.5;
figure('Unit', figureUnits, 'Position', [1, 1, figureWidth, figureHeight]);
leftMargin = 0.25;
rightMargin = 0.03;
bottomMargin = 0.27;
topMargin = 0.125;
hSpace = 0.36;
vSpace = 0.42;
availableWidth = figureWidth - leftMargin - rightMargin;
availableHeight = figureHeight - topMargin - bottomMargin;
subWidth = (availableWidth - hSpace) / 2;
subHeight = (availableHeight - 2*vSpace) / 3; 

subplot(3,2,1);
plot(time, signal, 'LineWidth', 1);
hold on
plot(time, STFT_RS, 'LineWidth', 1,"Color","#443295");
xlabel('Time (ms)', 'FontSize', fontSize);
ylabel('Amplitude', 'FontSize', fontSize);
title('(a)', 'FontSize', fontSize);
xlim([0 max(time)]);
ylim([-1.05 1.05]);
text(0.75,0.01,sprintf('MAE = %.4f\n', MAE_STFT),'Unit', figureUnits,'FontSize', fontSize, 'Fontname', 'Times New Roman','Color', 'red','FontWeight','bold');
pos1 = [leftMargin, bottomMargin+2*subHeight+2*vSpace, subWidth, subHeight];
set(gca,'Unit', figureUnits,'FontSize', fontSize, 'Fontname', 'Times New Roman', 'Position', pos1);

subplot(3,2,2);
plot(time, signal, 'LineWidth', 1);
hold on
plot(time, SST_RS, 'LineWidth', 1,"Color","#443295");
xlabel('Time (ms)', 'FontSize', fontSize);
ylabel('Amplitude', 'FontSize', fontSize);
title('(b)', 'FontSize', fontSize);
xlim([0 max(time)]);
ylim([-1.05 1.05]);
text(0.75,0.01,sprintf('MAE = %.4f\n', MAE_SST),'Unit', figureUnits,'FontSize', fontSize, 'Fontname', 'Times New Roman','Color', 'red','FontWeight','bold');
pos2 = [leftMargin + subWidth + hSpace, bottomMargin+2*subHeight+2*vSpace, subWidth, subHeight];
set(gca,'Unit', figureUnits,'FontSize', fontSize, 'Fontname', 'Times New Roman', 'Position', pos2);

subplot(3,2,3);
plot(time, signal, 'LineWidth', 1);
hold on
plot(time, SET_RS, 'LineWidth', 1,"Color","#443295");
xlabel('Time (ms)', 'FontSize', fontSize);
ylabel('Amplitude', 'FontSize', fontSize);
title('(c)', 'FontSize', fontSize);
xlim([0 max(time)]);
ylim([-1.05 1.05]);
text(0.75,0.01,sprintf('MAE = %.4f\n', MAE_SET),'Unit', figureUnits,'FontSize', fontSize, 'Fontname', 'Times New Roman','Color', 'red','FontWeight','bold');
pos1 = [leftMargin, bottomMargin+1*subHeight+1*vSpace, subWidth, subHeight];
set(gca,'Unit', figureUnits,'FontSize', fontSize, 'Fontname', 'Times New Roman', 'Position', pos1);

subplot(3,2,4);
plot(time, signal, 'LineWidth', 1);
hold on
plot(time, MSST_RS, 'LineWidth', 1,"Color","#443295");
xlabel('Time (ms)', 'FontSize', fontSize);
ylabel('Amplitude', 'FontSize', fontSize);
title('(d)', 'FontSize', fontSize);
xlim([0 max(time)]);
ylim([-1.05 1.05]);
text(0.75,0.01,sprintf('MAE = %.4f\n', MAE_MSST),'Unit', figureUnits,'FontSize', fontSize, 'Fontname', 'Times New Roman','Color', 'red','FontWeight','bold');
pos2 = [leftMargin + subWidth + hSpace, bottomMargin+1*subHeight+1*vSpace, subWidth, subHeight];
set(gca,'Unit', figureUnits,'FontSize', fontSize, 'Fontname', 'Times New Roman', 'Position', pos2);

subplot(3,2,5);
plot(time, signal, 'LineWidth', 1);
hold on
plot(time, SRT_RS, 'LineWidth', 1,"Color","#443295");
xlabel('Time (ms)', 'FontSize', fontSize);
ylabel('Amplitude', 'FontSize', fontSize);
title('(e)', 'FontSize', fontSize);
xlim([0 max(time)]);
ylim([-1.05 1.05]);
text(0.75,0.01,sprintf('MAE = %.4f\n', MAE_srt),'Unit', figureUnits,'FontSize', fontSize, 'Fontname', 'Times New Roman','Color', 'red','FontWeight','bold');
pos1 = [leftMargin, bottomMargin, subWidth, subHeight];
set(gca,'Unit', figureUnits,'FontSize', fontSize, 'Fontname', 'Times New Roman', 'Position', pos1);

subplot(3,2,6);
plot(time, signal, 'LineWidth', 1);
hold on
plot(time, PM_RS, 'LineWidth', 1,"Color","#443295");
xlabel('Time (ms)', 'FontSize', fontSize);
ylabel('Amplitude', 'FontSize', fontSize);
title('(f)', 'FontSize', fontSize);
xlim([0 max(time)]);
ylim([-1.05 1.05]);
text(0.75,0.01,sprintf('MAE = %.4f\n', MAE_pm),'Unit', figureUnits,'FontSize', fontSize, 'Fontname', 'Times New Roman','Color', 'red','FontWeight','bold');
pos2 = [leftMargin + subWidth + hSpace, bottomMargin, subWidth, subHeight];
set(gca,'Unit', figureUnits,'FontSize', fontSize, 'Fontname', 'Times New Roman', 'Position', pos2);

% Residual error curve Fig. 7
figureUnits = 'inches';
figureWidth = 3.5;
figureHeight = 1.35;
fontSize = 7.5;
figure('Unit', figureUnits, 'Position', [1, 1, figureWidth, figureHeight]);
leftMargin = 0.32;
rightMargin = 0.05;
bottomMargin = 0.27;
topMargin = 0.15;
plotWidth = (figureWidth - leftMargin - rightMargin);
plotHeight = figureHeight - topMargin - bottomMargin;
plot(time, signal, 'LineWidth', 1);
hold on
plot(time, signal_stft, 'LineWidth', 1,"Color","#FEB40B");
plot(time, signal_sst, 'LineWidth', 1,"Color","#443295");
plot(time, signal_set, 'LineWidth', 1,"Color","#F57F4B");
plot(time, signal_msst, 'LineWidth', 1,"Color","#518CD8");
plot(time, signal_srt, 'LineWidth', 1,"Color","#994487");
plot(time, signal_pm, 'LineWidth', 1,"Color","#6DC354");
hold off
xlabel('Time (ms)', 'FontSize', fontSize);
ylabel('Amplitude', 'FontSize', fontSize);
xlim([0 max(time)]);
ylim([-1.05 1.05]);
legend("Input signal",'STFT','SST','SET','MSST','SRT','PM','NumColumns', 4,'orientation', 'horizontal','FontSize', 5.5,"Location",'north');
pos1 = [leftMargin, bottomMargin, plotWidth, plotHeight];
set(gca,'Unit', figureUnits,'FontSize', fontSize, 'Fontname', 'Times New Roman', 'Position', pos1);