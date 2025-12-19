function [tfr_pm,tfr_stft] = PM_W(s,w_length)
% =========================================================================
% This function implements the frequency-direction post-processing method
% proposed for improving the energy concentration of time-frequency representation (TFR)
% and enabling accurate signal reconstruction. The core idea is to 
% reassign the diffused TF coefficients towards the nearest local amplitude 
% maximum (ridge points) along the frequency axis, based on the amplitude gradient of STFT.
% -------------------------------------------------------------------------
% Input:
%   s        - Input signal (column vector, length N)
%   w_length - Length of the Gaussian window (default: round(N/5))
% Output:
%   tfr_pm   - Reassigned (post-processed) TFR of the proposed method [frequency × time]
%   tfr_stft - Original TFR of STFT [frequency × time]
% -------------------------------------------------------------------------
% Reference:
%   Yabo Wang et al., "A Frequency-Direction Post-Processing Method for
%   Improving Energy Concentration and Signal Reconstruction"
% Author: Yabo Wang, School of Electronic Engineering, Xidian University
% GitHub: https://github.com/wwwyyybbbcn/A-novel-frequency-direction-post-processing-method.git
% E-mail: wybcn@stu.xidian.edu.cn
% =========================================================================

% Input checking and parameter setup
[s_row,s_col] = size(s);
if (s_col~=1)
    error('S must be column vector');
end

if (nargin < 1)
    error('Please set the pending signal as the first parameter');
end

if (nargin < 2)
    w_length=round(s_row/5);
end

% Gaussian window generation (ensure odd length)
w_length=w_length+1-rem(w_length,2);
w_t = linspace(-0.5,0.5,w_length);w_t=w_t';
w = exp(-pi/0.32^2*w_t.^2);
[w_row,~]=size(w); w_hl=(w_row-1)/2;

% Compute STFT and its frequency-derivative
N_t=s_row;
t=1:s_row;
tfr1= zeros (N_t,s_row) ;
tfr2= zeros (N_t,s_row) ;
for icol=1:s_row
    ti= t(icol); tau=-min([round(N_t/2)-1,w_hl,ti-1]):min([round(N_t/2)-1,w_hl,s_row-ti]);
    index= rem(N_t+tau,N_t)+1;
    rSig = s(ti+tau,1);
    tfr1(index,icol)=rSig.*conj(w(w_hl+1+tau));
    tfr2(index,icol)=-1i*2*pi*w_t(w_hl+1+tau).*rSig.*conj(w(w_hl+1+tau));
end

% FFT, amplitude, and gradient calculation
tfr1=fft(tfr1);
tfr2=fft(tfr2);
tfr1=tfr1(1:round(N_t/2),:);
tfr2=tfr2(1:round(N_t/2),:);
% Separate real and imaginary parts
re1=real(tfr1); im1=imag(tfr1);
re2=real(tfr2); im2=imag(tfr2);
tfr_abs=abs(tfr1);
rof=tfr_abs;
% Gradient of the amplitude in frequency direction:
drof=(re2.*re1+im2.*im1)./rof;

% Compute differences of TFR and its gradient
d1_tfr=tfr_abs-[tfr_abs(2:end,:); zeros(1,N_t)];
d2_tfr=tfr_abs-[zeros(1,N_t); tfr_abs(1:end-1,:)];
d1_drof=drof-[drof(2:end,:); zeros(1,N_t)];
d2_drof=drof-[zeros(1,N_t); drof(1:end-1,:)];

% Reassignment operator calculation
% Cases (corresponding to Eq.12 in the paper):
%   Case A: local amplitude peak (ridge point) – no reassignment needed
%   Case B: positive slope – reassign upwards to nearest ridge
%   Case C: negative slope – reassign downwards to nearest ridge
omega=ones(round(N_t/2),s_row); % Reassignment operator
flag=ones(round(N_t/2),s_row); % State flag for each TF points
for j=1:N_t
    for i=1:round(N_t/2)
        if i>1 && i< round(N_t/2)
            % Case A: Detect a local amplitude peak (ridge)
            if d2_drof(i,j)<0 && d1_drof(i,j)>0 && drof(i-1,j)>0 && drof(i+1,j)<0 && d2_tfr(i,j)>0 && d1_tfr(i,j)>0
                Index2=Index1; % Previous ridge location
                Index1=i; % Current ridge location (peak)
                flag(i,j)=0;
                omega(i,j)=Index1;
                % Reassign all points between previous ridge and current ridge
                % to the current ridge (for positive slope region)
                Index3=find(flag(Index2:i,j)==1)+Index2-1;
                omega(Index3,j)=Index1;
            % Case C: Negative slope (drof < 0)
            elseif drof(i,j)<0
                flag(i,j)=-1;
                omega(i,j)=Index1;
            % Case B: Positive slope (drof > 0)
            elseif drof(i,j)>0
                flag(i,j)=1;
            end
        % Frequency lower boundary
        elseif i==1 
            Index1=1; % Reset ridge location
            if d1_tfr(i,j)>0
                flag(i,j)=1;
                omega(i,j)=Index1;
            else
                flag(i,j)=1;
            end
        % Frequency upper boundary
        else 
            if d2_tfr(i,j)>0
                Index2=Index1;
                Index1=i;
                flag(i,j)=0;
                omega(i,j)=Index1;
                Index3=find(flag(Index2:i,j)==1)+Index2-1;
                omega(Index3,j)=Index1;
            else
                flag(i,j)=-1;
                omega(i,j)=Index1;
            end
        end
    end
end

% Synchrosqueezing operation: reassign TF coefficients according to reassignment operator
tfr_pm=zeros(round(N_t/2),s_row);
thr = 0.15 * mean(rof(:,:)); % Time-adaptive threshold
for t_s=1:N_t
    for f_s=1:round(N_t/2)
        if rof(f_s,t_s) > thr(t_s)
            omega_temp = omega(f_s,t_s);
            if omega_temp>=1 && omega_temp<=round(N_t/2)
                tfr_pm(omega_temp,t_s) = tfr_pm(omega_temp,t_s) + tfr1(f_s,t_s);
            end
        end
    end
end

% Amplitude normalization
tfr_pm=tfr_pm/(s_row/2); % Normalize reassigned TFR
tfr_stft=tfr1/(sum(w)/2); % Normalize original STFT
end