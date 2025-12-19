function [SRT tt ff] = SRT(x,hlength,fs);

% input parameters:
% x denotes a discrete signal data,a row or colomn vector
% hlength denotes the window size
% fs denotes the sampling frequency of the discrete signal

% output parameters:
% SRT denotes the TFR that only includes the local maximums of the
% original-TFR
% tt denotes the discrete time of the SRT
% ff denotes the discrete frequency of the SRT

% writen by Li Miaofen at Tsinghua University, Sept. 2020. 
% any question, please contact by the email: limf18@tsinghua.org.cn

%%
[rn cn]=size(x);
if cn>=2
    x=x';
end
[xrow,xcol] = size(x);
N=xrow;
tend=N/fs; 
t=1:N;
[trow,tcol] = size(t);


%% Gaussian window

hlength=hlength+1-rem(hlength,2);
ht = linspace(-0.5,0.5,hlength);ht=ht';
h = exp(-pi/0.32.^2*ht.^2);


%% doing FFT

[hrow,hcol]=size(h); Lh=(hrow-1)/2;
tfr1= zeros (N,tcol);
tfr2= zeros (N,tcol);
for icol=1:tcol,
ti= t(icol); tau=-min([round(N/2)-1,Lh,ti-1]):min([round(N/2)-1,Lh,xrow-ti]);
indices= rem(N+tau,N)+1;
rSig = x(ti+tau,1);
tfr1(indices,icol)=rSig.*conj(h(Lh+1+tau));
tfr2(indices,icol)=-1i*2*pi*ht(Lh+1+tau).*rSig.*conj(h(Lh+1+tau)); %
end;

tfr1=fft(tfr1);
tfr2=fft(tfr2);

tfr1=tfr1(1:round(N/2),:);   re1=real(tfr1); im1=imag(tfr1);
tfr2=tfr2(1:round(N/2),:);   re2=real(tfr2); im2=imag(tfr2);
tfr11=abs(tfr1); %取模
rof=tfr11;
drof=(re2.*re1+im2.*im1)./rof; drof1=abs(drof); %

SRT=zeros(round(N/2),tcol);
fac = min([fs/N, N/fs]);
%% extracting the local maximums

for i=2:round(N/2)-1   %frequency 
    for j=1:N    %time
        %===========================================
        % three-step selecting rule for extracting the local maximums
        if drof1(i,j)<1000 & drof(i,j)<drof(i-1,j) & drof(i,j)>drof(i+1,j)
            if drof1(i,j)<drof(i-1,j) & drof1(i,j)<drof1(i+1,j)
                if tfr11(i,j)>tfr11(i-1,j) & tfr11(i,j) > tfr11(i+1,j)
                     SRT(i-floor(fac),j)=tfr1(i,j);
                end
            end
        end
        %============================================
    end
end
%% drawing the TFR

SRT=SRT/(sum(h)/2);   
[Nf Nt]=size(SRT);
tt=(1:Nt)*tend/Nt;
ff=(1:Nf)*fs/2/Nf;
figure;
imagesc(tt,ff,abs(SRT));
xlabel('Time (s)');ylabel('Freq (Hz)');
axis xy;
colormap(1-hot(128));
end