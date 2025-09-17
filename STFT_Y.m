function [tfr] = STFT_Y(x,hlength);
%   Short-time Fourier transform (STFT)
%	x       : Signal.
%	hlength : Window length.

%	IF   : Instantaneous frequency representation.
%   Te   : SET result.
%   tfr  : STFT result
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
%
%   Written by YuGang in Shandong University at 2016.5.13.

[xrow,xcol] = size(x);

N=xrow;

if (xcol~=1),
 error('X must be column vector');
end;

if (nargin < 2),
 hlength=round(xrow/8);
end;

t=1:N;
ft = 1:round(N/2);

[trow,tcol] = size(t);

hlength=hlength+1-rem(hlength,2);
ht = linspace(-0.5,0.5,hlength);ht=ht';

% Gaussian window
h = exp(-pi/0.32^2*ht.^2);


[hrow,hcol]=size(h); Lh=(hrow-1)/2;

tfr1= zeros (N,tcol);

for icol=1:tcol,
ti= t(icol); tau=-min([round(N/2)-1,Lh,ti-1]):min([round(N/2)-1,Lh,xrow-ti]);
indices= rem(N+tau,N)+1;
rSig = x(ti+tau,1);
tfr1(indices,icol)=rSig.*conj(h(Lh+1+tau));
end;

tfr1=fft(tfr1);
tfr1=tfr1(1:round(N/2),:);
% tfr1=tfr1/(xrow/2);%the amplitude of tfr result has been pre-rectified.
tfr=tfr1/(sum(h)/2);
end