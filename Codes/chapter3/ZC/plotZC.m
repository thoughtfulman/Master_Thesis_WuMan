clear all;close all;clc;
nSC = 128;
n = [1:nSC]';
q = 0;
r = 1;
ZC = exp(-1i*2*pi*r*(n.*n/2+q*n)/nSC);%.* conj(repmat(powerWet,1,nSignal)');    % original ZC sequence of 1 cycle
%% Generate Hermitian symmetric sequence
ZC_freq = zeros(1,2*nSC);
ZC_freq(1) = real(ZC(1));
ZC_freq(2:nSC) = ZC(2:nSC);
ZC_freq(nSC+1) = imag(ZC(1));
ZC_freq(nSC+2:2*nSC) = conj(flipud(ZC(2:nSC)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Correspondent time domain sequence
ZC_time =ifft(ZC_freq)*sqrt(2*nSC);
% ZC_time =ZC_time/max(abs(ZC_time));
plot(0:2*nSC-1,abs(ZC_freq),'LineWidth',1.3);
axis([0,2*nSC-1,0,1.2]);
xlabel('Subcarrier Index');
ylabel('Magnitude');
grid on;
figure(2);
plot(0:2*nSC-1,ZC_time,'LineWidth',1.2);
axis([0,2*nSC-1,-2,2]);
xlabel('Subcarrier Index');
ylabel('Magnitude');
grid on;