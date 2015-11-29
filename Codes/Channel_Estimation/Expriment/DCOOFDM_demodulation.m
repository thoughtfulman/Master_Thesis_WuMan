%% program statement %%%%%%%%%%%%%%%%%%%%%%%
% This program is for DCO-OFDM demodulation, supplies various number of 
% subcarriers(nSC), number of subcarriers in cluster(nSC_clt), number of CP
%(CP), modulation order on each cluster(modOrd) and cluster power weight.
% allocation(powerWet).
% all arguments are defined in the file 'modulation_initial_file.txt'

% author: Man Wu
% build time: 2014-7-31
% revised by: Man Wu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clear all;
% close all;
% clc;
%% load files
% xlLoadChipScopeData('200M2.prn');
load(Output_file_name_mat);
rx = ADC4_D';

%% demodulation handles
% modulation handles
%256QAM
hDemod256 = comm.RectangularQAMDemodulator('ModulationOrder',256,'BitOutput',true,'NormalizationMethod',...
     'Average power');
%128QAM
hDemod128 = comm.RectangularQAMDemodulator('ModulationOrder',128,'BitOutput',true,'NormalizationMethod',...
     'Average power');
%64QAM 
hDemod64 = comm.RectangularQAMDemodulator('ModulationOrder',64,'BitOutput',true,'NormalizationMethod',...
     'Average power');
%32QAM 
hDemod32 = comm.RectangularQAMDemodulator('ModulationOrder',32,'BitOutput',true,'NormalizationMethod',...
     'Average power');
%16QAM
hDemod16 = comm.RectangularQAMDemodulator('ModulationOrder',16,'BitOutput',true,'NormalizationMethod',...
     'Average power');
%8QAM
hDemod8 = comm.RectangularQAMDemodulator('ModulationOrder',8,'BitOutput',true,'NormalizationMethod',...
     'Average power');
%16QAM 
hDemod4 = comm.RectangularQAMDemodulator('ModulationOrder',4,'BitOutput',true,'NormalizationMethod',...
     'Average power');
%4QAM 
hDemod2 = comm.BPSKDemodulator;


% calculate correlation
L = length(rx);
L_cor = L - (2*nSC) + 1;
cor = zeros(1,L_cor);
for i=1:L_cor
    cor(i) = rx(i:i+2*nSC-1)*ZC_time';
end
cor = cor/max(cor);

% find peak position （head of ZC）
threshold = 0.9;
peak_find = cor>threshold;
peak_position_p = find(peak_find);
pp = 2;
peak_position(1) = peak_position_p(1);
for p = 2:length(peak_position_p)
    if peak_position_p(p) == peak_position_p(p-1)+1
        peak_position_p(p) = 0;
    else
        peak_position(pp) = peak_position_p(p);
        pp = pp+1;
    end
end
% move to left some positions
peak_position = peak_position - nShiftLeftStep;
%% receive data process
nFrame_rec = length(peak_position) - 1;
dataTD_rec = zeros((2*nSC+nCP)*nSignal,nFrame_rec);
ZC_rec = zeros(2*nSC,nFrame_rec);

for i=1:nFrame_rec
ZC_rec(:,i) = rx(peak_position(i):peak_position(i)+2*nSC-1);
dataTD_rec(:,i) = rx(peak_position(i)+2*nSC : peak_position(i)+2*nSC+(2*nSC+nCP)*nSignal-1);
end


%  estimate channel
ZC_freq_rec = fft(ZC_rec) / sqrt(double(2*nSC));
% ZC_freq = ZC_freq ./ [conj(reshape(repmat(powerWet,nSC_clt,1),nSC,1)'),fliplr((reshape(repmat(powerWet,nSC_clt,1),nSC,1)'))];
% ZC_freq = ;
H = conj(ZC_freq_rec') ./ repmat(ZC_freq,nFrame_rec,1);
H_term = H;
if nFrame_rec>1
  H = sum(H) / nFrame_rec;
end

%% estimate noise %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ZC_esti = ZC_freq .* H;
ZC_noise = conj(ZC_freq_rec') - repmat(ZC_esti,nFrame_rec,1);
ZC_noise_power = mean( abs(ZC_noise) .^2);
ZC_sig_power = abs(ZC_esti) .^2;
SNR = 10*log10(ZC_sig_power ./ ZC_noise_power) ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nBit_sig = uint32(nBit_sig);
% equlization and demodulation process frame by frame
source_rec = zeros(nBit_sig*nSignal,nFrame_rec);
nError = zeros(nSC,nFrame_rec);
for i=1:nFrame_rec
    data = dataTD_rec(:,i);
    % reshape the data to remove cp
    data = reshape(data,2*nSC+nCP,nSignal);
    % remove cp
    data = data(nCP+1:2*nSC+nCP,:);
    % fft
    data_freq = fft(data) / sqrt(double(2*nSC));
    % equalization
    %% 只用1个ZC序列估计信道
%     H = H_term(i,:);
    data_freq = conj(data_freq') ./ repmat(H,nSignal,1);

    % deal with the first subcarriers, in modulation,Hermitian process 
    % has put the real of the first number of data in dataMod to the first
    % subcarrier and the image put to (1+nSC)-th subcarrier. 
    % data_freq(:,1) = data_freq(:,1)+1i*data_freq(:,1+nSC); -> complex
    % I do this on purpose to make sure that data_freq(1,:) is complex
    % double, so MATLAB will not report error"Complexity mismatch with input 1; 
    % expected complex, got real." as MATLAB regard 0+0j as a real number.
%      data_freq(:,1) = 1i;
    data_freq(:,1:nIdleLF) = 1i;
    data_freq(:,nSC-nIdleHF+1:nSC) = 1i;
    % remove the right part of frequency domain data
    data_freq = data_freq(:,1:nSC);
    % reverse it for conventient
    data_freq = conj(data_freq');
    % power allocation recovery
    error_vector = data_freq - dataMod;
    error_vector_mag(:,i) = mean(abs(error_vector).^2,2);
    data_freq = data_freq ./ repmat(powerWet,1,nSignal);
    bTH = 1;  % use to label the bit order in demodulated data
    for ns = 1:nSignal        
        for nc = nIdleLF+1:nSC-nIdleHF;
            % choose demodulation handle
              switch modOrd(nc)
                  case 256
                       hDemod = hDemod256;
                  case 128
                       hDemod = hDemod128;
                  case 64
                       hDemod = hDemod64;
                  case 32
                       hDemod = hDemod32;   
                  case 16
                       hDemod = hDemod16;
                  case 8
                       hDemod = hDemod8;
                  case 4
                       hDemod = hDemod4;
                  case 2
                        hDemod = hDemod2; 
                  otherwise
                        continue;
              end   
              
              % demodulation
             source_rec(bTH:bTH+nModBit(nc)-1,i) = step(hDemod,data_freq(nc,ns));
              % count errors
             nError(nc,i) = nError(nc,i) + sum(abs(source_rec(bTH:bTH+nModBit(nc)-1,i) - source_ini(bTH:bTH+nModBit(nc)-1)));
             bTH = bTH+nModBit(nc);                 
        end
    end
    
end
SNR_evm = 10*log10(powerWet.^2  ./ mean(error_vector_mag(:,i),2));
ModbitArray = nModBit;
BER_sc= sum(nError,2)/nFrame_rec ./double(ModbitArray)/double(nSignal);
BER_sc(1:nIdleLF)=1;
BER_sc(nSC-nIdleHF+1:nSC)=1;
% modOrd_sc = reshape(repmat(modOrd,nSC_clt,1),nSC,1);
BER_avg = sum(sum(nError(nIdleLF+1:nSC-nIdleHF,:))) /nFrame_rec /double(nSignal) /...
    double(nBit_sig - sum(ModbitArray(1:nIdleLF)) - sum(ModbitArray(nSC-nIdleHF+1:nSC)));
Var_ber = var(BER_sc(nIdleLF+1:nSC-nIdleHF));

%% write .xlsx file
res = zeros(nSC,9);
res(:,1) = 0:nSC-1;
res(:,2) = modOrd;
res(:,3) = powerWet;
res(:,4) = abs(H(1:nSC));
res(:,5) = ModbitArray * nSignal;
res(:,6) = sum(nError,2) /nFrame_rec;
res(:,7) = SNR(1:nSC);
res(:,8) = SNR_evm(1:nSC);
res(:,9) = BER_sc;

A = {'#SubCarrier','Modulation Order(QAM)','Power weight','channel magnitude(H)','Num of total bits','Num of error bits(average)','SNR(dB)','SNR_evm(dB)','BER_sc','BER_avg','BER_var'};
xlswrite(resFileName,A,1,'A1');
xlswrite(resFileName,res,1,'A2');
xlswrite(resFileName,BER_avg,1,'J2');
xlswrite(resFileName,Var_ber,1,'K2');



 