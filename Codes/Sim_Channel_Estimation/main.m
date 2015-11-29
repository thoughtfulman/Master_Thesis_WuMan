clear all;
close all;
clc;

%% %%%%%%%%%%%%%%%%%%% modulation 

% modulation handles
%256QAM
hMod256 = comm.RectangularQAMModulator('ModulationOrder',256,'BitInput',true,'NormalizationMethod',...
'Average power');
%64QAM
hMod64 = comm.RectangularQAMModulator('ModulationOrder',64,'BitInput',true,'NormalizationMethod',...
'Average power');
%16QAM 
hMod16 = comm.RectangularQAMModulator('ModulationOrder',16,'BitInput',true,'NormalizationMethod',...
'Average power');
%4QAM 
hMod4 = comm.RectangularQAMModulator('ModulationOrder',4,'BitInput',true,'NormalizationMethod',...
'Average power');
hDemod256 = comm.RectangularQAMDemodulator('ModulationOrder',256,'BitOutput',true,'NormalizationMethod',...
     'Average power');
%64QAM
hDemod64 = comm.RectangularQAMDemodulator('ModulationOrder',64,'BitOutput',true,'NormalizationMethod',...
     'Average power');
%16QAM 
hDemod16 = comm.RectangularQAMDemodulator('ModulationOrder',16,'BitOutput',true,'NormalizationMethod',...
     'Average power');
%4QAM 
hDemod4 = comm.RectangularQAMDemodulator('ModulationOrder',4,'BitOutput',true,'NormalizationMethod',...
     'Average power');
% the number of bits in an OFDM signal.
nBit_sig = sum(log2(modOrd)) * nSC_clt;
% generate the source bit
% if FixSeedEnable == 1
% rng(Seed);   % set the seed ;
% end

nFrame = 100;
source_rec = zeros(nBit_sig*nSignal,nFrame);
nError = zeros(nSC,nFrame);
for nf = 1:nFrame
    source = randi([0,1],nBit_sig*nSignal,1);

    % save source in source_ini.
    source_ini = source;  
    source = reshape(source,nBit_sig,nSignal);
    dataMod = zeros(nSC,nSignal);
    % modulation bits by bits
    for ns = 1:nSignal
        bTH =1;  % b-th bit of a signal bit vector.
        for nc = 1:nClt
            switch modOrd(nc)
                case 256
                   hMod = hMod256;
                case 64
                   hMod = hMod64;
                case 16
                   hMod = hMod16;
                case 4
                   hMod = hMod4;               
            end      
            for nt = 1:nSC_clt
                dataMod( (nc-1)*nSC_clt +nt, ns) = step(hMod,source(bTH:bTH+nModBit(nc)-1,ns));
                bTH = bTH+nModBit(nc);              
            end
        end
    end
    %% Power allocation
    dataMod = dataMod .* repmat(reshape(repmat(powerWet,nSC_clt,1),nSC,1),1,nSignal);
    %% configure virtual carriers
    dataMod(1:nIdleLF,:) =0;
    dataMod(nSC-nIdleHF+1:nSC,:) =0;

    %% Hermitian Symmetry process
    dataHS = zeros(nSC*2,nSignal);
    dataHS(1,:) = real(dataMod(1,:));
    dataHS(2:nSC,:) = dataMod(2:nSC,:);
    dataHS(nSC+1,:) = imag(dataMod(1,:));
    dataHS(nSC+2:nSC*2,:) = flipud(conj(dataMod(2:nSC,:)));

    %% IFFT
    dataTD = ifft(dataHS)*sqrt(2*nSC); 

    %% Add CP
    dataTD_cp = zeros(nCP+2*nSC,nSignal);
    dataTD_cp(1:nCP,:) = dataTD(2*nSC-nCP+1:2*nSC,:);
    dataTD_cp(nCP+1:nCP+2*nSC,:) = dataTD;

    %% reshape to a column
    dataTD_col = reshape(dataTD_cp,(nCP+2*nSC)*nSignal,1);

    %% ZC pilot
    n = 1:nSC;
    q = 0;
    r = 1;
    ZC = exp(-1i*2*pi*r*(n.*n/2+q*n)/nSC).* conj(reshape(repmat(powerWet,nSC_clt,1),nSC,1)');    % original ZC sequence of 1 cycle

    %% configure virtual carriers
    ZC(1:nIdleLF) =0;
    ZC(nSC-nIdleHF+1:nSC) =0;

    %% Generate Hermitian symmetric sequence
    ZC_freq = zeros(1,2*nSC);
    ZC_freq(1) = real(ZC(1));
    ZC_freq(2:nSC) = ZC(2:nSC);
    ZC_freq(nSC+1) = imag(ZC(1));
    ZC_freq(nSC+2:2*nSC) = conj(fliplr(ZC(2:nSC)));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Correspondent time domain sequence
    ZC_time =ifft(ZC_freq)*sqrt(2*nSC);
    % ZC_time =ZC_time/max(abs(ZC_time));

    ZC_w_CP = zeros(1,nCP+2*nSC);
    ZC_w_CP(nCP+1:nCP+2*nSC) = ZC_time;
    ZC_w_CP(1:nCP) = ZC_time(2*nSC-nCP+1:2*nSC);
    %% preCode coefficient
    preCof = 1;
    %% Add ZC
    dataTD_zc = zeros((nCP+2*nSC)*(nSignal+1),1);
    dataTD_zc(1:nCP+2*nSC) = preCof*ZC_w_CP;
    dataTD_zc(nCP+2*nSC+1:(nCP+2*nSC)*(nSignal+1)) = dataTD_cp;

    %% Pass channel
    % 1. LP channel.
    load Channel.mat;
     y_lp = conv(dataTD_zc,h_power_norm);
%     % 2. Nonlinear channel
%     bias = 4;
%     val_upclip = 8;
%     y_bias = y_lp+bias;
%     y_bias = y_bias .* (y_bias>=0 & y_bias<=val_upclip) + val_upclip*(y_bias>val_upclip);
% %     x_nl = 0:6;
% %     y_nl = [0,1,2,3,4,4.5,5];
%     x = y_bias;
%        p1 =    -0.00102  ;
%        p2 =    0.007479  ;
%        p3 =   -0.009178 ;
%        p4 =      0.9882  ;
%        p5 =    0.005828  ;
% 
%     y_nl = p1*x.^4 + p2*x.^3 + p3*x.^2 + p4*x+p5 -bias;
% %     y_nl = y_lp;
%     % 3. add noise
     SNR = 32;
     noise = randn(length(y_lp),1)/sqrt(10^(SNR/10));
%     y = y_nl + noise;

    %% demodulation 
    y = y_lp+ noise;
    rx = y';
    %256QAM


    L = length(rx);
    L_cor = L - (2*nSC) + 1;
    cor = zeros(1,L_cor);
    for i=1:L_cor
        cor(i) = rx(i:i+2*nSC-1)*ZC_time';
    end
    cor = cor/max(cor);
    threshold = 0.9;
    nShiftLeftStep = 0;
    peak_position_p = find(cor>threshold);
    peak_position = peak_position_p-nShiftLeftStep;

    ZC_rec = rx(peak_position:peak_position+2*nSC-1);
    dataTD_rec = rx(peak_position+2*nSC : peak_position+2*nSC+(2*nSC+nCP)*nSignal-1);

    %  estimate channel
    ZC_freq_rec = fft(ZC_rec) / sqrt(double(2*nSC));
    H = ZC_freq_rec./ ZC_freq;
    
    %%
    data = dataTD_rec;
    % reshape the data to remove cp
    data = reshape(data,2*nSC+nCP,nSignal);
    % remove cp
    data = data(nCP+1:2*nSC+nCP,:);
    % fft
    data_freq = fft(data) / sqrt(double(2*nSC));
    % equalization
    %% 只用1个ZC序列估计信道
    data_freq = conj(data_freq') ./ repmat(H,nSignal,1);

    % deal with the first subcarriers, in modulation,Hermitian process 
    % has put the real of the first number of data in dataMod to the first
    % subcarrier and the image put to (1+nSC)-th subcarrier. 
    % data_freq(:,1) = data_freq(:,1)+1i*data_freq(:,1+nSC); -> complex
    % I do this on purpose to make sure that data_freq(1,:) is complex
    % double, so MATLAB will not report error"Complexity mismatch with input 1; 
    % expected complex, got real." as MATLAB regard 0+0j as a real number.
    % data_freq(:,1) = 1i;
     data_freq(:,1) = data_freq(:,1)+data_freq(:,1+nSC)*1i;
    %     data_freq(:,1:nIdleLF) = complex(data_freq(:,1:nIdleLF));
    %     data_freq(:,nSC-nIdleHF+1:nSC) = complex(data_freq(:,nSC-nIdleHF+1:nSC));
    % remove the right part of frequency domain data
    data_freq = data_freq(:,1:nSC);
    % reverse it for conventient
    data_freq = conj(data_freq');
    % power allocation recovery
    data_freq = data_freq ./ repmat(reshape(repmat(powerWet,nSC_clt,1),nSC,1),1,nSignal);
    bTH = 1;  % use to label the bit order in demodulated data
    for ns = 1:nSignal        
        for nc = 1:nSC_clt
            % choose demodulation handle
              switch modOrd(nc)
                  case 256
                       hDemod = hDemod256;
                  case 64
                       hDemod = hDemod64; 
                  case 16
                       hDemod = hDemod16;
                  case 4
                       hDemod = hDemod4;
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
    H_temp(nf,:) =  H;
    ZC_freq_rec_temp(nf,:) = ZC_freq_rec;

end
%% estimate noise %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

H = mean(H_temp);
ZC_esti = ZC_freq .* H;
ZC_noise = ZC_freq_rec_temp - repmat(ZC_esti,nFrame,1);
ZC_noise_power = mean( abs(ZC_noise) .^2);
ZC_sig_power = abs(ZC_esti) .^2;
SNR = 10*log10(ZC_sig_power ./ ZC_noise_power) ;%10*log10

ModbitArray = reshape(repmat(nModBit,nSC_clt,1),nSC,1);
BER_sc = sum(nError,2)/nFrame ./double(ModbitArray)/double(nSignal);
% modOrd_sc = reshape(repmat(modOrd,nSC_clt,1),nSC,1);
BER_avg = sum(sum(nError(nIdleLF+1:nSC-nIdleHF,:))) /nFrame /double(nSignal) /...
    double(nBit_sig - sum(ModbitArray(1:nIdleLF)) - sum(ModbitArray(nSC-nIdleHF+1:nSC)));
Var_ber = var(BER_sc(nIdleLF+1:nSC-nIdleHF));
plot(SNR/max(SNR))
hold on;
plot(abs(H).^2/max(abs(H))^2,'r');
H_org = fft(h_power_norm,256);
plot(abs(H_org).^2/max(abs(H_org))^2,'k');

%% write .xlsx file
res = zeros(nSC,8);
res(:,1) = 0:nSC-1;
res(:,2) = reshape(repmat(modOrd,nSC_clt,1),nSC,1);
res(:,3) = reshape(repmat(powerWet,nSC_clt,1),nSC,1);
res(:,4) = abs(H(1:nSC));
res(:,5) = uint32(ModbitArray) * nSignal;
res(:,6) = sum(nError,2) /nFrame;
res(:,7) = SNR(1:nSC);
res(:,8) = BER_sc;
resFileName ='result_4QAM_h.xlsx';
A = {'#SubCarrier','Modulation Order(QAM)','Power weight','channel magnitude(H)','Num of total bits','Num of error bits(average)','SNR(dB)','BER_sc','BER_avg','BER_var'};
xlswrite(resFileName,A,1,'A1');
xlswrite(resFileName,res,1,'A2');
xlswrite(resFileName,BER_avg,1,'I2');
xlswrite(resFileName,Var_ber,1,'J2');