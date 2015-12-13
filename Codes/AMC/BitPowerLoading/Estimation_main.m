%% loading parameters
clear all;close all; clc;
%[M_LMMSE,M_SVD] = calc_Matrix();
ReadInitialFile;
% modulation handles
%256QAM
hMod256 = comm.RectangularQAMModulator('ModulationOrder',256,'BitInput',true,'NormalizationMethod',...
   'Average power');
%64QAM
hMod128 = comm.RectangularQAMModulator('ModulationOrder',128,'BitInput',true,'NormalizationMethod',...
   'Average power');
%64QAM
hMod64 = comm.RectangularQAMModulator('ModulationOrder',64,'BitInput',true,'NormalizationMethod',...
   'Average power');
%32QAM
hMod32 = comm.RectangularQAMModulator('ModulationOrder',32,'BitInput',true,'NormalizationMethod',...
   'Average power');
%16QAM 
hMod16 = comm.RectangularQAMModulator('ModulationOrder',16,'BitInput',true,'NormalizationMethod',...
   'Average power');
%4QAM
hMod8 = comm.RectangularQAMModulator('ModulationOrder',8,'BitInput',true,'NormalizationMethod',...
   'Average power');
%4QAM 
hMod4 = comm.RectangularQAMModulator('ModulationOrder',4,'BitInput',true,'NormalizationMethod',...
   'Average power');
%BPSK
hMod2 = comm.BPSKModulator;

%% demodulation handles
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
% the number of bits in an OFDM signal.(������������ز�)
% generate the source bit
if FixSeedEnable == 1
  rng(Seed);   % set the seed ;
end
SNR_range = 10:2:30;
BER_Total_avg = zeros(length(SNR_range),4);
MSE_Total = zeros(length(SNR_range),4);
for snr_index = 1:length(SNR_range)
    SNR = SNR_range(snr_index);
    %save file
    filename = ['BER_SNR_',num2str(SNR),'dB.mat'];
    save(filename,'SNR');
    %channel estimation method
%     [cem, names] = enumeration('ChannelEstiMethod');
    cem = ChannelEstiMethod.IDEAL;
    for channelEstiMethod = cem' 
        switch(channelEstiMethod)
            case ChannelEstiMethod.IDEAL
                varname_BER_sc = ['BER_SNR_',num2str(SNR),'dB_IDEAL_sc'];
                varname_BER_avg = ['BER_SNR_',num2str(SNR),'dB_IDEAL_avg'];
                varname_MSE = ['MSE_SNR_',num2str(SNR),'dB_IDEAL'];
            case ChannelEstiMethod.LS
                varname_BER_sc = ['BER_SNR_',num2str(SNR),'dB_LS_sc'];
                varname_BER_avg = ['BER_SNR_',num2str(SNR),'dB_LS_avg'];
                varname_MSE = ['MSE_SNR_',num2str(SNR),'dB_LS'];
            case ChannelEstiMethod.LMMSE
                varname_BER_sc = ['BER_SNR_',num2str(SNR),'dB_LMMSE_sc'];
                varname_BER_avg = ['BER_SNR_',num2str(SNR),'dB_LMMSE_avg'];
                varname_MSE = ['MSE_SNR_',num2str(SNR),'dB_LMMSE'];
            case ChannelEstiMethod.SVD
                varname_BER_sc = ['BER_SNR_',num2str(SNR),'dB_SVD_sc'];
                varname_BER_avg = ['BER_SNR_',num2str(SNR),'dB_SVD_avg'];
                varname_MSE = ['MSE_SNR_',num2str(SNR),'dB_SVD'];
        end   
        nError_sum = zeros(nSC,1);
        MSE = zeros(nFrame,1);
        for nf = 1:nFrame
            source = randi([0,1],nBit_sig*nSignal,1);
            % save source in source_ini.
            source_ini = source;  
            source = reshape(source,nBit_sig,nSignal);
            dataMod = zeros(nSC,nSignal);
            % modulation bits by bits
            for ns = 1:nSignal
                bTH =1;  % b-th bit of a signal bit vector.
                for nc = nIdleLF+1:nSC-nIdleHF;  
                    switch modOrd(nc)
                      case 256
                           hMod = hMod256;
                      case 128
                           hMod = hMod128;
                      case 64
                           hMod = hMod64;
                      case 32
                           hMod = hMod32;  
                      case 16
                           hMod = hMod16;
                      case 8
                           hMod = hMod8;
                      case 4
                           hMod = hMod4;
                      case 2
                           hMod = hMod2;    
                      otherwise
                           continue;
                    end            
                    dataMod(nc, ns) = step(hMod,source(bTH:bTH+nModBit(nc)-1,ns));
                    bTH = bTH+nModBit(nc);         
                end
            end

            %% Power allocation
            dataMod = dataMod .* repmat(powerWet,1,nSignal);

            %% Hermitian Symmetry process
            dataHS = zeros(nSC*2,nSignal);
            dataHS(1,:) = real(dataMod(1,:));
            dataHS(2:nSC,:) = dataMod(2:nSC,:);
            dataHS(nSC+1,:) = imag(dataMod(1,:));
            dataHS(nSC+2:nSC*2,:) = flipud(conj(dataMod(2:nSC,:))); 

            %% IFFT
            dataTD = ifft(dataHS)*sqrt(double(2*nSC)); 

            %% Add CP
            dataTD_cp = zeros(nCP+2*nSC,nSignal);
            dataTD_cp(1:nCP,:) = dataTD(2*nSC-nCP+1:2*nSC,:);
            dataTD_cp(nCP+1:nCP+2*nSC,:) = dataTD;

            %% reshape to a column
            dataTD_col = reshape(dataTD_cp,(nCP+2*nSC)*nSignal,1);

            %% ZC pilot
            n = [1:nSC]';
            q = 0;
            r = 1;
            ZC = exp(-1i*2*pi*r*(n.*n/2+q*n)/nSC).*powerWet;%.* conj(repmat(powerWet,1,nSignal)');    % original ZC sequence of 1 cycle


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

            ZC_w_CP = zeros(1,nCP+2*nSC);
            ZC_w_CP(nCP+1:nCP+2*nSC) = ZC_time;
            ZC_w_CP(1:nCP) = ZC_time(2*nSC-nCP+1:2*nSC);

            %% Add ZC 
            dataTD_zc = zeros((nCP+2*nSC)*(nSignal+1),1);
            dataTD_zc(1:nCP+2*nSC) = ZC_w_CP;
            dataTD_zc(nCP+2*nSC+1:(nCP+2*nSC)*(nSignal+1)) = dataTD_cp;

            %% Pass Channel
            load Channel.mat;
            y_lp = conv(dataTD_zc,h_power_norm);
            noise = randn(length(y_lp),1)/sqrt(10^(SNR/10));
            y = y_lp+ noise;

            %% demodulation
            rx = y';
            L = length(rx);
            L_cor = L - (2*nSC) + 1;
            cor = zeros(1,L_cor);
            for i=1:L_cor
                cor(i) = rx(i:i+2*nSC-1)*ZC_time';
            end
            cor = cor/max(cor);
            threshold = 0.9;
            nShiftLeftStep = 4;
            peak_position_p = find(cor>threshold);
            peak_position = peak_position_p-nShiftLeftStep;

            ZC_rec = rx(peak_position:peak_position+2*nSC-1);
            dataTD_rec = rx(peak_position+2*nSC : peak_position+2*nSC+(2*nSC+nCP)*nSignal-1);

            %  estimate channel
            ZC_freq_rec = fft(ZC_rec) / sqrt(double(2*nSC));
            H = ZC_freq_rec./ ZC_freq;

            %% Channel Estimation Method Choice %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            load H_ideal.mat;
            switch(channelEstiMethod)
                case ChannelEstiMethod.IDEAL
                    H = H_ideal;
                case ChannelEstiMethod.LS
                    %H=H
                case ChannelEstiMethod.LMMSE
                    H = (M_LMMSE * H.').';
                case ChannelEstiMethod.SVD
                    H = (M_SVD * H.').';
            end    
            MSE(nf) = mean(abs(H_ideal-H).^2);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            nBit_sig = uint32(nBit_sig);
            % equlization and demodulation process frame by frame
            source_rec = zeros(nBit_sig*nSignal,1);
            nError = zeros(nSC,1);

            data = dataTD_rec;
            % reshape the data to remove cp
            data = reshape(data,2*nSC+nCP,nSignal);
            % remove cp
            data = data(nCP+1:2*nSC+nCP,:);
            % fft
            data_freq = fft(data) / sqrt(double(2*nSC));
            % equalization
            %% ֻ��1��ZC���й����ŵ�
            %     H = H_term(i,:);
            data_freq = conj(data_freq') ./ repmat(H,nSignal,1);
            data_freq(:,1) = data_freq(:,1)-data_freq(:,1+nSC)*1i;
            % remove the right part of frequency domain data
            data_freq = data_freq(:,1:nSC);
            % reverse it for conventient
            data_freq = conj(data_freq');
            % power allocation recovery
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
                     source_rec(bTH:bTH+nModBit(nc)-1) = step(hDemod,data_freq(nc,ns));
                      % count errors
                     nError(nc) = nError(nc) + sum(abs(source_rec(bTH:bTH+nModBit(nc)-1) - source_ini(bTH:bTH+nModBit(nc)-1)));
                     bTH = bTH+nModBit(nc);                 
                end
            end
            nError_sum = nError_sum + nError;
        end
        BER_sc= nError_sum/nFrame ./double(nModBit)/double(nSignal);
        BER_sc(1:nIdleLF)=1;
        BER_sc(nSC-nIdleHF+1:nSC)=1;
        BER_avg = sum(nError_sum)/sum(nModBit(nIdleLF+1:nSC-nIdleHF))/nFrame/nSignal;
        MSE_avg = mean(MSE);
        eval([varname_BER_sc,'=[',num2str(BER_sc'),']']); 
        eval([varname_BER_avg,'=[',num2str(BER_avg),']']); 
        eval([varname_MSE,'=[',num2str(MSE_avg),']']); 
        save(filename,varname_BER_sc,varname_BER_avg,varname_MSE,'-append');
        BER_Total_avg(snr_index,int32(channelEstiMethod))=BER_avg;
        MSE_Total(snr_index,int32(channelEstiMethod))=MSE_avg;
    end
end
save('BER_MSE_Total.mat','BER_Total_avg','MSE_Total');