clear all;close all; clc;
% load Channel.mat;
ReadInitialFile;
% H = fft(h_power_norm,256);
% SNR = 10*log10(abs(H).^2/Sigma);
% plot(1:nSC,SNR(1:nSC));
% modulation handles
%1024QAM
hMod2048 = comm.RectangularQAMModulator('ModulationOrder',2048,'BitInput',true,'NormalizationMethod',...
   'Average power');
%1024QAM
hMod1024 = comm.RectangularQAMModulator('ModulationOrder',1024,'BitInput',true,'NormalizationMethod',...
   'Average power');
%512QAM
hMod512 = comm.RectangularQAMModulator('ModulationOrder',512,'BitInput',true,'NormalizationMethod',...
   'Average power');
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
%2048QAM
hDemod2048 = comm.RectangularQAMDemodulator('ModulationOrder',2048,'BitOutput',true,'NormalizationMethod',...
     'Average power');
%1024QAM
hDemod1024 = comm.RectangularQAMDemodulator('ModulationOrder',1024,'BitOutput',true,'NormalizationMethod',...
     'Average power');
%512QAM
hDemod512 = comm.RectangularQAMDemodulator('ModulationOrder',512,'BitOutput',true,'NormalizationMethod',...
     'Average power');
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
    totalPower = nSC;
    targetRate = (nSC-nIdleLF-nIdleHF)*6;
    ALGO = {'Hughers_Hartogs';'Chow';'Fischer'};
    numChannelTest = 1;
    nErrorChannelTest = zeros(nSC,3);
    for nCT = 1:numChannelTest
        H = genChannel();
        for algorithmName = ALGO'
            if(strcmp(algorithmName,'Hughers_Hartogs'))
               [loadedBit,loadedPower]=...
                BitPowerLoadingAlgorithm.Hughers_Hartogs(...
                    H,SNR,totalPower,targetRate);
            elseif(strcmp(algorithmName,'Chow'))
               [loadedBit,loadedPower]=...
                BitPowerLoadingAlgorithm.Chow(...
                    H,SNR,totalPower,targetRate);
            elseif(strcmp(algorithmName,'Fischer'))
               [loadedBit,loadedPower]=...
                BitPowerLoadingAlgorithm.Fischer(...
                    H,SNR,totalPower,targetRate);
            else
                error('no such algorithmName');

            end
            nModBit = loadedBit;
            powerWet = loadedPower;
            modOrd = 2.^nModBit;
            nBit_sig = sum(log2(modOrd(nIdleLF+1:nSC-nIdleHF)));
            %save file
%             !mkdir ./resTarget_4;
%             filename = ['./resTarget_4/Result_',char(algorithmName),'_SNR_',num2str(SNR),'dB.mat'];
%             save(filename,'SNR');
            %save(filename,'loadedBit','loadedPower','-append');
            %channel estimation method
            cem = ChannelEstiMethod.IDEAL;
            for channelEstiMethod = cem' 
                switch(channelEstiMethod)
                    case ChannelEstiMethod.IDEAL
                        varname_BER_sc = ['berOnCarriers_IDEAL_',char(algorithmName),'_SNR_',num2str(SNR),'dB'];
                        varname_BER_avg =['berAverage_IDEAL_',char(algorithmName),'_SNR_',num2str(SNR),'dB'];
                        varname_MSE = ['MSE_IDEAL_',char(algorithmName),'_SNR_',num2str(SNR),'dB'];
                        varname_loadedBit = ['loadedBit_IDEAL_',char(algorithmName),'_SNR_',num2str(SNR),'dB'];
                        varname_loadedPower = ['loadedPower_IDEAL_',char(algorithmName),'_SNR_',num2str(SNR),'dB'];
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
                              case 2048
                                   hMod = hMod2048;
                              case 1024
                                   hMod = hMod1024;
                              case 512
                                   hMod = hMod512;
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


                    %% Pass Channel
                    Zsigma = sqrt(1/2/(10^(SNR/10)));  
                    Z = Zsigma*(randn(nSC,nSignal)+randn(nSC,nSignal)*1i);
                    data_freq = dataMod .* repmat(H(1:nSC).',1,nSignal)+ Z;
                    %% demodulation
                    data_freq = data_freq ./repmat(H(1:nSC).',1,nSignal) ./ repmat(powerWet,1,nSignal);
                    nBit_sig = uint32(nBit_sig);
                    % equlization and demodulation process frame by frame
                    source_rec = zeros(nBit_sig*nSignal,1);
                    nError = zeros(nSC,1);
                    bTH = 1;  % use to label the bit order in demodulated data
                    for ns = 1:nSignal        
                        for nc = nIdleLF+1:nSC-nIdleHF;
                            % choose demodulation handle
                              switch modOrd(nc)
                                  case 2048
                                      hDemod = hDemod2048;
                                  case 1024
                                      hDemod = hDemod1024;
                                  case 512
                                      hDemod = hDemod512;
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
            end
            if(strcmp(algorithmName,'Hughers_Hartogs'))
                nErrorChannelTest(:,1) = nErrorChannelTest(:,1)+nError_sum;
            elseif(strcmp(algorithmName,'Chow'))
                nErrorChannelTest(:,2) = nErrorChannelTest(:,2)+nError_sum;
            elseif(strcmp(algorithmName,'Fischer'))
                nErrorChannelTest(:,3) = nErrorChannelTest(:,3)+nError_sum;
            else
                error('no such algorithmName');

            end
        end
%             BER_sc= nErrorChannelTest/nFrame/numChannelTest ./double(nModBit)/double(nSignal);
%             BER_sc(1:nIdleLF)=1;
%             BER_sc(nSC-nIdleHF+1:nSC)=1;
%             BER_avg = sum(nErrorChannelTest)/sum(nModBit(nIdleLF+1:nSC-nIdleHF))/numChannelTest/nFrame/nSignal;
%             MSE_avg = mean(MSE);
%             eval([varname_BER_sc,'=[',num2str(BER_sc'),']']); 
%             eval([varname_BER_avg,'=[',num2str(BER_avg),']']); 
%             eval([varname_MSE,'=[',num2str(MSE_avg),']']); 
%             eval([varname_loadedBit,'=[',num2str(loadedBit'),']']);
%             eval([varname_loadedPower,'=[',num2str(loadedPower'),']']);
%             save(filename,varname_BER_sc,varname_BER_avg,...
%             varname_MSE,varname_loadedBit,varname_loadedPower,'-append');
%             if(strcmp(algorithmName,'Hughers_Hartogs'))
%                 BER_Total_avg(snr_index,1)=BER_avg;
%             elseif(strcmp(algorithmName,'Chow'))
%                 BER_Total_avg(snr_index,2)=BER_avg;
%             elseif(strcmp(algorithmName,'Fischer'))
%                 BER_Total_avg(snr_index,3)=BER_avg;
%             else
%                 error('no such algorithmName');
% 
%             end
    end
    BER_Total_avg(snr_index,1:3) = sum(nErrorChannelTest)/sum(nModBit(nIdleLF+1:nSC-nIdleHF))/numChannelTest/nFrame/nSignal;
end
