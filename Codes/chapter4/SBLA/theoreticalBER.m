clear all;
close all;
clc;

nSignal = 5000; % number of signals to simulation
nSC = 128;
SNR = 0:2:40;

berAll = zeros(length(SNR),10);
for b=1:10
    M = 2^b;
    hMod = comm.RectangularQAMModulator('ModulationOrder',M,'BitInput',false,'NormalizationMethod',...
       'Average power');
    hDemod = comm.RectangularQAMDemodulator('ModulationOrder',M,'BitOutput',false,'NormalizationMethod',...
         'Average power');
    source = randi(M-1,nSignal*nSC,1);
    data_Mod = step(hMod,source);
    data_Mod = reshape(data_Mod,nSC,nSignal);
    dataHS = zeros(nSC*2,nSignal);
    dataHS(1,:) = real(data_Mod(1,:));
    dataHS(2:nSC,:) = data_Mod(2:nSC,:);
    dataHS(nSC+1,:) = imag(data_Mod(1,:));
    dataHS(nSC+2:nSC*2,:) = flipud(conj(data_Mod(2:nSC,:))); 
    % ifft
    dataTD = ifft(dataHS)*sqrt(double(2*nSC)); 

    BER = zeros(length(SNR),1);
    for i = 1:length(SNR)
    n = sqrt(1  / 10^(SNR(i)/10) ) * randn(2*nSC,nSignal) ;
    dataTD_rec = dataTD + n;
    dataHS_rec = fft(dataTD_rec) / sqrt(double(2*nSC));
    dataHS_rec(1,:) = dataHS_rec(1,:) + 1i * dataHS_rec(1+nSC,:);
    dataHS_rec = dataHS_rec(1:nSC,:);
    source_rec = step(hDemod,reshape(dataHS_rec,nSignal*nSC,1));
    [Nerr,BER(i)] = biterr(source,source_rec);
    end
    berAll(:,b) = BER;
end
save('berTheoretical.mat','berAll');
figure(1);
% semilogy(SNR,BER_256);
% 
% % 64 QAM
% M = 64;
% hMod = comm.RectangularQAMModulator('ModulationOrder',M,'BitInput',false,'NormalizationMethod',...
%    'Average power');
% hDemod = comm.RectangularQAMDemodulator('ModulationOrder',M,'BitOutput',false,'NormalizationMethod',...
%      'Average power');
% source = randi(M-1,nSignal*nSC,1);
% data_Mod = step(hMod,source);
% data_Mod = reshape(data_Mod,nSC,nSignal);
% dataHS = zeros(nSC*2,nSignal);
% dataHS(1,:) = real(data_Mod(1,:));
% dataHS(2:nSC,:) = data_Mod(2:nSC,:);
% dataHS(nSC+1,:) = imag(data_Mod(1,:));
% dataHS(nSC+2:nSC*2,:) = flipud(conj(data_Mod(2:nSC,:))); 
% % ifft
% dataTD = ifft(dataHS)*sqrt(double(2*nSC)); 
% 
% BER = zeros(length(SNR),1);
% for i = 1:length(SNR)
% n = sqrt(1  / 10^(SNR(i)/10) ) * randn(2*nSC,nSignal) ;
% dataTD_rec = dataTD + n;
% dataHS_rec = fft(dataTD_rec) / sqrt(double(2*nSC));
% dataHS_rec(1,:) = dataHS_rec(1,:) + 1i * dataHS_rec(1+nSC,:);
% dataHS_rec = dataHS_rec(1:nSC,:);
% source_rec = step(hDemod,reshape(dataHS_rec,nSignal*nSC,1));
% [Nerr,BER(i)] = biterr(source,source_rec);
% end
% BER_64 = BER;
% 
% % 16 QAM
% M = 16;
% hMod = comm.RectangularQAMModulator('ModulationOrder',M,'BitInput',false,'NormalizationMethod',...
%    'Average power');
% hDemod = comm.RectangularQAMDemodulator('ModulationOrder',M,'BitOutput',false,'NormalizationMethod',...
%      'Average power');
% source = randi(M-1,nSignal*nSC,1);
% data_Mod = step(hMod,source);
% data_Mod = reshape(data_Mod,nSC,nSignal);
% dataHS = zeros(nSC*2,nSignal);
% dataHS(1,:) = real(data_Mod(1,:));
% dataHS(2:nSC,:) = data_Mod(2:nSC,:);
% dataHS(nSC+1,:) = imag(data_Mod(1,:));
% dataHS(nSC+2:nSC*2,:) = flipud(conj(data_Mod(2:nSC,:))); 
% % ifft
% dataTD = ifft(dataHS)*sqrt(double(2*nSC)); 
% 
% BER = zeros(length(SNR),1);
% for i = 1:length(SNR)
% n = sqrt(1  / 10^(SNR(i)/10) ) * randn(2*nSC,nSignal) ;
% dataTD_rec = dataTD + n;
% dataHS_rec = fft(dataTD_rec) / sqrt(double(2*nSC));
% dataHS_rec(1,:) = dataHS_rec(1,:) + 1i * dataHS_rec(1+nSC,:);
% dataHS_rec = dataHS_rec(1:nSC,:);
% source_rec = step(hDemod,reshape(dataHS_rec,nSignal*nSC,1));
% [Nerr,BER(i)] = biterr(source,source_rec);
% end
% BER_16 = BER;
% 
% % 4 QAM
% M = 4;
% hMod = comm.RectangularQAMModulator('ModulationOrder',M,'BitInput',false,'NormalizationMethod',...
%    'Average power');
% hDemod = comm.RectangularQAMDemodulator('ModulationOrder',M,'BitOutput',false,'NormalizationMethod',...
%      'Average power');
% source = randi(M-1,nSignal*nSC,1);
% data_Mod = step(hMod,source);
% data_Mod = reshape(data_Mod,nSC,nSignal);
% dataHS = zeros(nSC*2,nSignal);
% dataHS(1,:) = real(data_Mod(1,:));
% dataHS(2:nSC,:) = data_Mod(2:nSC,:);
% dataHS(nSC+1,:) = imag(data_Mod(1,:));
% dataHS(nSC+2:nSC*2,:) = flipud(conj(data_Mod(2:nSC,:))); 
% % ifft
% dataTD = ifft(dataHS)*sqrt(double(2*nSC)); 
% 
% BER = zeros(length(SNR),1);
% for i = 1:length(SNR)
% n = sqrt(1 / 10^(SNR(i)/10) ) * randn(2*nSC,nSignal) ;
% dataTD_rec = dataTD + n;
% dataHS_rec = fft(dataTD_rec) / sqrt(double(2*nSC));
% dataHS_rec(1,:) = dataHS_rec(1,:) + 1i * dataHS_rec(1+nSC,:);
% dataHS_rec = dataHS_rec(1:nSC,:);
% source_rec = step(hDemod,reshape(dataHS_rec,nSignal*nSC,1));
% [Nerr,BER(i)] = biterr(source,source_rec);
% end
% BER_4 = BER;
% 
% % plot the result
% linewidth = 2;
% semilogy(SNR,BER_4,'b','LineWidth',linewidth);
% hold on;
% 
% semilogy(SNR,BER_16,'r','LineWidth',linewidth);
% 
% semilogy(SNR,BER_64,'k','LineWidth',linewidth);
% 
% semilogy(SNR,BER_256,'g','LineWidth',linewidth);
% grid on;
% xlabel('SNR (dB)');
% ylabel('BER');
% axis([0,35,10^-5,1]);
% legend('4 QAM','16 QAM','64 QAM','256 QAM');
