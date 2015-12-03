close all;clear all;clc;

load BER_MSE_Total_nF10k1.mat;
BER_Total_avg_nf10k1 = BER_Total_avg;
MSE_Total_nf10k1 = MSE_Total;
load BER_MSE_Total_nF10k2.mat;
BER_Total_avg_nf10k2 = BER_Total_avg;
MSE_Total_nf10k2 = MSE_Total;
load BER_MSE_Total_nF10k5.mat;
BER_Total_avg_nf10k5 = BER_Total_avg;
MSE_Total_nf10k5 = MSE_Total;
load BER_MSE_Total_nF10k256.mat;
BER_Total_avg_nf10k256 = BER_Total_avg;
MSE_Total_nf10k256 = MSE_Total;
SNR = 15:5:35;
linewidth = 1.3;
% ideal
semilogy(SNR,BER_Total_avg_nf10k2(2:6,1),'->','LineWidth',linewidth);
hold on;
% LMMSE
% semilogy(SNR,BER_Total_avg_nf10k2(2:6,3),'k');
semilogy(SNR,BER_Total_avg_nf10k1(:,4),'-<r','LineWidth',linewidth);
semilogy(SNR,BER_Total_avg_nf10k5(:,4),'-sk','LineWidth',linewidth);
semilogy(SNR,BER_Total_avg_nf10k256(:,4),'-+m','LineWidth',linewidth);
grid on;
xlabel('SNR (dB)');
ylabel('BER')
legend('Ideal','SVD, M=1','SVD, M=5','LMMSE');


figure(2);
% semilogy(SNR,MSE_Total_nf10k2(2:6,1),'->','LineWidth',linewidth);
% semilogy(SNR,MSE_Total_nf10k2(2:6,3),'k');
semilogy(SNR,MSE_Total_nf10k1(:,4),'-<r','LineWidth',linewidth);
hold on;
semilogy(SNR,MSE_Total_nf10k5(:,4),'-sk','LineWidth',linewidth);
semilogy(SNR,MSE_Total_nf10k256(:,4),'-+m','LineWidth',linewidth);
legend('SVD, M=1','SVD, M=5','LMMSE');
xlabel('SNR (dB)');
ylabel('MSE')
grid on;