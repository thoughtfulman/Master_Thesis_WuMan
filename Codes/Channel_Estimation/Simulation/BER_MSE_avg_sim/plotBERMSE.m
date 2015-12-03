clear all;close all;clc;
load BER_MSE_Total_svd_nF10k2.mat;
% load BER_MSE_Total.mat;
SNR = 15:5:35;
linewidth = 1.3;
semilogy(SNR,BER_Total_avg(2:6,1),'-<','LineWidth',linewidth);
hold on;
semilogy(SNR,BER_Total_avg(2:6,2),'->m','LineWidth',linewidth);
semilogy(SNR,BER_Total_avg(2:6,3),'--sk','LineWidth',linewidth);
semilogy(SNR,BER_Total_avg(2:6,4),'-+r','LineWidth',linewidth);
grid on;
xlabel('SNR (dB)');
ylabel('BER');
legend('Ideal','LS','LMMSE','SVD');

figure(2);

semilogy(SNR,MSE_Total(2:6,2),'->m','LineWidth',linewidth);
hold on;
semilogy(SNR,MSE_Total(2:6,3),'--sk','LineWidth',linewidth);
semilogy(SNR,MSE_Total(2:6,4),'-+r','LineWidth',linewidth);
grid on;
xlabel('SNR (dB)');
ylabel('MSE');
legend('LS','LMMSE','SVD');