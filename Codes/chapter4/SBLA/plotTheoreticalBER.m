clear all;close all; clc;

load theoreticalBER.mat;
figure(1);
linewidth = 1;
SNR = 0:2:40;
a = berAll(:,1);
semilogy(SNR,berAll(:,1),'-*b','LineWidth',linewidth);
hold on;
semilogy(SNR,berAll(:,2),'-+r','LineWidth',linewidth);
semilogy(SNR,berAll(:,3),'-<k','LineWidth',linewidth);
semilogy(SNR,berAll(:,4),'->g','LineWidth',linewidth);
semilogy(SNR,berAll(:,5),'-^c','LineWidth',linewidth);
semilogy(SNR,berAll(:,6),'-vm','LineWidth',linewidth);
semilogy(SNR,berAll(:,7),'-sb','LineWidth',linewidth);
semilogy(SNR,berAll(:,8),'-or','LineWidth',linewidth);
semilogy(SNR,berAll(:,9),'-xk','LineWidth',linewidth);
semilogy(SNR,berAll(:,10),'-dc','LineWidth',linewidth);
xlabel('SNR (dB)');
ylabel('BER');
grid on;
axis([0,44,10^-5,1]);
legend('BPSK','4QAM','8QAM','16QAM','32QAM','64QAM',...
    '128QAM','256QAM','512QAM','1024QAM');