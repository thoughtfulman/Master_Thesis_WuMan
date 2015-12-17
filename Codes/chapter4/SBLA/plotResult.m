clear all;close all;clc;


load BER_SBLA.mat;
SNR = 10:2:30;
figure(1);
semilogy(SNR,BER_Total_avg(:,6),'-om');
hold on;
semilogy(SNR,BER_Total_avg(:,5),'-*k');
semilogy(SNR,BER_Total_avg(:,4),'-+');
semilogy(SNR,BER_Total_avg(:,7),'-sk');
semilogy(SNR,BER_Total_avg(:,2),'-^r');

grid on;
xlabel('SNR (dB)');
ylabel('BER');
legend('Fixed Modulation','Allocation in Table 3.2','SBLA','Improved-SBLA','Chow');

load Result_Chow_SNR_25dB.mat;
load Result_ImprovedSBLA_SNR_25dB.mat;
load Result_SBLA_SNR_25dB.mat;

figure(2);
plot(loadedBit_IDEAL_Chow_SNR_25dB,'-+r');
hold on;
plot(loadedBit_IDEAL_SBLA_SNR_25dB,'-sb');
% plot(loadedBit_IDEAL_ImprovedSBLA_SNR_25dB,'-k');
grid on;
xlabel('Subcarrier Index');
ylabel('Bit Allocation');
legend('Chow','SBLA, Improved-SBLA');
figure(3);
plot(loadedPower_IDEAL_Chow_SNR_25dB,'-+r');
hold on;
plot(loadedPower_IDEAL_SBLA_SNR_25dB,'-sb');
plot(loadedPower_IDEAL_ImprovedSBLA_SNR_25dB,'--k');
grid on;
xlabel('Subcarrier Index');
ylabel('Bit Allocation');
legend('Chow','SBLA','Improved-SBLA');