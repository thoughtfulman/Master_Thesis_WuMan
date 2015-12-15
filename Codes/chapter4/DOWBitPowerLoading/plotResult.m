clear all;close all;clc;


load BER.mat;
SNR = 10:2:30;
figure(3);
semilogy(SNR,BER_Total_avg(:,4),'--m');
hold on;
semilogy(SNR,BER_Total_avg(:,1),'-*');
semilogy(SNR,BER_Total_avg(:,2),'-sk');
semilogy(SNR,BER_Total_avg(:,3),'-r+');

grid on;
xlabel('SNR (dB)');
ylabel('BER');
legend('Allocation in Table 3.2','Hughers-Hartogs','P.S.Chow','Fischer');
