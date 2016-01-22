 clear all;
 close all;
 clc;
 
data1=xlsread('E:\WuMan\TestData\20150421\BPSK_20sig\H5_200M_80cm_da60db_ad710db_agr65_1.xlsx');
data2=xlsread('E:\WuMan\TestData\20150421\BPSK_20sig\H5_200M_80cm_da60db_ad710db_agr65_2.xlsx');

data3=xlsread('E:\WuMan\TestData\20150421\BPSK_20sig_1\H5_200M_80cm_da60db_ad710db_agr65_1.xlsx');
data4=xlsread('E:\WuMan\TestData\20150421\BPSK_20sig_1\H5_200M_80cm_da60db_ad710db_agr65_2.xlsx');

data5=xlsread('E:\WuMan\TestData\20150421\BPSK_20sig_2\H5_200M_80cm_da60db_ad710db_agr65_1.xlsx');
data6=xlsread('E:\WuMan\TestData\20150421\BPSK_20sig_2\H5_200M_80cm_da60db_ad710db_agr65_2.xlsx');

data7=xlsread('E:\WuMan\TestData\20150421\BPSK_20sig_3\H5_200M_80cm_da60db_ad710db_agr65_1.xlsx');
data8=xlsread('E:\WuMan\TestData\20150421\BPSK_20sig_3\H5_200M_80cm_da60db_ad710db_agr65_2.xlsx');

data9=xlsread('E:\WuMan\TestData\20150421\BPSK_20sig_4\H5_200M_80cm_da60db_ad710db_agr65_1.xlsx');
data10=xlsread('E:\WuMan\TestData\20150421\BPSK_20sig_4\H5_200M_80cm_da60db_ad710db_agr65_2.xlsx');

SNR_evm_bpsk = data1(:,8)+data2(:,8)+data3(:,8)+data4(:,8)+data5(:,8)+data6(:,8)+data7(:,8)+data8(:,8)+data9(:,8)+data10(:,8);
SNR_evm_bpsk = SNR_evm_bpsk/10;
plot((5:124),SNR_evm_bpsk(5:124),'-s');
hold on;
SNR_zc_bpsk = data1(:,7)+data2(:,7)+data3(:,7)+data4(:,7)+data5(:,7)+data6(:,7)+data7(:,7)+data8(:,7)+data9(:,7)+data10(:,7);
SNR_zc_bpsk = SNR_zc_bpsk/10;
plot((5:124),SNR_zc_bpsk(5:124),'-+k');

data11=xlsread('E:\WuMan\TestData\20150421\QAM4_20sig_0\H5_200M_80cm_da60db_ad710db_agr65_1.xlsx');
data12=xlsread('E:\WuMan\TestData\20150421\QAM4_20sig_0\H5_200M_80cm_da60db_ad710db_agr65_2.xlsx');

data13=xlsread('E:\WuMan\TestData\20150421\QAM4_20sig_1\H5_200M_80cm_da60db_ad710db_agr65_1.xlsx');
data14=xlsread('E:\WuMan\TestData\20150421\QAM4_20sig_1\H5_200M_80cm_da60db_ad710db_agr65_2.xlsx');

data15=xlsread('E:\WuMan\TestData\20150421\QAM4_20sig_2\H5_200M_80cm_da60db_ad710db_agr65_1.xlsx');
data16=xlsread('E:\WuMan\TestData\20150421\QAM4_20sig_2\H5_200M_80cm_da60db_ad710db_agr65_2.xlsx');

data17=xlsread('E:\WuMan\TestData\20150421\QAM4_20sig_3\H5_200M_80cm_da60db_ad710db_agr65_1.xlsx');
data18=xlsread('E:\WuMan\TestData\20150421\QAM4_20sig_3\H5_200M_80cm_da60db_ad710db_agr65_2.xlsx');

data19=xlsread('E:\WuMan\TestData\20150421\QAM4_20sig_4\H5_200M_80cm_da60db_ad710db_agr65_1.xlsx');
data20=xlsread('E:\WuMan\TestData\20150421\QAM4_20sig_4\H5_200M_80cm_da60db_ad710db_agr65_2.xlsx');

% SNR_evm_4QAM = data11(:,8)+data12(:,8)+data13(:,8)+data14(:,8)+data15(:,8)+data16(:,8)+data17(:,8)+data18(:,8)+data19(:,8)+data20(:,8);
% SNR_evm_4QAM = SNR_evm_4QAM/10;
% plot((5:124)*100/128,SNR_evm_4QAM(5:124),'-*r');
% hold on;
% SNR_zc_4QAM = data11(:,7)+data12(:,7)+data13(:,7)+data14(:,7)+data15(:,7)+data16(:,7)+data17(:,7)+data18(:,7)+data19(:,7)+data20(:,7);
% SNR_zc_4QAM = SNR_zc_4QAM/10;
% plot((5:124)*100/128,SNR_zc_4QAM(5:124),'-*g');

grid on;
xlabel('Frequency (MHz)');
ylabel('SNR (dB)');
legend('EVM','ZC');