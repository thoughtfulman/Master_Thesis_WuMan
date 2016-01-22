clear all;
close all;
clc;

addpath('E:\WuMan\Program\Matlab\DCO_OFDM_Sim_v1');
xlLoadChipScopeData('H5_200M_80cm_da60db_ad710db_agr65_1.prn');
resFileName = 'H5_200M_80cm_da60db_ad710db_agr65_1.xlsx';
LMMSE_Enable = 0;
nShiftLeftStep = 12;   % in channel estimation

ReadInitialFile;
DCOOFDM_demodulation;

if pilotEnable == 1
figure(1);
plot(abs(H(1:256)));
hold on;
plot(abs(H_term(1:256)),'r')
legend('Pilot','ZC');
figure(2);
plot(phase(H(1:256)));
hold on;
plot(phase(H_term(1:256)),'r')
legend('Pilot','ZC');
end







