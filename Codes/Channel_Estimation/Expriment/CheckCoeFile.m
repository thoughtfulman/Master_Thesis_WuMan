%% This program is for reading data in .coe files
% the purpose of reading data in .coe files is to check
% the correction of demodulation program.
% author: Man Wu
% build time: 2014-8-1
% revised by: Man Wu
clear all;
close all;
clc;
% addpath('E:\WuMan\Program\Matlab\SourceCode');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file process program just copy from modulation part
% open the file.
LMMSE_Enable = 0;
nShiftLeftStep = 12;

ReadInitialFile;
fid_h = fopen(Output_file_name_h);

inData1_coe = textscan(fid_h, '%*[^\n]',2,'commentStyle','\\');   
if pilotEnable ==1
    inData2_coe =  textscan(fid_h, '%d,',(2*nSC+nCP)*(nSignal+2),'commentStyle','\\');
else
    inData2_coe =  textscan(fid_h, '%d,',(2*nSC+nCP)*(nSignal+1),'commentStyle','\\');
end
    
fclose(fid_h);

if pilotEnable ==1
    Data_coe = zeros((2*nSC+nCP)*(nSignal+2)+(2*nSC+nCP),1);
    for i=1:(2*nSC+nCP)*(nSignal+2)
      Data_coe(i) = inData2_coe{1}(i);
    end
    Data_coe((2*nSC+nCP)*(nSignal+2)+1:(2*nSC+nCP)*(nSignal+2)+(2*nSC+nCP)) = Data_coe(1:2*nSC+nCP);
else
    Data_coe = zeros((2*nSC+nCP)*(nSignal+1)+(2*nSC+nCP),1);
    for i=1:(2*nSC+nCP)*(nSignal+1)
      Data_coe(i) = inData2_coe{1}(i);
    end
   Data_coe((2*nSC+nCP)*(nSignal+1)+1:(2*nSC+nCP)*(nSignal+1)+(2*nSC+nCP)) = Data_coe(1:2*nSC+nCP);
end
ADC4_D = Data_coe;
%% 加噪声，试试EVM估计准不准
% SNR = 30; %dB
% ADC4_D(281:length(ADC4_D)) = ADC4_D(281:length(ADC4_D)) + sqrt( 10^(-SNR/10))*2^13 * randn(length(ADC4_D)-280,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
resFileName = 'checkResult.xlsx';
DCOOFDM_demodulation;