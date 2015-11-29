clear all;
close all;
clc;

SNR = 20;
%BER
BER = zeros(length(SNR),1);
% M-QAM
M = 64; 
% 星座点归一化参数
NormNumber = sqrt(3/2/(M-1));
% 子载波数
SC = 128;
%调制器
hMod = comm.RectangularQAMModulator(M); 
%解调器
hDemod = comm.RectangularQAMDemodulator(M); 

for i=1:length(SNR)
    snr = SNR(i);
    %每个SNR下仿真的帧数
    NF = 2000;
    sumerrbit = 0;
    sumbit = log2(M)*SC*NF;
    for j = 1:NF
        % 每帧的符号数
        NS = 20;
        % 信源
        source = randi(M,SC,1)-1;
        % 归一化星座点
        X = NormNumber * step(hMod,source);
        % 噪声方差
        % Zsigma = sqrt(1/log2(M)/2/(10^(snr/10)));
        Zsigma = sqrt(1/2/(10^(snr/10)));
        % 频域噪声
        Z = Zsigma*(randn(SC,1)+randn(SC,1)*1i);
        
        % 过信道
        Y = X+ Z;
        
        % 解调
        y = step(hDemod,Y/NormNumber);
        % 统计错误比特数
        [nerr,rerr] = biterr(source,y);
        sumerrbit = sumerrbit + nerr;
    end
    BER(i) = sumerrbit/sumbit;
end
semilogy(SNR,BER);
grid on;