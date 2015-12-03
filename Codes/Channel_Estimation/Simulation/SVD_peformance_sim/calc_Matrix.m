function [M_LMMSE,M_SVD] = calc_Matrix()

    %% loading parameters
    clear all;close all; clc;
    ReadInitialFile;
    nFrame = 10;
    nSignal = 1;
    H_cache = zeros(nFrame,2*nSC);
    for nf = 1:nFrame
        %% ZC pilot
        n = [1:nSC]';
        q = 0;
        r = 1;
        ZC = exp(-1i*2*pi*r*(n.*n/2+q*n)/nSC).*powerWet;%.* conj(repmat(powerWet,1,nSignal)');    % original ZC sequence of 1 cycle


        %% Generate Hermitian symmetric sequence
        ZC_freq = zeros(1,2*nSC);
        ZC_freq(1) = real(ZC(1));
        ZC_freq(2:nSC) = ZC(2:nSC);
        ZC_freq(nSC+1) = imag(ZC(1));
        ZC_freq(nSC+2:2*nSC) = conj(flipud(ZC(2:nSC)));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Correspondent time domain sequence
        ZC_time =ifft(ZC_freq)*sqrt(2*nSC);
          % ZC_time =ZC_time/max(abs(ZC_time));

        ZC_w_CP = zeros(1,nCP+2*nSC);
        ZC_w_CP(nCP+1:nCP+2*nSC) = ZC_time;
        ZC_w_CP(1:nCP) = ZC_time(2*nSC-nCP+1:2*nSC);

        %% Add ZC 
        dataTD_zc = zeros((nCP+2*nSC)*(nSignal+1),1);
        dataTD_zc(1:nCP+2*nSC) = ZC_w_CP;
       

        %% Pass Channel
        load Channel.mat;
        y_lp = conv(dataTD_zc,h_power_norm);
        SNR = 30;
        noise = randn(length(y_lp),1)/sqrt(10^(SNR/10));
        y = y_lp+ noise;

        %% demodulation
        rx = y';
        L = length(rx);
        L_cor = L - (2*nSC) + 1;
        cor = zeros(1,L_cor);
        for i=1:L_cor
            cor(i) = rx(i:i+2*nSC-1)*ZC_time';
        end
        cor = cor/max(cor);
        threshold = 0.9;
        nShiftLeftStep = 4;
        peak_position_p = find(cor>threshold);
        peak_position = peak_position_p-nShiftLeftStep;

        ZC_rec = rx(peak_position:peak_position+2*nSC-1);
        dataTD_rec = rx(peak_position+2*nSC : peak_position+2*nSC+(2*nSC+nCP)*nSignal-1);

        %  estimate channel
        ZC_freq_rec = fft(ZC_rec) / sqrt(double(2*nSC));
        H = ZC_freq_rec./ ZC_freq;
        H_cache(nf,:) = H;
    end
    R_HH = zeros(2*nSC);
    for i = 1:nFrame
        H = H_cache(i,:).';
        R_HH = R_HH + H*H';
    end
    R_HH = R_HH/nFrame;
    SNR = 25;
    snr = 1/10/log10(SNR);
    M_LMMSE = R_HH * inv(R_HH+snr*diag(ones(1,2*nSC)));
    [v,d] = eigs(R_HH,256);
    for k=1:256
        if k<=256
           d(k,k) = d(k,k)/(d(k,k)+snr);
        else
            d(k,k) =0;
        end
    end
    M_SVD = v*d*v';
end