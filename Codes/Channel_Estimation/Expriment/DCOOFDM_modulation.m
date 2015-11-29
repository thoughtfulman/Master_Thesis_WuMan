%% program statement %%%%%%%%%%%%%%%%%%%%%%%
% This program is for DCO-OFDM modulation, supplies various number of 
% subcarriers(nSC), number of subcarriers in cluster(nSC_clt), number of CP
%(CP), modulation order on each cluster(modOrd) and cluster power weight.
% allocation(powerWet).
% all arguments are defined in the file 'modulation_initial_file.txt'

% author: Man Wu
% build time: 2014-7-14
% revised by: Man Wu
% revised time: 2014-7-16
% revised statement: add modulation code.
% revised time: 2014-7-30
% revised statement: add ZC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%% modulation 

% modulation handles
%256QAM
hMod256 = comm.RectangularQAMModulator('ModulationOrder',256,'BitInput',true,'NormalizationMethod',...
   'Average power');
%64QAM
hMod128 = comm.RectangularQAMModulator('ModulationOrder',128,'BitInput',true,'NormalizationMethod',...
   'Average power');
%64QAM
hMod64 = comm.RectangularQAMModulator('ModulationOrder',64,'BitInput',true,'NormalizationMethod',...
   'Average power');
%32QAM
hMod32 = comm.RectangularQAMModulator('ModulationOrder',32,'BitInput',true,'NormalizationMethod',...
   'Average power');
%16QAM 
hMod16 = comm.RectangularQAMModulator('ModulationOrder',16,'BitInput',true,'NormalizationMethod',...
   'Average power');
%4QAM
hMod8 = comm.RectangularQAMModulator('ModulationOrder',8,'BitInput',true,'NormalizationMethod',...
   'Average power');
%4QAM 
hMod4 = comm.RectangularQAMModulator('ModulationOrder',4,'BitInput',true,'NormalizationMethod',...
   'Average power');
%BPSK
hMod2 = comm.BPSKModulator;

% the number of bits in an OFDM signal.(除掉了虚拟子载波)
% generate the source bit
if FixSeedEnable == 1
  rng(Seed);   % set the seed ;
end
source = randi([0,1],nBit_sig*nSignal,1);
  
% save source in source_ini.
source_ini = source;  
source = reshape(source,nBit_sig,nSignal);
dataMod = zeros(nSC,nSignal);
% modulation bits by bits
for ns = 1:nSignal
    bTH =1;  % b-th bit of a signal bit vector.
    for nc = nIdleLF+1:nSC-nIdleHF;  
        switch modOrd(nc)
          case 256
               hMod = hMod256;
          case 128
               hMod = hMod128;
          case 64
               hMod = hMod64;
          case 32
               hMod = hMod32;  
          case 16
               hMod = hMod16;
          case 8
               hMod = hMod8;
          case 4
               hMod = hMod4;
          case 2
               hMod = hMod2;    
          otherwise
               continue;
        end            
        dataMod(nc, ns) = step(hMod,source(bTH:bTH+nModBit(nc)-1,ns));
        bTH = bTH+nModBit(nc);         
    end
end

%% Power allocation
dataMod = dataMod .* repmat(powerWet,1,nSignal);
  
%% Hermitian Symmetry process
dataHS = zeros(nSC*2,nSignal);
dataHS(1,:) = real(dataMod(1,:));
dataHS(2:nSC,:) = dataMod(2:nSC,:);
dataHS(nSC+1,:) = imag(dataMod(1,:));
dataHS(nSC+2:nSC*2,:) = flipud(conj(dataMod(2:nSC,:))); 

%% IFFT
dataTD = ifft(dataHS)*sqrt(double(2*nSC)); 
  
%% Add CP
dataTD_cp = zeros(nCP+2*nSC,nSignal);
dataTD_cp(1:nCP,:) = dataTD(2*nSC-nCP+1:2*nSC,:);
dataTD_cp(nCP+1:nCP+2*nSC,:) = dataTD;
  
%% reshape to a column
dataTD_col = reshape(dataTD_cp,(nCP+2*nSC)*nSignal,1);
  
%% ZC pilot
n = [1:nSC]';
q = 0;
r = 1;
ZC = exp(-1i*2*pi*r*(n.*n/2+q*n)/nSC).*powerWet;%.* conj(repmat(powerWet,1,nSignal)');    % original ZC sequence of 1 cycle

%% configure virtual carriers
ZC(1:nIdleLF) =0;
ZC(nSC-nIdleHF+1:nSC) =0;
  
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
dataTD_zc(nCP+2*nSC+1:(nCP+2*nSC)*(nSignal+1)) = dataTD_cp;

%% DA format
% 归一化
dataTD_zc = dataTD_zc * 2^13;
depth=length(dataTD_zc);                 %存储单元数 对应数据长度;
width=16;                    %数据宽度 对应DAC位宽
%%% .h file
fidd=fopen(Output_file_name_h,'wt');  %以"wt"的形式打开,\n为换行
fprintf(fidd,'#define dataLength %d\n',depth);   % # define dataLength
fprintf(fidd,'unsigned int data[%d]={\n',depth);
for i=1:depth
    data = floor(dataTD_zc(i));
    if (data>2^(width-1)-1)
        data = 2^(width-1)-1;
    elseif(data<-2^(width-1))
        data = -2^(width-1);
    end

    if (i==depth)
        fprintf(fidd,'%d\n};',data);
    else
        fprintf(fidd,'%d,\n', data);
    end
end
fclose(fidd);

% clear something to save;
clear data_f;
clear data_fh;
clear data_t;
clear maxd;
% save source for count BER
save(Output_file_name_mat);
  
  
  