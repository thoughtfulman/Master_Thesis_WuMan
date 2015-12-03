%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%把所有的数据都做成了浮点的格式，避免运行中报错。
% open the file.
 fid = fopen('modulation_initial_file.txt'); 
% use textscan to input for format data.
 inData1 = textscan(fid, '%s =%f',11,'commentStyle','\\');   
 % assign values in the file to the program.
 nSC = inData1{2}(1);           % number of subcarriers
 nClt = inData1{2}(2);       % number of clusters
 nCP = inData1{2}(3);           % length of CP
 nFrame = inData1{2}(4);        % frames of data
 nSignal = inData1{2}(5);       % number of OFDM signals in a frame
 FixSeedEnable = inData1{2}(6); % fix the rand generator or not 
 Seed = inData1{2}(7);          % Seed, must not a nonzero integer.
 pilotEnable = inData1{2}(8);
 nIdleLF = inData1{2}(9);     % number of idle subcarriers in low frequency end
 nIdleHF = inData1{2}(10);     % number of idle subcarriers in low frequency end
 ModOrdPowerWetInputMod =inData1{2}(11); %modulation order and power weight input method 1->define in this .txt file 0-> define in a .xlsx 
 
 nSC_clt = nSC/nClt;            %number of subcarriers in a cluster 
 inData2 = textscan(fid,'%s ={%s}',1,'commentStyle','\\');
 modOrd_str = inData2{2}{1,1};
 modOrd_str = modOrd_str(1:length(modOrd_str)-1);
 inData3 = textscan(fid,'%s ={%s}',1,'commentStyle','\\');
 powerWet_str = inData3{2}{1,1};
 powerWet_str = powerWet_str(1:length(powerWet_str)-1);
 %%% Define Output .coe file name in initial file.
 inData4 = textscan(fid,'%s ="%s"',1,'commentStyle','\\');
 input_file_name_ModPow  = inData4{2}{1,1};
 input_file_name_ModPow = input_file_name_ModPow(1:length(input_file_name_ModPow)-1);
 %%% Define Output .mat file name in initial file.
 inData5 = textscan(fid,'%s ="%s"',1,'commentStyle','\\');
 Output_file_name_mat = inData5{2}{1,1};
 Output_file_name_mat = Output_file_name_mat(1:length(Output_file_name_mat)-1);
 %%% Define Output .h file name in initial file.
 inData6 = textscan(fid,'%s ="%s"',1,'commentStyle','\\');
 Output_file_name_h = inData6{2}{1,1};
 Output_file_name_h = Output_file_name_h(1:length(Output_file_name_h)-1);
 % close the file
fclose(fid);  

if(ModOrdPowerWetInputMod)
     modOrd_t1 = textscan(modOrd_str,'%f','Delimiter',',');
     modOrd_clt = modOrd_t1{1};
     modOrd = reshape(repmat(modOrd_clt',nSC_clt,1),nSC,1);
     powerWet_t1 = textscan(powerWet_str,'%f','Delimiter',',');
     powerWet_clt = powerWet_t1{1};
     powerWet = reshape(repmat(powerWet_clt',nSC_clt,1),nSC,1);
else
     modpow_info=xlsread(input_file_name_ModPow);
     modOrd = modpow_info(:,2);
     powerWet = modpow_info(:,3);
end
nModBit = log2(modOrd);
nBit_sig = sum(log2(modOrd(nIdleLF+1:nSC-nIdleHF)));

  % file process end
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%