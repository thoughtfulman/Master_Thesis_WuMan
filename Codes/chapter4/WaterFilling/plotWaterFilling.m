nSubChannel             = 16;
totalPower              = 16;           %  -20 dBm
channelStateInformation = random('rayleigh',1/0.6552,1,nSubChannel);
bandwidth               = 1e6;            %  1 MHz
noiseDensity            = 1e-5;          % -80 dBm
[Capacity PowerAllocated] = ofdmwaterfilling(...
   nSubChannel,totalPower,channelStateInformation,bandwidth,noiseDensity)