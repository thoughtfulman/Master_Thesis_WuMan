 clear all;
 close all;
 clc;
 
 load ChannelExample.mat
 linewidth = 1.3;
 t = 3:3:48;
 plot(t,h_Exp_profile+1-max(h_Exp_profile),'-+r','LineWidth',linewidth);
 hold on;
 plot(t,h_Ceil_profile+1-max(h_Ceil_profile),'-^k','LineWidth',linewidth);
 plot(t,h_profile+1-max(h_profile),'-*b','LineWidth',linewidth);
 stem(t,h_Example+1-max(h_Example),'LineWidth',linewidth);
 grid on;
 xlabel('Propagation Time (ns)');
 ylabel('Normalized Implus Response h(t) (s^{-1})');
 axis([0,48,0,1]);
 legend('Profile of exponential decay model','Profile of ceiling bounce model',...
     'Profile of the mixed model','A channel example');
 
 H_Exp_profile = fft(h_Exp_profile,300);
 H_Ceil_profile = fft(h_Ceil_profile,300);
 H_profile = fft(h_profile,300);
 H_Example = fft(h_Example,300);
 
 linewidth = 1;
 nStep = 5;
 figure(2)
 plot(0:nStep:149,20*log10(abs(H_Exp_profile(1:nStep:150)))-...
     max(20*log10(abs(H_Exp_profile(1:150)))),...
    '-+r','LineWidth',linewidth);
 hold on; 
 plot(0:nStep:149,20*log10(abs(H_Ceil_profile(1:nStep:150)))-...
     max(20*log10(abs(H_Ceil_profile(1:150)))),...
    '-^k','LineWidth',linewidth);
 plot(0:nStep:149,20*log10(abs(H_profile(1:nStep:150)))-...
     max(20*log10(abs(H_profile(1:150)))),...
     '-*b','LineWidth',linewidth);
 
 plot(0:nStep:149,20*log10(abs(H_Example(1:nStep:150)))-...
     max(20*log10(abs(H_Example(1:150)))),...
     '-<m','LineWidth',linewidth);
 ylabel('Normalized Gain (dB)');
 xlabel('Frequency (MHz)');
  legend('Profile of exponential decay model','Profile of ceiling bounce model',...
     'Profile of the mixed model','A channel example');
 grid on;