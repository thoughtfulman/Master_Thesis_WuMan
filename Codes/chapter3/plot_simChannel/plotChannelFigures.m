%用于画第二章中的信道图
clear all;
close all;
clc;

%导入信道时域数据
load Channel.mat;
linewidth = 1.3;
%画出时域冲激响应图
N = 19;
plot((0:N)*5,h_power_norm,'-*','LineWidth',linewidth);
grid on;
xlabel('Propagation time (ns)');
ylabel('Normalized impluse response h(t) (s^{-1})');

%画出频率冲激响应图
figure(2);
H = fft(h_power_norm,400);
plot(0:100,20*log10(abs(H(1:101))/max(abs(H))),'LineWidth',linewidth);
grid on;
set(gca,'xtick',(0:20:100))
xlabel('Frequency (MHz)');
ylabel('Normalized gain (dB)');

