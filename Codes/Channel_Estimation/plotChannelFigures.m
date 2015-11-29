%用于画第二章中的信道图
clear all;
close all;
clc;

%导入信道时域数据
load Channel.mat;
linewidth = 1.3;
%画出时域冲激响应图
plot([0:55]*5,h_max1_norm,'LineWidth',linewidth);
grid on;
xlabel('Propagation time (ns)');
ylabel('Normalized impluse response h(t) (s^{-1})');

%画出频率冲激响应图
figure(2);
H = fft(h_max1_norm,400);
plot(0:150,20*log10(abs(H(1:151)))-max(20*log10(abs(H))),'LineWidth',linewidth);
grid on;
set(gca,'xtick',(0:30:150))
xlabel('Frequency (MHz)');
ylabel('Normalized gain (dB)');

figure(3);
H = fft(h_max1_norm,400);
plot(0:150,angle(H(1:151))/pi*180,'LineWidth',linewidth);
grid on;
set(gca,'xtick',(0:30:150))
xlabel('Frequency (MHz)');
ylabel('Phase (degree)');
