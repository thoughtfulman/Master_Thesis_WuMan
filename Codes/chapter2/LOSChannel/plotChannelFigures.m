%���ڻ��ڶ����е��ŵ�ͼ
clear all;
close all;
clc;

%�����ŵ�ʱ�����
load Channel.mat;
linewidth = 1.3;
%����ʱ��弤��Ӧͼ
plot([0:55]*5,h_max1_norm,'LineWidth',linewidth);
grid on;
xlabel('Propagation Time (ns)');
ylabel('Normalized Impluse Response h(t) (s^{-1})');

%����Ƶ�ʳ弤��Ӧͼ
figure(2);
H = fft(h_max1_norm,400);
plot(0:150,20*log10(abs(H(1:151)))-max(20*log10(abs(H))),'LineWidth',linewidth);
grid on;
set(gca,'xtick',(0:30:150))
xlabel('Frequency (MHz)');
ylabel('Normalized Gain (dB)');

figure(3);
H = fft(h_max1_norm,400);
plot(0:150,angle(H(1:151))/pi*180,'LineWidth',linewidth);
grid on;
set(gca,'xtick',(0:30:150))
xlabel('Frequency (MHz)');
ylabel('Phase (degree)');
