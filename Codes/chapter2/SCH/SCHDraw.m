clear all;
close all;
clc;

load SCH.mat;

subplot(3,1,1);
stem([0:N-1],x_original,'b');
axis([0,31,-2,2]);
% title('Original optical OFDM time domain signal ');
title('(a)Bipolar optical OFDM time domain signals');
xlabel('n');
ylabel('x(n)');


subplot(3,1,2);
x14 = x_0clip(14);
x_0clip(13) = 0;
stem(0:12, x_0clip(1:13),'b');
hold on;
stem(14:28, x_0clip(15:29),'b');
stem(30:31, x_0clip(31:32),'b');
line([0,35],[x_clip_level,x_clip_level],'LineStyle','--','Color','b');
stem(13,x14,'--sr','MarkerFaceColor','g');
stem(29,0,'--sr','MarkerFaceColor','g');
axis([0,31,0,2]);
% title('ACO-OFDM time domain signal ');
title('(b)ACO-OFDM time domain signals');
xlabel('n');
ylabel('x_c(n)');

text(-1.3,x_clip_level+0.1,' \eta_c')

subplot(3,1,3);
x14 = x_0clip(14);
x_0clip(13) = 0;
stem(0:12, x_0clip(1:13),'b');
hold on;
stem(14:28, x_0clip(15:29),'b');
stem(30:31, x_0clip(31:32),'b');
line([0,35],[x_clip_level,x_clip_level],'LineStyle','--','Color','b');
stem(13,x_upclip(14),'--sr','MarkerFaceColor','g');
stem(29,x_upclip(30),'--sr','MarkerFaceColor','g');
axis([0,31,0,2]);
% title('ACO-OFDM time domain signal ');
title('(c)Recoverable upper clipped ACO-OFDM time domain signals');
xlabel('n');
ylabel('x_c(n)');

text(-1.3,x_clip_level+0.1,' \eta_c')


% subplot(3,1,3);
% stem(0:8, x_upclip(1:9),'k');
% hold on;
% stem(10:24, x_upclip(11:25),'k');
% stem(26:31, x_upclip(27:32),'k');
% line([0,35],[x_clip_level,x_clip_level],'LineStyle','--','Color','k');
% stem(9,x_upclip(10),'--sk');
% stem(25,x_upclip(26),'--sk');
% text(-1.3,x_clip_level+0.1,' \eta_c')
% axis([0,31,0,2]);
% % title('Recoverable up-clipped ACO-OFDM time domain signal ');
% title('(c)Recoverable upper clipped ACO-OFDM time domain signals');
% xlabel('n');
% ylabel('x_{uc}(n)');