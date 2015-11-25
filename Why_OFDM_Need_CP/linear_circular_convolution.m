clear all;
close all;
clc;

x = 1:8;
h = [1,2];

% linear convolution
y_l = conv(x,h);
% check the frequency domain
X16 = fft(x,16);
H16 = fft(h,16);
Y16 = fft(y_l,16);
E16 = Y16-X16.*H16;

figure(1);
plot(abs(E16));

% circular convolution
y_c = cconv(x,h,8);
X8 = fft(x,8);
H8 = fft(h,8);
Y8 = fft(y_c,8);
E8 = Y8-X8.*H8;
figure(2);
plot(abs(E8));

% use linear convolution to calculate circular convolution
% CP is needed
x_cp = [x(7:8),x];
y_l_cp = conv(x_cp,h);
figure(3);
stem(y_l_cp);
hold on;
stem(3:10,y_c,'r*');
legend('linear conv with cp','circular conv');
