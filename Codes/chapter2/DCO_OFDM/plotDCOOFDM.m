clear all;close all;clc;

hMod = comm.RectangularQAMModulator('ModulationOrder',16,...
'BitInput',false,'NormalizationMethod','Average power');

N = 4;

data = [7,9,8,7]';
dataMod = zeros(2*N,1);
dataMod(1:N) = step(hMod,data);
dataMod(N+1) = imag(dataMod(1));
dataMod(1) = real(dataMod(1));
dataMod(N+2:2*N) = conj(flipud(dataMod(2:N)));

subplot(3,1,1);
stem([0:2*N-1],abs(dataMod));
axis([-1,8,0,1.8]);
biasX = -0.2;
biasY = 0.3;
text(-0.2+biasX,0.3+biasY,'R\{X(0)\}');
text(1+biasX,1+biasY,'X(1)');
text(2+biasX,1.3+biasY,'X(2)');
text(3+biasX,0.4+biasY,'X(3)');
text(4-0.2+biasX,0.3+biasY,'I\{X(0)\}');
text(5+biasX,0.4+biasY,'X*(3)');
text(6+biasX,1.3+biasY,'X*(2)');
text(7+biasX,1+biasY,'X*(1)');
xlabel('(a)k (frequency domain subcarriers)');
ylabel('$\tilde{X}(k)$','Interpreter','latex');

%%
t = ifft(dataMod);
subplot(3,1,2);
stem([0:2*N-1],t);
axis([-1,8,-0.6,0.6]);
xlabel('(b) n (time domain samples)');
ylabel('x(n)');
%%
tb = t + 0.5;
subplot(3,1,3);
stem([0:2*N-1],tb);
axis([-1,8,0,1.0]);
line([-1,8],[0.5,0.5],'LineStyle','--');
text(-1,0.6,'bias')
xlabel('(c) n (time domain samples)');
ylabel('x(n)+bias');