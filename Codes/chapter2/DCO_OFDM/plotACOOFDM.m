clear all;close all;clc;

hMod = comm.RectangularQAMModulator('ModulationOrder',16,...
'BitInput',false,'NormalizationMethod','Average power');

N = 4;

dataMod = zeros(2*N,1);
dataMod(2) = 2*(0.0211983690570943 + 0.678591264876017i);
dataMod(4) = 1.8*(1.13416565020406 + 0.245179359324752i);
dataMod(N+1) = imag(dataMod(1));
dataMod(1) = real(dataMod(1));
dataMod(N+2:2*N) = conj(flipud(dataMod(2:N)));

subplot(3,1,1);
stem([0:2*N-1],abs(dataMod));
axis([-1,8,0,2.8]);
text(0.8,1.8,'X(0)');
text(2.8,2.5,'X(1)');
text(4.8,2.5,'X*(1)');
text(6.8,1.8,'X*(0)');
xlabel('(a)k (frequency domain subcarriers)');
ylabel('$\tilde{X}(k)$','Interpreter','latex');

%%
t = ifft(dataMod);
subplot(3,1,2);
stem([0:4:2*N-1],t(1:4:2*N),'*b');
hold on;
stem([1:4:2*N-1],t(2:4:2*N),'ob');
stem([2:4:2*N-1],t(3:4:2*N),'sb');
stem([3:4:2*N-1],t(4:4:2*N),'^b');
axis([-1,8,-0.8,0.8]);
xlabel('(b) n (time domain samples)');
ylabel('x(n)');
%%
tb = t .* (t>0);
subplot(3,1,3);
stem(0,tb(1),'*b');
hold on;
stem(3,tb(4),'^b');
stem(5,tb(6),'ob');
stem(6,tb(7),'sb');
stem(1,0,'ob');
stem(2,0,'sb');
stem(4,0,'*b');
stem(7,0,'^b');
axis([-1,8,0,0.8]);

xlabel('(c) n (time domain samples)');
ylabel('x_r(n)');