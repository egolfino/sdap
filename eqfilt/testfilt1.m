function testfilt1
x=randn(1e4,1);
load EQF_008R   % gives h
% h=h(257:end);

figure;
fs=44100;
N=length(h);
H=fft(h);
semilogx(0:(fs/N):(fs/N*(N-1)),20*log10(H),'b'); hold on;

h=resample(h,48000,44100);


fs=48000;


[j,idxmax]=max(h);
h=h(idxmax-20:idxmax+20);
H=fft(h,82);
N=length(h);

semilogx(0:(fs/N):(fs/N*(N-1)),20*log10(H),'r');

% x=conv(x,h);

% X=fft(x);
% figure;
% plot(abs(20*log10(X)));
format long;
h=h';
[j,idxmax]=max(h);
h(idxmax:idxmax+40)
format
return
