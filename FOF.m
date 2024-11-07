c = 1;
d = 0.5;
a = 0.3;
b = 0.7;
fs = 440;
t = 1;
x = [1 zeros(1,fs*t-1)];
B = c*[1 d];
A = [1 a b];
[z,p,~] = tf2zp(B,A);
y = filter(B,A,x);
[h, w] = freqz(B,A,fs);
[b,a] = invfreqz(h,w,1,2);
% figure
% n = length(x);
% tn = 0:t/n:t-1/n;
% plot(tn,y)
figure
fn = [0:1/n:1-1/n]*fs;
% plot(fn,abs(fft(y)))
plot(w/pi*fs/2,abs(h))
% xlim([0 fs/2])

% alpha = -1/2*log(b);
% omega = angle(p(1));
% Phi = asin(sin(omega*exp(-alpha))/(d-a-cos(omega*exp(-alpha))));
% G = c/sin(Phi);
% k = 0:fs;
% s = G * exp(-alpha*k).*sin(omega*k+Phi);
% figure
% plot(real(s))
