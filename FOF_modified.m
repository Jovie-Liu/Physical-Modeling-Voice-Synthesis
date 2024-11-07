% modified FOF
t = 0.002;
beta = pi/t;
alpha = 30*pi;
f = 441;
omega = 2*pi*f;
Phi =0;
fs = 44100;
kn1 = [0:fs*t-1]/fs;
gain = 1;
s1 = gain*(1 - cos(beta*kn1))/2.*exp(-alpha*kn1).*sin(omega*kn1+Phi);

for ts = 1:5000
    st = gain*exp(-alpha*t*ts);
    if st < 0.005
        break
    end
end

kn2 = [t*fs:t*fs*ts]/fs;
s2 = gain*exp(-alpha*kn2).*sin(omega*kn2 + Phi);
s = [s1 s2];
figure
n = length(s);
% plot(s)
tn = [0:n-1]/fs;
plot(tn,s)
% hold on
% scatter(round(fs/f*2),s(round(fs/f*2)))
% figure
% fn = [0:1/n:1-1/n]*fs;
% S = 20*log10(abs(fft(s)));
% S_0db = S - max(S);
% plot(fn,S_0db)
% xlim([0 5000])

% % % S = fft(s);
% % % w = [1:round(n/2)]/(fs/2)*pi;
% % % [b,a] = invfreqz(S(1:round(n/2)),w,1,2);
% % % 
% % % [h, w] = freqz(b,a,'whole',n);
% % % figure
% % % S = 20*log10(abs(fft(h)));
% % % S_0db = S - max(S);
% % % plot(w/(2*pi)*fs,S)
% % % xlim([0 5000])

% s = conv(s,rand(1,length(s)+1));
% s = conv(s,s);
% s = s(10:end-50);
% n = length(s);
% figure
% tn = [0:n-1]/fs;
% plot(s)
% 
% figure
% fn = [0:1/n:1-1/n]*fs;
% S = 20*log10(abs(fft(s)));
% S_0db = S - max(S);
% plot(fn,S_0db)
% xlim([0 5000])



% db = -6;
% P_1 = find(S_0db(1:round(n/2)) >= db, 1);
% dif1 = abs(S_0db(P_1-1)-db)/(S_0db(P_1) - S_0db(P_1-1));
% bw1 = P_1-1 + dif1;
% 
% P_2 = find(S_0db(P_1:round(n/2)) <= db, 1)-1+P_1;
% dif2 = abs(S_0db(P_2-1)-db)/abs(S_0db(P_2) - S_0db(P_2-1));
% bw2 = P_2-1 + dif2;
% 
% bw = fs/n*(bw2 - bw1);

% vq = interp1([1:n]/n*fs,S_0db,1:fs,'spline');
% P_1 = find(vq(1:fs/2) >= db, 1);
% P_2 = find(vq(P_1:fs/2) <= db, 1)-1+P_1;
% bw = (P_2 - P_1)/2;