% LF-model voice synthesis
Tp = 0.003;
To = 0.008;
Tc = 0.005;
Te = 0.004;
Fe = 1/Te;
Ee = 1.1;
alpha = 444.9445; %443
eps = 3000;

wg = pi/Tp;
Ta = (1-exp(-eps*(Tc-Te)))/eps;
E0 = -Ee/(exp(alpha*Te)*sin(wg*Te));
fs = 44100;
t = Tc;
tn = 0:1/fs:t-1/fs;
g1 = E0*exp(alpha*tn(1:round(Te*fs))).*sin(wg*tn(1:round(Te*fs)));
g2 = -Ee/(eps*Ta)*(exp(-eps*(tn(round(Te*fs)+1:end)-Te)) - exp(-eps*(Tc-Te)));
g = [g1 g2];
% integral = (sum(g(1:end-1))+g(end)/2)/fs;
% figure
% plot(tn,g)

% G2 = -Ee/(eps*Ta)*((1-exp(-eps*(Tc-Te)))/eps - exp(-eps*(Tc-Te))*(Tc-Te));
% G1 = -Ee/(exp(alpha*Te)*sin(wg*Te))/(1+alpha^2/wg^2)*(1/wg*(1-exp(alpha*Te)*cos(wg*Te)) + alpha/wg^2*exp(alpha*Te)*sin(wg*Te));
% G = G1 + G2;
% G_sum = cumsum(g);
    
% figure
% plot(tn,G_sum)
% xlim([0 Tc])

% generate sequence
g_p = [g zeros(1,round((To-Tc)*fs))];
% g_pp = [g zeros(1,fs)];
% g_pp = [1 zeros(1,round((To-Tc)*fs))];
% n = length(g_pp);
% tn = [0:n-1]/fs;
% fn = [0:1/n:1-1/n]*fs;
% Y = fft(g_pp);
% plot(fn, abs(fft(g_pp)))
% xlim([0 fs/2])

g_sq = [];
for i = 1:125
    g_sq = [g_sq g_p];
end

% g_sq1 = [];
% for i = 1:125
%     To = randi([5 10])*0.001;
%     g_p = [g zeros(1,round((To-Tc)*fs))];
%     g_sq1 = [g_sq1 g_p];
% end

% g_sq2 = [];
% f = 10;
% for i = 1:125
%     T = 0.005*sin(2*pi*f*[1:125]/fs);
%     g_p = [g zeros(1,round(T(i)*fs))];
%     g_sq2 = [g_sq2 g_p];
% end
n = length(g_sq);
tn = [0:n-1]/fs;
fn = [0:1/n:1-1/n]*fs;
figure
plot(tn, g_sq)
figure
plot(fn, abs(fft(g_sq)))
xlim([0 fs/2])