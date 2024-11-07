ng = length(g);
ngp = length(g_p);
n = length(g_sq);
nslc = ngp-ng;
ni = n/ngp;
G = zeros(1,n+ngp);
for i = 1:ni
    L = randi([min(nslc,ng),ngp]);
    hanning = hann(L)';
    gaussian = randn(1,L);
    lag = 0;
    G(ngp*(i-1)+ng-round(L/2)+lag:ngp*(i-1)+ng-round(L/2)+lag+L-1) = G(ngp*(i-1)+ng-round(L/2)+lag:ngp*(i-1)+ng-round(L/2)+lag+L-1) + hanning.*gaussian;
end
g_sqnoise = g_sq + 0.1*G(1:n);
% plot(g_sqnoise)

% all gaussian noise
Gaussian = randn(1,n+ngp*2);
for i = 1:ni
    L = randi([ngp,max(nslc,ng)*2]);
    hanning = hann(L)';
    lag = randi(50);
    Gaussian(ngp*(i-1)+ng-round(L/2)+lag:ngp*(i-1)+ng-round(L/2)+lag+L-1) = Gaussian(ngp*(i-1)+ng-round(L/2)+lag:ngp*(i-1)+ng-round(L/2)+lag+L-1).*hanning;
end
g_sqnoise1 = g_sq + 0.1*Gaussian(1:n);
% plot(g_sqnoise1)

noise = filter(1,[1 -0.9],Gaussian);
noise = diff(noise);
g_sqnoise2 = g_sq + 0.1*noise(1:n);
% plot(g_sqnoise2)

% gaussian noise for varied sequence
g_sq2ns = [];
f = 10;
Tsin = 0.004*abs(sin(2*pi*f*[1:ni]/fs));
for i = 1:ni
    g_p = [g zeros(1,round(Tsin(i)*fs))];
    g_sq2 = [g_sq2 g_p];
end