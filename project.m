Tp = 0.003;
To = 0.008;
Tc = 0.005;
Te = 0.004;
Ta = 3.3e-04;

Ee = 1.1;
wg = pi/Tp;

% Newton iteration
e = 2000;
for i = 1:10000
    y = 1-exp(-e*(Tc-Te)) - e*Ta;
    dy = exp(-e*(Tc-Te))*(Tc-Te)-Ta;
    if abs(y) < 1e-06
        break
    else
        e = e - y/dy;
    end
end
eps = e;
G2 = -Ee/(eps*Ta)*((1-exp(-eps*(Tc-Te)))/eps - exp(-eps*(Tc-Te))*(Tc-Te));

alpha = 300;
for i = 1:10000
    y = (1/wg*(1-exp(alpha*Te)*cos(wg*Te)) + alpha./wg^2.*exp(alpha*Te)*sin(wg*Te)) - G2*(1+alpha.^2/wg^2).*(exp(alpha*Te)*sin(wg*Te))/Ee;
    dy = -Te/wg*exp(alpha*Te)*cos(wg*Te) + (1+alpha*Te)/wg^2*exp(alpha*Te)*sin(wg*Te) - G2/Ee*(exp(alpha*Te)*sin(wg*Te))*(2*alpha/wg^2 + (1+alpha^2/wg^2)*Te);
    
    if abs(y) < 1e-06
        break
    else
        alpha = alpha - y/dy;
    end
end

E0 = -Ee/(exp(alpha*Te)*sin(wg*Te));

fs = 44100;
t = Tc;
tn = 0:1/fs:t-1/fs;
g1 = E0*exp(alpha*tn(1:round(Te*fs))).*sin(wg*tn(1:round(Te*fs)));
g2 = -Ee/(eps*Ta)*(exp(-eps*(tn(round(Te*fs)+1:end)-Te)) - exp(-eps*(Tc-Te)));
g = [g1 g2];

% generate sequence
ni = 250;

g_p = [g zeros(1,round((To-Tc)*fs))];
g_sq = [];
for i = 1:ni
    g_sq = [g_sq g_p];
end

% g_sq1 = [];
% for i = 1:ni
%     Tr = rand/2*0.01;
%     g_p = [g zeros(1,round(Tr*fs))];
%     g_sq1 = [g_sq1 g_p];
% end
% 
g_sq2 = [];
f = 10;
Tsin = 0.004*abs(sin(2*pi*f*[1:ni]/fs));
for i = 1:ni
    g_p = [g zeros(1,round(Tsin(i)*fs))];
    g_sq2 = [g_sq2 g_p];
end
% 
% g_sq3 = [];
% f = 100;
% Tcos = 0.004*abs(cos(2*pi*f*[1:ni]/fs));
% for i = 1:ni
%     g_p = [g zeros(1,round(Tcos(i)*fs))];
%     g_sq3 = [g_sq3 g_p];
% end
% 
% g_sq4 = [];
% Tlu = linspace(0,1,ni)*0.004;
% for i = 1:ni
%     g_p = [g zeros(1,round(Tlu(i)*fs))];
%     g_sq4 = [g_sq4 g_p];
% end
% 
% g_sq5 = [];
% Tld = linspace(1,0,ni)*0.004;
% for i = 1:ni
%     g_p = [g zeros(1,round(Tld(i)*fs))];
%     g_sq5 = [g_sq5 g_p];
% end
% 
% g_sq6 = [];
% a = 0.005;
% Tep = (exp(a*[0:ni-1])-1)*0.004;
% for i = 1:ni
%     g_p = [g zeros(1,round(Tep(i)*fs))];
%     g_sq6 = [g_sq6 g_p];
% end
% 
% g_sq7 = [];
% a = 0.01;
% Tep = exp(-a*[0:ni-1])*0.004;
% for i = 1:ni
%     g_p = [g zeros(1,round(Tep(i)*fs))];
%     g_sq7 = [g_sq7 g_p];
% end

% n = length(g_p);
% tn = [0:n-1]/fs;
% fn = [0:1/n:1-1/n]*fs;
% figure
% plot(tn, g_p)
% figure
% plot(fn, abs(fft(g_p)))
% xlim([0 5e+03])

Ra = Ta/To;
Rg = To/(2*Tp);
Rk = (Te - Tp)/Tp;
Rd = (1/0.11)*(0.5+1.2*Rk)*(Rk/(4*Rg)+Ra);
% Rap = (-1+4.8*Rd)/100;
% Rkp = (22.4+11.8*Rd)/100;
% Rgp = Rkp/(4*((0.11*Rd/(0.5+1.2*Rkp))-Rap));
% To = 0.008;
% Tc = 0.005;
% Ta = Rap*To;
% Tp = To/(2*Rgp);
% Te = Rkp*Tp + Tp;
% OQ = Te/To;