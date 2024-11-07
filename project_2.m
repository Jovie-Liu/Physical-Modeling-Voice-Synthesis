% ni = 200;
% g_sq = [];
% Rd_v = linspace(1,2.5,ni);
% f = 10;
% Tsin = 0.008+0.002*abs(sin(2*pi*f*[1:ni]/fs));
% for j = 1:ni
%     Rd = Rd_v(ni);
    Rd = 0.85;
    Rap = (-1+4.8*Rd)/100;
    Rkp = (22.4+11.8*Rd)/100;
    Rgp = Rkp/(4*((0.11*Rd/(0.5+1.2*Rkp))-Rap));
    To = 0.008;
    Ta = Rap*To;
    Tp = To/(2*Rgp);
    Te = Rkp*Tp + Tp;
    OQ = Te/To;
    Tc = min(Te + Ta*2,To);

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

%     n = length(g);
%     tn = [0:n-1]/fs;
%     fn = [0:1/n:1-1/n]*fs;
%     figure
%     plot(tn, g)
%     figure
%     plot(tn, cumsum(g))

    % generate sequence
ni = 250;

g_p = [g zeros(1,round((To-Tc)*fs))];
g_sq = [];
for i = 1:ni
    g_sq = [g_sq g_p];
end

% g_sq2 = [];
% f = 10;
% Tsin = 0.004*abs(sin(2*pi*f*[1:ni]/fs));
% for i = 1:ni
%     g_p = [g zeros(1,round(Tsin(i)*fs))];
%     g_sq2 = [g_sq2 g_p];
% end
%     g_p = [g zeros(1,round((To-Tc)*fs))];
%     g_sq = [g_sq g_p];
% end
