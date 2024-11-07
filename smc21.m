% 'raw' area functions as per Brad Story JASA 1996 
% sampling interval: D = 0.396 cm, corresponds to fs/2 = c/(2D) 

fvowel = ['A-bart.txt'; ...  % 1, aa
          'u-food.txt'; ...  % 2
          'i-beet.txt'; ...  % 3
          'U-foot.txt']; ... % 4

fvowel44 = ['A-bart44100.txt'; ...  % 1, aa
            'u-food44100.txt'; ...  % 2
            'i-beet44100.txt'; ...  % 3
            'U-foot44100.txt']; ... % 4

STORY = '1994_Vowels_BStory/';
faudio = ['1994-08-02_Supine_aa.wav '; ... % 1, aa: A-bart,father,hot
          '1994-08-02_Upright_aa.wav']; 

VOWEL = 1;          
AUDIO = 1;

%% calculate vocal tract TF using piecewise cylinrical model

areafile = strip(fvowel44(VOWEL, :));
display(strcat('vowel:', ' ', areafile));

Nb = 1024*4;                      % half-bandwidth (fs/2), bins
wT = pi*[0:Nb-1]'/Nb;             % half-bandwidth frequency axis, rad/sample

S = load(areafile);               % vocal tract area for vowel
M = length(S);                     % number of sections (length of S)
Nj = M - 1;                       % number of junctions
rL = sqrt((S(end)*1e-4)./pi);     % radius at opening, m (for RL)

% tube shape to be modeled: 
r = sqrt(S./pi)./sqrt(max(S)/pi); % normalized radii 
S = pi.*r.^2;                     % normalized cross-sectional area

% piecewise cylindrical model
k = (S(1:end-1) - S(2:end))./...
    (S(1:end-1) + S(2:end));
A = zeros(Nb, 4);
P = [ones(Nb, 1) zeros(Nb, 1) ...
     zeros(Nb, 1) ones(Nb, 1)];
for m = 1:Nj
    A(:,1) = exp(j*wT);
    A(:,2) = k(m)*exp(-j*wT);
    A(:,3) = k(m)*exp(j*wT);
    A(:,4) = exp(-j*wT);
    
    P = matrixMultiply(P, A./(1 + k(m)));
end

%% Vocal Tract Frequency Responses H0 (glottis) and HL (mouth)

R0 = .8;                      % reflection at glottis 
RL = -1;                      % lossless (inverting) reflection at mouth
bL = -[0.2162 0.2171 0.0545]; % coefficients for reflection at mouth
aL = [1 -0.6032 0.0910];
RLw = freqz(bL, aL, Nb);      % frequency response for reflection at mouth
TLw = freqz(bL + aL, aL, Nb); % complementary TL = 1 + RLw (for non-lossless)


HLw1 = exp(-j*wT)./...
       (P(:,1) + P(:,2).*RL - R0.*(P(:,3) + P(:,4).*RL).*exp(-2*j*wT));


GLw1 = exp(j*wT).* ...         % inverse filter frequency response 1./HLw1
       (P(:,1) + P(:,2).*RL - R0.*(P(:,3) + P(:,4).*RL).*exp(-2*j*wT));


HLw2 = TLw.*exp(-j*wT)./...
       (P(:,1) + P(:,2).*RLw - R0.*(P(:,3) + P(:,4).*RLw).*exp(-2*j*wT));


GLw2 = exp(j*wT).* ...         % inverse filter frequency response 1./HLw2 (TL omitted)
       (P(:,1) + P(:,2).*RLw - R0.*(P(:,3) + P(:,4).*RLw).*exp(-2*j*wT));


%% Calculate coefficient vectors

cN = [1; 0; k(1)*k(2); 0];
dN = [k(1); 0; k(2); 0];
for n = 3:Nj
    cNm1 = cN; % store value of cN before updating
    cN = [cN; 0; 0] + k(n)*[0; flip(eye(2*(n-1)))*dN; 0];
    dN = [dN; 0; 0] + k(n)*[0; flip(eye(2*(n-1)))*cNm1; 0];
end

%% Transfer Functions as ratio of polynomial functions B(z)/A(z)

R = [1 -1 -R0 -R0*(-1)]';                  % lossless boundary reflections
Rv = [aL; bL; -R0*aL; -R0*bL];             % freq-dependent reflection coefficients

JN = flip(eye(2*Nj));
CN = [ [cN; 0; 0; 0] [0; JN*dN; 0; 0] ...
       [0; 0; dN; 0] [0; 0; 0; JN*cN] ];

% A coefficients for HL
AM = CN*R;                                 % lossless case
AM2 = [CN; zeros(2, 4)]*Rv(:,1) + ...      % frequency-dependent case
      [zeros(1, 4); CN; zeros(1,4)]*Rv(:,2) + ...
      [zeros(2, 4); CN]*Rv(:,3);

% B coefficients for HL in lossless case
B = [zeros(1, M) prod(1+k)];               

% B coefficients for HL in frequency-dependent case: 
% B(z)*AL(z)*TL 
%   = B(z)*AL(z)*(BL(z) + AL(z))/AL(z) 
%   = B(z)*(BL(z) + AL(z)
B2T = conv(B, bL + aL);                % conv. of poly. coefficients (w TL)
B2 = conv(B, aL);                      % coefficients (w/o TL)

HLz1 = freqz(B, AM', Nb);              % lossless HL
GLz1 = freqz(AM', B, Nb);              % inverse filter (lossless)

HLz2 = freqz(B2T, AM2', Nb);           % HL with reflection functions
HLz2b = freqz(B2, AM2', Nb).*TLw;      % debug: check is equal to HLz2
GLz2 = freqz(AM2', B2, Nb);            % TL omitted b/c ill-conditioned in its inverse

%GLcy = [GLcy; flipud(conj(GLcy))];

%% Load audio files 

brad = audioread(strcat(STORY, faudio(1,:))); % Brad Story
Udev = audioread('Udev5.wav');                % devansh5.m
ydev = audioread('ydev5.wav');

ix = 1000 + [1:Nb];
ysnip = ydev(ix);
%bsnip = hanning(Nb).*brad(19000+[1:Nb], 1);
bsnip = brad(19000+[1:Nb], 1);

GL = [GLz2; flipud(conj(GLz2))];

% Estimate flow using inverse filter which, without TL, is really derivative
Udhatf = real(ifft(fft(ysnip, 2*Nb).*GL)); % frequency domain implementation
Udhatf = Udhatf(1:Nb);

Udhatt = filter(AM2', B2(M+1:end), ysnip); % time domain (remove pure delay in denominator)     
Udhatt = [Udhatt(M+1:end); zeros(M, 1)];   % restore pure delay with sample advance

% Estimate true flow by applying inverse transmission filter
Uhatf = filter(aL, bL + aL, Udhatf);
Uhatfn = (Uhatf - min(Uhatf))./(max(Uhatf) - min(Uhatf)); % normalized

Uhatt = filter(aL, bL + aL, Udhatt);
Uhattn = (Uhatt - min(Uhatt))./(max(Uhatt) - min(Uhatt)); % normalized

% Compare to flow from devansh5.m
plot([Udev(ix)./max(Udev(ix)) Uhatfn Uhattn]); pause; 


% Flow from Brad Story audio 
UdhatB = real(ifft(fft(bsnip, 2*Nb).*GL));                % flow derivative
UhatB = filter(aL, bL + aL, UdhatB(1:Nb));                % apply inverse TL
UhatBn = (UhatB - min(UhatB))./(max(UhatB) - min(UhatB)); % normalize

plot(UdhatB); 
title('Flow Derivative (Brad Story)'); 
pause; 
plot(UhatB);
title('Flow (Brad Story)');

%% TODO: add standard LPC flow estimation


%% Estimate RL (IGNORE FOR NOW)

HL = HLcy;
RLhat = -(HL.*( P(:,1) - R0*P(:,3).*exp(-2*j*wT) ) - exp(-j*wT))./ ...
        (HL.*( P(:,2) - R0*P(:,4).*exp(-2*j*wT) ) - exp(-j*wT));

theta1 = angle(P(:,1));
theta2 = angle(P(:,2));
G1oG2 = abs(P(:,1))./abs(P(:,2));

b0 = G1oG2.*exp(-j*theta2).*exp(j.*theta1); %P(:,1)./P(:,2);
b1 = 1./(HL.*P(:,2));
b2 = R0.*exp(-j.*2*theta2); %R0.*conj(P(:,2))./P(:,2);
a2 = G1oG2.*R0.*exp(-j.*theta2).*exp(-j.*theta1); %R0.*conj(P(:,1))./P(:,2);
a2 = b0.*R0.*exp(-j*2*theta1);

if 0
b0 = P(:,1)./P(:,2);
b1 = 1./(HL.*P(:,2));
b2 = R0.*exp(-j.*2*theta2); %R0.*conj(P(:,2))./P(:,2);
a2 = R0.*conj(P(:,1))./P(:,2);
end
RLhat2 = -(b0 - b1.*exp(-j*wT) - b2.*exp(-2*j*wT))./...
         (1 - b1.*exp(-j*wT) - a2.*exp(-2*j*wT));


if 0

k = 88200*wT/340;
xsi = 2.2;       % scalar, near one, that alters transition
%behaviour
xsi = 1.2;
ra = .01; % radius, m (1 cm)
RL = -1./(1 + 2*j*k*ra/xsi);

 %f = (fs/2)*[0:nbins]'/nbins;
f = 88200*[0:Nbins-1]'/Nbins;
f0 = 340/(2*pi*ra);
ka = 1.15*f/f0; 
%ka = f/f0;

%form Levine-Schwinger (theoretical) reflection filter
kaX = [0.01 0.03 0.1 0.3  1.0 1.55 3.0 4.10 5.0 6.44 7.82 10.0]';
X = exp(-[5.10 4.02 2.84 1.76 0.62 0.35 1.14 2.12 1.92 2.71 2.28 2.97])';

kaR = [0.01 0.03  0.1 0.3 1.0 2.36 3.0  4.10 5.0  5.66 7.10 8.90 10.0]';
R = exp(-[10.63 8.41 6.11 3.82 1.47 0.00 -0.065 0.11 -0.052 0.00 0.079 -0.032 0.00])';

%Xi = interp1([0; kaX], [0; X], ka);
Xi = interp1([0; kaX], [0; X], ka);
Ri = spline([0; kaR], [0; R], ka);
rhoYt = (Ri + sqrt(-1)*Xi - 1) ./ (Ri + sqrt(-1)*Xi + 1);

semilogx([0:Nbins-1]*88200/Nbins/1000, 20*log10(abs([RL RLm rhoYt]))); 
grid; set(gca, 'XLim', [0.1 10], 'YLim', [-40 5]);

end


%% multiply two 2X2 matrices for which entries are vectors
function y = matrixMultiply(A1, A2)

y(:,1) = A1(:,1).*A2(:,1) + A1(:, 2).*A2(:,3);
y(:,2) = A1(:,1).*A2(:,2) + A1(:, 2).*A2(:,4);
y(:,3) = A1(:,3).*A2(:,1) + A1(:, 4).*A2(:,3);
y(:,4) = A1(:,3).*A2(:,2) + A1(:, 4).*A2(:,4);

end

    
