% [attack,fs] = audioread('attack_vocal.wav');
% y = attack(26*fs:39*fs);
% [Tokyo,fs] = audioread('Tokyo_vocal.wav');
% y = Tokyo(29*fs:34*fs);
[NY,fs] = audioread('verse1_vocal.wav');
y = NY(0.5*fs:5.5*fs);
% [didadi,fs] = audioread('didadi_vocal.wav');
% y = didadi(78*fs:86*fs);
% [speech,fs] = audioread('speech.mp3');
% y = speech(56*fs:62*fs);
y = y/max(abs(y));

c1 = 80;
c2 = 350;
c3 = 1000;
freq = 11000;

nh = 1024;
wid = chebwin(nh);
sample = y(1.9*fs:1.9*fs+nh-1).*wid';
n = fs;
% figure
% tn = [0:n-1]/fs;
% plot(tn,sample)
figure
fn = [0:1/n:1-1/n]*fs;
S = 20*log10(abs(fft(sample,fs)));
S_0db = S - max(S(1:freq));
plot(fn,S_0db)
xlim([0 freq])

S_linear = abs(fft(sample,fs));

% figure
% plot(fn,S_linear)
% xlim([0 8000])

% hold on
% a = lpc(sample,20);
% [h,f] = freqz(1,a,fs,'whole',fs);
% plot(fn,20*log10(abs(h)))
% hold off

% order = 13;
% smooth = hann(order,'periodic');
% S_smt = zeros(1,n);
% S_smt(1:6) = S(1:6);
% S_smt(n-5:n) = S(n-5:n);
% for i = 7:n-6
%     S_smt(i) = S(i-6:i+6)*smooth;
% end
% plot(fn,S_smt)

order = 350;
num = round(freq/order);
points = zeros(1,2*num);
for i = 1:num
    piece = S_0db((i-1)*order+1:i*order);
    points(i) = find(piece == max(piece)) + (i-1)*order;
    points(i+num) = find(piece == min(piece)) + (i-1)*order;
end
% points = [1 points];
points = unique(points);
y_p = S_0db(points);
% hold on
% scatter(points,y_p,'*')

o = 5;
ma = [];
st = 1;
for i = 1:500
    cal = abs(y_p(st) - y_p);
    for j = st+1:num*2
        ind = find(cal(st:j) == max(cal(st:j)),1)+st-1;
        if ind + o <= num*2 && ind ~= j && ind == find(cal(st:ind+o) == max(cal(st:ind+o)),1)+st-1
            ma = [ma ind];
            st = ind;
            break
        elseif ind + o > num*2
            ma = [ma ind];
            break
        end
    end
    if ind + o > num*2 || j == num*2
        break
    end
end
hold on
scatter(points(ma),S_0db(points(ma)))

order = 550;
num = round(freq/order);
points2 = zeros(1,num);
for i = 1: num
    piece = S_0db((i-1)*order+1:i*order);
    points2(i) = find(piece == max(piece)) + (i-1)*order;
end
points2 = [1 points2];
y_p = S_0db(points2);
% hold on
% scatter(points2,S_0db(points2),'*')


pt = unique([points2 points(ma)]);
M =  [1 points(ma)];
yM = S_0db(M);
Parameter = zeros(4,ceil((length(M)-2)/2)); % center frequency; amplitude. bandwidth; skirt width
c = 0;
db = 6;
for i = 2:length(M)-1
    if yM(i) > yM(i-1) && yM(i) > yM(i+1)
        if S_linear(M(i))/max(S_linear) > 0.003
            c = c + 1;
            Parameter(1,c) = M(i);
            Parameter(2,c) = norm(S_linear(M(i-1):M(i+1)))/sqrt(M(i+1)-M(i-1)+1);%/(M(i+1)-M(i-1));

            st = find(pt == M(i-1));
            ed = find(pt == M(i+1));
            pt_band = pt(st:ed);
            maxima = find(pt_band == M(i));
            peak = S_0db(pt_band(maxima));
            y_band = S_0db(pt_band);
            fr = find(y_band(1:maxima) < peak - db,1,'last');
            nt = find(y_band(maxima:end) < peak - db,1)+maxima-1;
            if isempty(fr)
                fr = 1;
            end
            if isempty(nt)
                nt = length(y_band);
            end
            bw_i = [fr maxima nt];
            bw_sample = pt_band(bw_i);
            bw = y_band(bw_i);
            if bw(2) - bw(1) < 0 || bw(2) - bw(3) < 0
                error('not local maxima')
            end
            Parameter(3,c) = (db*(bw_sample(2) - bw_sample(1))/(bw(2) - bw(1)) + db*(bw_sample(3) - bw_sample(2))/(bw(2) - bw(3)))/2;
%             Parameter(3,c) = min(Parameter(3,c),-log(0.1)*Parameter(1,1));
            if Parameter(3,c) < 0
                error('not positive bandwidth')
            end
            Parameter(4,c) = M(i+1) - M(i-1); % max(M(i+1) - M(i),M(i) - M(i-1))*2
%             
%     %         vq = interp1(pt_band,S_0db(pt_band),M(i-1):M(i+1),'spline');
%     %         hold on
%     %         plot(M(i-1):M(i+1),vq)
        end
    end
end
[~,col] = size(Parameter);

index = [];
for i = 1:col
    if isequal(Parameter(:,i), zeros(4,1))
        index = [index i];
    end
end
Parameter(:,index) = [];
[~,col] = size(Parameter);
y_fof = zeros(1,fs);
for i = 1:col
    beta = Parameter(4,i);
    t = pi/beta;
    alpha = min(Parameter(3,i),c1)*pi;
%     alpha = 15*pi;
    f = Parameter(1,i);
    omega = 2*pi*f;
    Phi = 0;
    kn1 = [0:fs*t-1]/fs;
    s1 = (1 - cos(beta*kn1))/2.*exp(-alpha*kn1).*sin(omega*kn1+Phi);
    for ts = 1:3000
        st = exp(-alpha*t*ts);
        if st < 0.003
            break
        end
    end
    kn2 = [t*fs:t*fs*ts]/fs;
    s2 = exp(-alpha*kn2).*sin(omega*kn2 + Phi);
    s = [s1 s2];
    
%     if Parameter(3,i) > c2 && Parameter(3,i) <= c3
%         s = conv(s,rand(1,length(s)+1)*2-1);
%     elseif Parameter(3,i) > c3
%         s = conv(s,rand(1,length(s)+1)*2-1);
%         s = conv(s,rand(1,length(s)+1)*2-1);
%     end

    s = s/max(abs(s))*Parameter(2,i);
    y_fof(1:length(s)) = y_fof(1:length(s))+s;
end
% gau = randn(1,301);
% noise = filter(1,[1 -0.9],gau);
% noise = diff(noise);
% 
% y_fof(1:300) = y_fof(1:300)+noise.*hann(300)'*0.02;

hold on
fn = [0:1/n:1-1/n]*fs;
S = 20*log10(abs(fft(y_fof)));
S_0db = S - max(S);
plot(fn,S_0db)
title('Spectral Plot for One Frame')
xlabel('Frequency(Hz)')
ylabel('Amplitude(dB)')

% figure
% plot(fn,S_linear)
% xlim([0 8000])