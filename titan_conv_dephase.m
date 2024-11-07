% [attack,fs] = audioread('attack_vocal.wav');
% y = attack(26*fs:33.1*fs);
% [didadi,fs] = audioread('didadi_vocal.wav');
% y = didadi(78*fs:86*fs);
[Tokyo,fs] = audioread('Tokyo_vocal.wav');
y = Tokyo(29*fs:34*fs);
% [NY,fs] = audioread('verse1_vocal.wav');
% y = NY(0.5*fs:5.5*fs);
% [speech,fs] = audioread('speech.wav');
% y = speech(56*fs:62*fs);

c1 = 80; %80
c2 = 350; %300
c3 = 1000;

ny = length(y);
y = y/max(abs(y));
nh = 1024;
input = [zeros(1,nh/4) y zeros(1,nh)];
y_fof = zeros(1,length(input)*10);
n = length(input);
sst = nh/4 + 1;
window = chebwin(nh);
Phi = 0;

while sst < ny-nh/4
    sample = input(sst-nh/4:sst+nh*3/4-1).*window';
    S = 20*log10(abs(fft(sample,fs))+0.00001);
    S_0db = S - max(S(1:8000));
    S_linear = abs(fft(sample,fs));
    
%     figure
%     n = fs;
%     fn = [0:1/n:1-1/n]*fs;
%     plot(fn,S_0db)
%     xlim([0 8000])
    
    order = 300;
    num = round(8000/order);
    points = zeros(1,2*num);
    for i = 1:num
        piece = S_0db((i-1)*order+1:i*order);
        points(i) = find(piece == max(piece),1) + (i-1)*order;
        points(i+num) = find(piece == min(piece),1) + (i-1)*order;
    end
    points = unique(points);
    y_p = S_0db(points);
    
    o = 6;
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
    
%     hold on
%     scatter(points(ma),S_0db(points(ma)))
    
    order = 500;
    num = round(8000/order);
    points2 = zeros(1,num);
    for i = 1: num
        piece = S_0db((i-1)*order+1:i*order);
        points2(i) = find(piece == max(piece),1) + (i-1)*order;
    end
    points2 = [1 points2];
    y_p = S_0db(points2);
    
%     hold on
%     scatter(points2,S_0db(points2),'*')

    pt = unique([points2 points(ma)]);
    M =  [1 points(ma)];
    yM = S_0db(M);
    Parameter = [];
    if length(M) > 3
        Parameter = zeros(4,ceil((length(M)-2)/2)); % center frequency; amplitude; bandwidth; skirt width
        c = 0;
        for i = 2:length(M)-1
            if yM(i) > yM(i-1) && yM(i) > yM(i+1)
                c = c + 1;
                Parameter(1,c) = M(i);
%                 Parameter(1,c) = M(i)*1.5; % pitch up
%                 Parameter(1,c) = M(i)/3*2; % pitch down
                
%                 Parameter(2,c) = log10(norm(S_linear(M(i-1):M(i+1)))/sqrt(M(i+1)-M(i-1)+1)+1);
                Parameter(2,c) = norm(S_linear(M(i-1):M(i+1)))/sqrt(M(i+1)-M(i-1)+1);
%                 Parameter(2,c) = log10(S_linear(M(i))+1);

                    st = find(pt == M(i-1));
                    ed = find(pt == M(i+1));
                    pt_band = pt(st:ed);
                    maxima = find(pt_band == M(i));
                    peak = S_0db(pt_band(maxima));
                    y_band = S_0db(pt_band);
                    fr = find(y_band(1:maxima) < peak - 6,1,'last');
                    nt = find(y_band(maxima:end) < peak - 6,1)+maxima-1;
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
                    Parameter(3,c) = (6*(bw_sample(2) - bw_sample(1))/(bw(2) - bw(1)) + 6*(bw_sample(3) - bw_sample(2))/(bw(2) - bw(3)))/2;
%                     Parameter(3,c) = min(Parameter(3,c),-log(0.1)*Parameter(1,1));
                    if Parameter(3,c) < 0
                        error('not positive bandwidth')
                    end
                Parameter(4,c) = max(max(M(i+1) - M(i),M(i) - M(i-1))*2,1500);

        %         vq = interp1(pt_band,S_0db(pt_band),M(i-1):M(i+1),'spline');
        %         hold on
        %         plot(M(i-1):M(i+1),vq)
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
     
        if col > 0
            for i = 1:col  % min(col,4)
                beta = Parameter(4,i);
                t = pi/beta;
                alpha = min(Parameter(3,i),c1)*pi;
                f = Parameter(1,i);
                omega = 2*pi*f;
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
                
                if Parameter(3,i) > c2 && Parameter(3,i) <= c3
                    s = conv(s,rand(1,length(s)+1)*2-1);
                elseif Parameter(3,i) > c3
                    s = conv(s,rand(1,length(s)+1)*2-1);
                    s = conv(s,rand(1,length(s)+1)*2-1);
                end
                
                s = s/max(abs(s))*Parameter(2,i);
                y_fof(sst:sst+length(s)-1) = y_fof(sst:sst+length(s)-1)+s;
%                 y_fof(round(sst/2):round(sst/2)+length(s)-1) = y_fof(round(sst/2):round(sst/2)+length(s)-1)+s;
%                 y_fof(round(sst*2):round(sst*2)+length(s)-1) = y_fof(round(sst*2):round(sst*2)+length(s)-1)+s;
            end
           
        end
    end
    
    if ~isempty(Parameter)
        f0 = Parameter(1,1);
        step = round(fs/f0);
        sst = sst + step;
        Phi = mod((step/fs*f0 - floor(step/fs*f0))*2*pi + Phi,2*pi);
    else
        sst = sst + nh;
        %Phi = 0;
        
        
%         noise_n = randi([100,500]);
%         gau = randn(1,noise_n+1);
%         noise = filter(1,[1 -0.9],gau);
%         noise = diff(noise);
%         deviation = randi(noise_n)-round(noise_n/2);
%         y_fof(sst-round(noise_n/2)+deviation:sst-round(noise_n/2)+deviation+noise_n-1) = y_fof(sst-round(noise_n/2)+deviation:sst-round(noise_n/2)+deviation+noise_n-1) + noise.*hann(noise_n)'*noise_gain;
    end
    
%     y_p = S_0db(points);
%     f0 = 90;
%     for i = 2:length(points)
%         if y_p(i) > y_p(i-1) && y_p(i) > y_p(i+1)
%            f0 = points(i);
%            break
%         end
%     end
%     sst = sst + round(fs/f0);
end

%     hold on
%     S = 20*log10(abs(fft(y_fof(1:fs))));
%     S_0db = S - max(S);
%     plot(fn,S_0db)
% y_fof = y_fof(nh/4:ny+nh/4-1);
y_fof = y_fof(find(y_fof));
y_fof = y_fof/max(abs(y_fof));

% figure
% window = 1024;
% [S,fc,t1] = melSpectrogram(y',fs, ...
%                    "Window",hann(window,'periodic'),...
%                    "OverlapLength",window/2,...
%                    'NumBands',128, ...
%                    'FrequencyRange',[0,8000]);
% S = 20*log10(S+0.0001);
% S = S/max(max(abs(S)))+1;
% imagesc(t1, fc, S);
% set(gca, 'YDir', 'normal'); % flip the Y Axis so lower frequencies are at the bottom
% title('Spectrogram')
% colorbar
% 
% figure
% window = 1024;
% [S,fc,t1] = melSpectrogram(y_fof',fs, ...
%                    "Window",hann(window,'periodic'),...
%                    "OverlapLength",window/2,...
%                    'NumBands',128, ...
%                    'FrequencyRange',[0,8000]);
% S = 20*log10(S+0.0001);
% S = S/max(max(abs(S)))+1;
% imagesc(t1, fc, S);
% set(gca, 'YDir', 'normal'); % flip the Y Axis so lower frequencies are at the bottom
% title('Spectrogram')
% colorbar

% tn = (0:ny-1)/fs;
% figure
% plot(tn,y)
% hold on 
% plot(tn,y_fof')

% figure
% n = 1024;
% fn = [0:1/n:1-1/n]*fs;
% S = 20*log10(abs(fft(y_fof(2.4*fs:2.4*fs+1023))));
% S_0db = S - max(S);
% plot(fn,S_0db)