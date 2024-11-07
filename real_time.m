fs = 44100;
frame = 1024;

afr = dsp.AudioFileReader('didadi_vocal.wav','ReadRange',[78*fs 86.5*fs],'SamplesPerFrame',frame);
% afr = dsp.AudioFileReader('speech.mp3','ReadRange',[56*fs 62.3*fs],'SamplesPerFrame',frame);
% afr = dsp.AudioFileReader('verse1_vocal.wav','ReadRange',[0.5*fs 5.5*fs],'SamplesPerFrame',frame);
% afr = dsp.AudioFileReader('Tokyo_vocal.wav','ReadRange',[29*fs 34*fs],'SamplesPerFrame',frame);

% deviceReader = audioDeviceReader( ...
%     'SamplesPerFrame',frame, ...
%     'SampleRate',44100);
% setup(deviceReader)

adw = audioDeviceWriter('SampleRate', afr.SampleRate,'SupportVariableSizeInput',true);

c1 = 80; %80
c2 = 350; %300
c3 = 1000;

buffer = zeros(1,fs);
sample_buffer = zeros(1,fs);
pos = 1;
sst = 1;
samples_read = 0;
% tic
ratio = [0.8 0.9 1.2 0.5];
while ~isDone(afr) || pos > 300
% while toc < 10
    if sst+frame-1 > samples_read && ~isDone(afr)
        audio = afr(); %[1024,2]
%         audio = deviceReader();
        audio = (audio(:,1)+audio(:,2))/2;
%         audio = audio(:,1);
        if pos > 0
            sample_buffer(pos:pos+frame-1) = audio';
        end
        pos = pos+frame;
        samples_read = samples_read + frame;
        
%         audio_acc = afr_acc();
%         audio_acc = (audio_acc(:,1)+audio_acc(:,2))/2;
%         sample_buffer(2,pos:pos+frame-1) = audio_acc';
    end
    
    sample = sample_buffer(1:frame);
    S = 20*log10(abs(fft(sample,fs))+0.00001);
    S_0db = S - max(S(1:8000));
    S_linear = abs(fft(sample,fs));
    
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
    
    order = 500;
    num = round(8000/order);
    points2 = zeros(1,num);
    for i = 1: num
        piece = S_0db((i-1)*order+1:i*order);
        points2(i) = find(piece == max(piece),1) + (i-1)*order;
    end
    points2 = [1 points2];
    y_p = S_0db(points2);
    
    pt = unique([points2 points(ma)]);
    M =  [1 points(ma)];
    yM = S_0db(M);
    Parameter = [];
    if length(M) > 3
        Parameter = zeros(4,4); % center frequency; amplitude; bandwidth; skirt width
        c = 0;
        for i = 2:length(M)-1
            if c <= 4
                if yM(i) > yM(i-1) && yM(i) > yM(i+1)
                    c = c + 1;
                    Parameter(1,c) = M(i);
%                     Parameter(1,c) = M(i)*1.2; % pitch up
%                     Parameter(1,c) = M(i)/1.2; % pitch down

%                     Parameter(2,c) = log10(norm(S_linear(M(i-1):M(i+1)))/sqrt(M(i+1)-M(i-1)+1)+1);
                    Parameter(2,c) = norm(S_linear(M(i-1):M(i+1)))/sqrt(M(i+1)-M(i-1)+1);
%                     Parameter(2,c) = log10(S_linear(M(i))+1);

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

                    Parameter(3,c) = (6*(bw_sample(2) - bw_sample(1))/(bw(2) - bw(1)) + 6*(bw_sample(3) - bw_sample(2))/(bw(2) - bw(3)))/2;
                    Parameter(4,c) = max(M(i+1) - M(i),M(i) - M(i-1))*2; 
                end
            end
        end

        index = [];
        for i = 1:4
            if isequal(Parameter(:,i), zeros(4,1))
                index = [index i];
            end
        end
        Parameter(:,index) = [];
        [~,col] = size(Parameter);
        
        if col > 0
            for i = 1:min(col,4)
                beta = Parameter(4,i);
                t = pi/beta;
                alpha = min(Parameter(3,i),c1)*pi;
                f = Parameter(1,i); % *ratio(i);
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
                
                if Parameter(3,i) > c2 && Parameter(3,i) <= c3
                    s = conv(s,rand(1,length(s)+1)*2-1);
                elseif Parameter(3,i) > c3
                    s = conv(s,rand(1,length(s)+1)*2-1);
                    s = conv(s,rand(1,length(s)+1)*2-1);
                end
                
                s = s/max(abs(s))*Parameter(2,i);
                buffer(1:length(s)) = buffer(1:length(s))+s;
            end
        end
    end
    
    if ~isempty(Parameter)
        f0 = Parameter(1,1);
        sst = sst + round(fs/f0);
        number = round(fs/f0);
        number_pitch = round(fs/(f0));
        
        output = buffer(1:number);
        output = output/25;

        adw(output');

        buffer(1:end-number_pitch) = buffer(number_pitch+1:end);
        buffer(end-number_pitch+1:end) = zeros(1,number_pitch);

        sample_buffer(1:end-number) = sample_buffer(number+1:end);
        sample_buffer(end-number+1:end) = zeros(1,number);
        pos = pos - number;
    else
        sst = sst + frame;
        number = frame;
        
%         output = buffer(1:number);
%         output = output/25;
%         adw(output');
%         buffer(1:end-number) = buffer(number+1:end);
%         buffer(end-number+1:end) = zeros(1,number);

        sample_buffer(1:end-number) = sample_buffer(number+1:end);
        sample_buffer(end-number+1:end) = zeros(1,number);
        pos = pos - number;
    end
end
release(afr); 
release(adw);
% % release(deviceReader)