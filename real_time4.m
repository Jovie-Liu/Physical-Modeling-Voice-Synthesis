fs = 44100;
frame = 128;
buffer_size = 8;

afr = dsp.AudioFileReader('didadi_vocal.wav','ReadRange',[78*fs 86.5*fs],'SamplesPerFrame',frame);
% afr = dsp.AudioFileReader('speech.mp3','ReadRange',[56*fs 62.5*fs],'SamplesPerFrame',frame);
% afr = dsp.AudioFileReader('verse1_vocal.wav','ReadRange',[0.5*fs 5.5*fs],'SamplesPerFrame',frame);
% afr = dsp.AudioFileReader('Tokyo_vocal.wav','ReadRange',[29*fs 34*fs],'SamplesPerFrame',frame);

adw = audioDeviceWriter('SampleRate', afr.SampleRate);

c1 = 80; %80
c2 = 350; %300
c3 = 1000;
Phi = 0;

buffer = zeros(1,2*fs);
sst = 1;

sample_buffer = zeros(1,buffer_size*frame);
for i = 1:buffer_size
    audio = afr();
    audio = (audio(:,1)+audio(:,2))/2;
    sample_buffer(1+(i-1)*frame:i*frame) = audio';
end

while ~isDone(afr)
    counter = 0;
    while sst < frame && counter < 3
        counter = counter + 1;
        sample = sample_buffer(sst:end);
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
            Parameter = zeros(4,ceil((length(M)-2)/2)); % center frequency; amplitude; bandwidth; skirt width
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
                f0 = Parameter(1,1);
                for i = 1:col
                    beta = Parameter(4,i);
                    t = pi/beta;
                    alpha = min(Parameter(3,i),c1)*pi;
    %                 alpha = 15*pi;
                    f = Parameter(1,i); % *ratio(i);
                    omega = 2*pi*f;
                    kn1 = [0:fs*t-1]/fs;
                    vibrato_width = 0.01*f;
                    vibrato_rate = randi(10);
    %                 vibrato_width = 0;
    %                 vibrato_rate = 0;
                    s1 = (1 - cos(beta*kn1))/2.*exp(-alpha*kn1).*sin(omega*kn1+vibrato_width*cos(2*pi*vibrato_rate*kn1)+Phi);
                    for ts = 1:3000
                        st = exp(-alpha*t*ts);
                        if st < 0.003
                            break
                        end
                    end
                    kn2 = [t*fs:t*fs*ts]/fs;
                    s2 = exp(-alpha*kn2).*sin(omega*kn2 + vibrato_width*cos(2*pi*vibrato_rate*kn2)+Phi);
                    s = [s1 s2];

                    if Parameter(3,i) > c2 && Parameter(3,i) <= c3
                        s = conv(s,rand(1,length(s)+1)*2-1);
                    elseif Parameter(3,i) > c3
                        s = conv(s,rand(1,length(s)+1)*2-1);
                        s = conv(s,rand(1,length(s)+1)*2-1);
                    end

                    s = s/max(abs(s))*Parameter(2,i);
                    buffer(sst:sst+length(s)-1) = buffer(sst:sst+length(s)-1)+s;
                    
                end
                step = round(fs/f0);
                sst = sst + step;
                if Phi > 1000
                    Phi = mod((step/fs*f0 - floor(step/fs*f0))*2*pi + Phi,2*pi);
                else
                    Phi = (step/fs*f0 - floor(step/fs*f0))*2*pi + Phi;
                end
                
            else
                sst = sst + frame;
            end
        else
            sst = sst + frame;
        end
        
    end

    sst = sst - frame;
    if sst < 0
        sst = sst + frame; % recover sst
        while sst < frame
            f0 = Parameter(1,1);
            buffer(sst:sst+length(s)-1) = buffer(sst:sst+length(s)-1)+s;
            sst = sst + round(fs/f0) + randi(20)-10;
        end
        sst = sst - frame;
    end
    if sst == 0
        sst = 1;
    end
    
    output = buffer(1:frame);
    output = output/25;

    adw(output');

    buffer(1:end-frame) = buffer(frame+1:end);
    buffer(end-frame+1:end) = zeros(1,frame);
    
    audio = afr(); %[1024,2]
    audio = (audio(:,1)+audio(:,2))/2;
%     audio = audio(:,1);

    sample_buffer(1:(buffer_size-1)*frame) = sample_buffer(1+frame:end);
    sample_buffer((buffer_size-1)*frame+1:end) = audio';
    
end
release(afr); 
release(adw);