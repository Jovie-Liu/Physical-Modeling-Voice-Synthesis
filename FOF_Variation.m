function y_fof_v = FOF_Variation(y,fs,property,note)
    ny = length(y);
    nh = 1024;
    input = [y zeros(1,nh)];
    y_fof_v = zeros(1,length(input)+fs);
    n = length(input);
    sst = 1;
    t_sin = 0;
    set = [62 67 69 71];
    flag = 1;
    pitch_true = 100;
    pitch_change = 100;
    midi_choose = set(1);
    while sst < ny
        sample = input(sst:sst+nh-1);
        S = 20*log10(abs(fft(sample,fs))+eps);
        S_0db = S - max(S);
        S_linear = abs(fft(sample,fs));

        order = 100;
        num = round(8000/order);
        points = zeros(1,num);
        for i = 1: num
            piece = S_0db((i-1)*order+1:i*order);
            points(i) = find(piece == max(piece),1) + (i-1)*order;
        end
        y_p = S_0db(points);

        o = 5;
        ma = [];
        st = 1;
        for i = 1:500
            cal = abs(y_p(st) - y_p);
            for j = st+1:num
                ind = find(cal(st:j) == max(cal(st:j)),1)+st-1;
                if ind + o <= num && ind ~= j && ind == find(cal(st:ind+o) == max(cal(st:ind+o)),1)+st-1
                    ma = [ma ind];
                    st = ind;
                    break
                elseif ind + o > num
                    ma = [ma ind];
                    break
                end
            end
            if ind + o > num || j == num
                break
            end
        end

        order = 300;
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
            Parameter = zeros(4,ceil((length(M)-2)/2)); % center frequency; amplitude. bandwidth; skirt width
            c = 0;
            for i = 2:length(M)-1
                if yM(i) > yM(i-1) && yM(i) > yM(i+1)
                    if S_linear(M(i)) > 0.005
                        c = c + 1;
                        Parameter(1,c) = M(i);
                        Parameter(2,c) = S_linear(M(i));

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
                        Parameter(3,c) = min(Parameter(3,c),-log(0.1)*Parameter(1,1));
                        if Parameter(3,c) < 0
                            error('not positive bandwidth')
                        end
                        Parameter(4,c) = M(i+1) - M(i-1); % max(M(i+1) - M(i),M(i) - M(i-1))*2

                %         vq = interp1(pt_band,S_0db(pt_band),M(i-1):M(i+1),'spline');
                %         hold on
                %         plot(M(i-1):M(i+1),vq)
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
                for i = 1:col
                    beta = Parameter(4,i);
                    t = pi/beta;
%                     alpha = Parameter(3,i)*pi;
                    alpha = min(Parameter(3,i),150)*pi;
                    f = Parameter(1,i);
                    omega = 2*pi*f;
                    Phi = 0;
                    kn1 = [0:fs*t-1]/fs;
                    s1 = (1 - cos(beta*kn1))/2.*exp(-alpha*kn1).*sin(omega*kn1+Phi);
                    for ts = 1:3000
                        st = exp(-alpha*t*ts);
                        if st < 0.005
                            break
                        end
                    end
                    kn2 = [t*fs:t*fs*ts]/fs;
                    s2 = exp(-alpha*kn2).*sin(omega*kn2 + Phi);
                    s = [s1 s2]*Parameter(2,i);
                    y_fof_v(sst:sst+length(s)-1) = y_fof_v(sst:sst+length(s)-1)+s;
                end
            end
        end

    % stable pitch
    %     if ~isempty(Parameter)
    %         f0 = Parameter(1,1);
    %         midi = round(log2(f0/440)*12+69);
    %         f = 440*2^((midi - 69)/12);
    %         sst = sst + round(fs/f);
    %     else
    %         sst = sst + nh/2;
    %     end

    % choose pitch
        switch property
            case 'melody'
                if ~isempty(Parameter)
                    f0 = Parameter(1,1);
                    midi = log2(f0/440)*12+69;
                    if midi > pitch_true - 1 && midi < pitch_true + 1
                        flag = 0;
                        midi_choose = pitch_change;% + midi - pitch_true;
                    else
                        flag = 1;
                    end
                    if flag == 1
                        indicator = 0;
                        for item = 1:length(set)
                            if round(midi) == set(item)
                                indicator = 1;
                                pitch_true = set(item);
                                pitch_change = set(randi(length(set)));
                                midi_choose = pitch_change; %+ midi - pitch_true;
                                break
                            end
                        end
                        if indicator == 0
                            midi_choose = pitch_change; %+ midi - pitch_true;
                        end
                    end
                    f_choose = 440*2^((midi_choose - 69)/12);
                    sst = sst + round(fs/f_choose);
                else
                    f_choose = 440*2^((midi_choose - 69)/12);
                    sst = sst + round(fs/f_choose);
                end
                
            case 'flat'
                midi_choose = note;
                f_choose = 440*2^((midi_choose - 69)/12);
                sst = sst + round(fs/f_choose);
        end
    end
    y_fof_v = y_fof_v/max(y_fof_v);
    y_fof_v = y_fof_v(1:ny);
end