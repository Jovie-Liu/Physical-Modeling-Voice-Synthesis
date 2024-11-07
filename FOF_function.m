function y_fof = FOF_function(y,fs)
    ny = length(y);
    nh = 1024;
    input = [y zeros(1,nh)];
    y_fof = zeros(1,length(input)+fs);
    n = length(input);
    sst = 1;
    alpha = 15*pi;
    while sst < ny
        sample = input(sst:sst+nh-1);
        S = 20*log10(abs(fft(sample,fs))+0.00001);
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

        M =  [1 points(ma)];
        yM = S_0db(M);
        Parameter = [];
        if length(M) > 3
            Parameter = zeros(3,ceil((length(M)-2)/2)); % center frequency; amplitude; skirt width
            c = 0;
            for i = 2:length(M)-1
                if yM(i) > yM(i-1) && yM(i) > yM(i+1)
                    if S_linear(M(i)) > 0.003
                        c = c + 1;
                        Parameter(1,c) = M(i);
                        Parameter(2,c) = S_linear(M(i));
                        Parameter(3,c) = M(i+1) - M(i-1); % max(M(i+1) - M(i),M(i) - M(i-1))*2
                    end
                end
            end
            [~,col] = size(Parameter);
    
            index = [];
            for i = 1:col
                if isequal(Parameter(:,i), zeros(3,1))
                    index = [index i];
                end
            end
            Parameter(:,index) = [];
            [~,col] = size(Parameter);
            
            if col > 0
                for i = 1:col
                    beta = Parameter(3,i);
                    t = pi/beta;
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
                    y_fof(sst:sst+length(s)-1) = y_fof(sst:sst+length(s)-1)+s;
                end
            end
        end
        if ~isempty(Parameter)
            f0 = Parameter(1,1);
            sst = sst + round(fs/f0);
        else
            sst = sst + nh;
        end
    end
    y_fof = y_fof/max(y_fof);
    y_fof = y_fof(1:ny);
end