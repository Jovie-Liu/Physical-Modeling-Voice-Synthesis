function Parameter = FOF_process(sample,fs,f_low,f_high,resolt)
    
    S = 20*log10(abs(fft(sample,fs))+0.00001);
    S_0db = S - max(S(f_low+1:f_high));
    S_linear = abs(fft(sample,fs));

    order = round(resolt);
    num = round((f_high-f_low)/order);
    points = zeros(1,2*num);
    for i = 1:num
        piece = S_0db((i-1)*order+1+f_low:i*order+f_low);
        points(i) = find(piece == max(piece),1) + (i-1)*order+f_low;
        points(i+num) = find(piece == min(piece),1) + (i-1)*order+f_low;
    end
    points = unique(points);
    y_p = S_0db(points);

    o = 5;
    ma = [];
    st = 1;
    for i = 1:500
        cal = abs(y_p(st) - y_p);
        ind = 1;
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
            break;
        end
    end

    order = round(resolt)+200;
    num = round((f_high-f_low)/order);
    points2 = zeros(1,num);
    for i = 1: num
        piece = S_0db((i-1)*order+1+f_low:i*order+f_low);
        points2(i) = find(piece == max(piece),1) + (i-1)*order+f_low;
    end
    points2 = [1 points2];

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
                Parameter(2,c) = norm(S_linear(M(i-1):M(i+1)))/sqrt(M(i+1)-M(i-1)+1);


                st = reshape(find(pt == M(i-1)),[1,1]);
                ed = reshape(find(pt == M(i+1)),[1,1]);
                pt_band = pt(st:ed);
                maxima = reshape(find(pt_band == M(i)),[1,1]);
                peak = S_0db(pt_band(maxima));
                y_band = S_0db(pt_band);
                
                
                fr = find(y_band(1:maxima) < peak - 6,1,'last');
                if isempty(fr)
                    fr = 1;
                end
                
                nt = find(y_band(maxima:end) < peak - 6,1)+maxima-1;
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
        [~,col] = size(Parameter);

        index = [];
        for i = 1:col
            if isequal(Parameter(:,i), zeros(4,1))
                index = [index i];
            end
        end
        Parameter(:,index) = [];
    end
end