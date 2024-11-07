function para_store = gen_parameter(audiofile)
[NY,fs] = audioread(audiofile);
y = NY(:,1)';

ny = length(y);
y = y/max(abs(y));
nh = 1024;
input = [zeros(1,nh/4) y zeros(1,nh)];
sst = nh/4 + 1;
window = chebwin(nh);
freq = 11000;

para_store = {};
while sst < ny-nh/4
    sample = input(sst-nh/4:sst+nh*3/4-1).*window';
    S = 20*log10(abs(fft(sample,fs))+0.00001);
    S_0db = S - max(S(1:freq));
    S_linear = abs(fft(sample,fs));
    
%     figure
%     n = fs;
%     fn = [0:1/n:1-1/n]*fs;
%     plot(fn,S_0db)
%     xlim([0 8000])
    
    order = 350;
    num = round(freq/order);
    points = zeros(1,2*num);
    for i = 1:num
        piece = S_0db((i-1)*order+1:i*order);
        points(i) = find(piece == max(piece),1) + (i-1)*order;
        points(i+num) = find(piece == min(piece),1) + (i-1)*order;
    end
    points = unique(points);
    y_p = S_0db(points);
    
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
    
%     hold on
%     scatter(points(ma),S_0db(points(ma)))
    
    order = 550;
    num = round(freq/order);
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

                Parameter(2,c) = norm(S_linear(M(i-1):M(i+1)))/sqrt(M(i+1)-M(i-1)+1);

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
    end
    para_store{end+1} = Parameter;
    if ~isempty(Parameter)
        f0 = Parameter(1,1);
        step = round(fs/f0);
        sst = sst + step;
    else
        sst = sst + nh;
    end
end