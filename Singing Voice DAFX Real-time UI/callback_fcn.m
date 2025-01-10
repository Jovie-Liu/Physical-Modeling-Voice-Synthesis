function callback_fcn(obj,event,para_store,p,set_speed_ratio,f1,f2,f3,f4,f5,a1,a2,a3,a4,a5,...
                        t1,t2,t3,t4,t5,autotune,set_freq,deviation,freqValueDisplay,midiValueDisplay,delay_stepsize,...
                        pitch_ratio,amplitude_ratio,vibrato_width,vibrato_rate,beat_ratio,x,f_speed_ratio,ct,x_f1,x_f2,x_f3,x_f4,x_f5,x2,...
                        x_a1,x_a2,x_a3,x_a4,x_a5,x_t1,x_t2,x_t3,x_t4,x_t5,x_auto,store_auto,x_dev,...
                        x_beat,x_vwidth,x_vrate,c1,c2,c3,sst,fs,y_fof,frame,len)
org_sst = sst.value;
while sst.value < org_sst + frame
    if p.value > len
        break
    end
    Parameter = para_store{p.value};
    p.value = p.value + 1;
    [~,col] = size(Parameter);
    drawnow limitrate
    
    if col > 0
        f0 = Parameter(1,1);
        
        if f_speed_ratio.value < set_speed_ratio.value
            f_speed_ratio.value = f_speed_ratio.value + set_speed_ratio.step;
        elseif f_speed_ratio.value > set_speed_ratio.value
            f_speed_ratio.value = f_speed_ratio.value - set_speed_ratio.step;
        end
        step = round(fs/(f0*f_speed_ratio.value));
        
%         freq_ratio_array = [f1.value,f2.value,f3.value,f4.value,f5.value];
%         amp_ratio_array = [a1.value,a2.value,a3.value,a4.value,a5.value];
%         time_delay_array = [t1.value,t2.value,t3.value,t4.value,t5.value];
        
%         x_array = [x_f1.value,x_f2.value,x_f3.value,x_f4.value,x_f5.value];
%         x2_array = [x_a1.value,x_a2.value,x_a3.value,x_a4.value,x_a5.value];
%         x3_array = [x_t1.value,x_t2.value,x_t3.value,x_t4.value,x_t5.value];
        
        if ~isnan(autotune.value)
            if store_auto.value ~= set_freq.value
                autotune.step = 0;
            end
            
            if autotune.step == 0 % set frequency
                
                if x_auto.value < set_freq.value
                    x_auto.value = x_auto.value + set_freq.step;
                elseif x_auto.value > set_freq.value
                    x_auto.value = x_auto.value - set_freq.step;
                end
                store_auto.value = set_freq.value;
                freqq = x_auto.value;
            elseif autotune.step == 1 % midi
                freqq = autotune.value;
                x_auto.value = autotune.value;
                store_auto.value = set_freq.value;
            end
                
            if f0 > 800
                r = 1;
            else
                
                if x_dev.value < deviation.value
                    x_dev.value = x_dev.value + deviation.step;
                elseif x_dev.value > deviation.value
                    x_dev.value = x_dev.value - deviation.step;
                end
                dev = x_dev.value;
                if freqq > f0
                    r = min(freqq,f0 + dev)/f0;
                else
                    r = max(freqq,f0 - dev)/f0;
                end
            end  
            
        else
            r = 1;
                
        end
        
        
        ct.value = ct.value + 1;
        if mod(ct.value,50) == 0 && f0 < 800
            freq_org = f0;
            midi_org = log2(freq_org/440)*12+69;
            set(freqValueDisplay,'String',num2str(freq_org));
            set(midiValueDisplay,'String',num2str(midi_org));
        end
   
        if x.value < pitch_ratio.value
            x.value = x.value + pitch_ratio.step;
        elseif x.value > pitch_ratio.value
            x.value = x.value - pitch_ratio.step;
        end
        
        if x2.value < amplitude_ratio.value
            x2.value = x2.value + amplitude_ratio.step;
        elseif x2.value > amplitude_ratio.value
            x2.value = x2.value - amplitude_ratio.step;
        end
    
        for i = 1:col
            
            beta = Parameter(4,i);
            t = pi/beta;
            alpha = min(Parameter(3,i),c1)*pi;
            
            f = Parameter(1,i)*r*x.value;
            switch i
                case 1 
                    if x_f1.value < f1.value
                        x_f1.value = x_f1.value + pitch_ratio.step;
                    elseif x_f1.value > f1.value
                        x_f1.value = x_f1.value - pitch_ratio.step;
                    end
                    f = f*x_f1.value; % *ratio(i);
                case 2
                    if x_f2.value < f2.value
                        x_f2.value = x_f2.value + pitch_ratio.step;
                    elseif x_f2.value > f2.value
                        x_f2.value = x_f2.value - pitch_ratio.step;
                    end
                    f = f*x_f2.value; % *ratio(i);
                case 3
                    if x_f3.value < f3.value
                        x_f3.value = x_f3.value + pitch_ratio.step;
                    elseif x_f3.value > f3.value
                        x_f3.value = x_f3.value - pitch_ratio.step;
                    end
                    f = f*x_f3.value; % *ratio(i);
                case 4
                    if x_f4.value < f4.value
                        x_f4.value = x_f4.value + pitch_ratio.step;
                    elseif x_f4.value > f4.value
                        x_f4.value = x_f4.value - pitch_ratio.step;
                    end
                    f = f*x_f4.value; % *ratio(i);
                case 5
                    if x_f5.value < f5.value
                        x_f5.value = x_f5.value + pitch_ratio.step;
                    elseif x_f5.value > f5.value
                        x_f5.value = x_f5.value - pitch_ratio.step;
                    end
                    f = f*x_f5.value; % *ratio(i);
            end
            
            omega = 2*pi*f;
            kn1 = [0:fs*t-1]/fs;
            
            if x_vwidth.value < vibrato_width.value
                x_vwidth.value = x_vwidth.value + vibrato_width.step;
            elseif x_vwidth.value > vibrato_width.value
                x_vwidth.value = x_vwidth.value - vibrato_width.step;
            end
            
            if x_vrate.value < vibrato_rate.value
                x_vrate.value = x_vrate.value + vibrato_rate.step;
            elseif x_vrate.value > vibrato_rate.value
                x_vrate.value = x_vrate.value - vibrato_rate.step;
            end
            
            s1 = (1 - cos(beta*kn1))/2.*exp(-alpha*kn1).*sin(omega*kn1+x_vwidth.value*f0*cos(2*pi*x_vrate.value*kn1));
            ts = ceil(-log(0.003)/(alpha*t));
            kn2 = [t*fs:t*fs*ts]/fs;
            s2 = exp(-alpha*kn2).*sin(omega*kn2 + x_vwidth.value*f0*cos(2*pi*x_vrate.value*kn2));
            s = [s1 s2];

            if Parameter(3,i) > c2 && Parameter(3,i) <= c3
                s = conv(s,rand(1,length(s)+1)*2-1);
            elseif Parameter(3,i) > c3
                s = conv(s,rand(1,length(s)+1)*2-1);
                s = conv(s,rand(1,length(s)+1)*2-1);
            end

            amp = Parameter(2,i)*x2.value;
            
            switch i
                case 1 
                    if x_a1.value < a1.value
                        x_a1.value = x_a1.value + amplitude_ratio.step;
                    elseif x_a1.value > a1.value
                        x_a1.value = x_a1.value - amplitude_ratio.step;
                    end
                    amp = amp*x_a1.value; % *ratio(i);
                case 2
                    if x_a2.value < a2.value
                        x_a2.value = x_a2.value + amplitude_ratio.step;
                    elseif x_a2.value > a2.value
                        x_a2.value = x_a2.value - amplitude_ratio.step;
                    end
                    amp = amp*x_a2.value;
                case 3
                    if x_a3.value < a3.value
                        x_a3.value = x_a3.value + amplitude_ratio.step;
                    elseif x_a3.value > a3.value
                        x_a3.value = x_a3.value - amplitude_ratio.step;
                    end
                    amp = amp*x_a3.value;
                case 4
                    if x_a4.value < a4.value
                        x_a4.value = x_a4.value + amplitude_ratio.step;
                    elseif x_a4.value > a4.value
                        x_a4.value = x_a4.value - amplitude_ratio.step;
                    end
                    amp = amp*x_a4.value;
                case 5
                    if x_a5.value < a5.value
                        x_a5.value = x_a5.value + amplitude_ratio.step;
                    elseif x_a5.value > a5.value
                        x_a5.value = x_a5.value - amplitude_ratio.step;
                    end
                    amp = amp*x_a5.value;
            end

            s = s/max(abs(s))*amp;
            
            % add dispersion
            dl = 0;
            
            switch i
                case 1 
                    if x_t1.value < t1.value
                        x_t1.value = x_t1.value + delay_stepsize.value;
                    elseif x_t1.value > t1.value
                        x_t1.value = x_t1.value - delay_stepsize.value;
                    end
                    dl = round(fs*x_t1.value);
                case 2
                    if x_t2.value < t2.value
                        x_t2.value = x_t2.value + delay_stepsize.value;
                    elseif x_t2.value > t2.value
                        x_t2.value = x_t2.value - delay_stepsize.value;
                    end
                    dl = round(fs*x_t2.value);
                case 3
                    if x_t3.value < t3.value
                        x_t3.value = x_t3.value + delay_stepsize.value;
                    elseif x_t3.value > t3.value
                        x_t3.value = x_t3.value - delay_stepsize.value;
                    end
                    dl = round(fs*x_t3.value);
                case 4
                    if x_t4.value < t4.value
                        x_t4.value = x_t4.value + delay_stepsize.value;
                    elseif x_t4.value > t4.value
                        x_t4.value = x_t4.value - delay_stepsize.value;
                    end
                    dl = round(fs*x_t4.value);
                case 5
                    if x_t5.value < t5.value
                        x_t5.value = x_t5.value + delay_stepsize.value;
                    elseif x_t5.value > t5.value
                        x_t5.value = x_t5.value - delay_stepsize.value;
                    end
                    dl = round(fs*x_t5.value);
            end
            
            start = sst.value + dl;
            y_fof.value(start:start+length(s)-1) = y_fof.value(start:start+length(s)-1)+s; 
        end
        sst.value = sst.value + step;
    else
        sst.value = sst.value + 1024;
    end
    
    if x_beat.value < beat_ratio.value
        x_beat.value = x_beat.value + beat_ratio.step;
    elseif x_beat.value > beat_ratio.value
        x_beat.value = x_beat.value - beat_ratio.step;
    end
end
% disp([event.Type ' executed '...
%     datestr(event.Data.time,'dd-mmm-yyyy HH:MM:SS.FFF')])