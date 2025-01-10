function [pitch_ratio,amplitude_ratio,time_delay,vibrato_width,vibrato_rate,freq_lowEnd,...
    freq_highEnd,resolution,f1,f2,f3,f4,f5,a1,a2,a3,a4,a5,t1,t2,t3,t4,t5] = left_panel(p1)

pitch_ratio = parameterRef;
pitch_ratio.name = "pitch ratio";
pitch_ratio.value = 1;
pitch_ratio.step = 0.0001;
parameterTuning_sub(p1,pitch_ratio,0.5,2,[3,415,60,23],[65,415,300,23],[370,415,40,23],[415,415,43.5,23])

amplitude_ratio = parameterRef;
amplitude_ratio.name = "ampltd ratio";
amplitude_ratio.value = 1;
amplitude_ratio.step = 0.001;
parameterTuning_linear_sub(p1,amplitude_ratio,0,2,[3,385,60,23],[65,385,300,23],[370,385,40,23],[415,385,43.5,23])

time_delay = parameterRef;
time_delay.name = "time delay";
time_delay.value = 0;
time_delay.step = 0.001;
parameterTuning_linear_sub(p1,time_delay,0,2,[3,355,60,23],[65,355,300,23],[370,355,40,23],[415,355,43.5,23])

vibrato_width = parameterRef;
vibrato_width.name = "vibrato wid";
vibrato_width.value = 1;
vibrato_width.step = 0.001;
parameterTuning_linear_sub(p1,vibrato_width,0,2,[3,325,60,23],[65,325,300,23],[370,325,40,23],[415,325,43.5,23])

vibrato_rate = parameterRef;
vibrato_rate.name = "vibrato rate";
vibrato_rate.value = 1;
vibrato_rate.step = 0.001;
parameterTuning_linear_sub(p1,vibrato_rate,0,2,[3,295,60,23],[65,295,300,23],[370,295,40,23],[415,295,43.5,23])

freq_lowEnd = parameterRef;
freq_lowEnd.name = "freq lowEnd";
freq_lowEnd.value = 0;
freq_lowEnd.step = 0.001;
parameterTuning_linear_sub(p1,freq_lowEnd,0,2,[3,265,60,23],[65,265,300,23],[370,265,40,23],[415,265,43.5,23])

freq_highEnd = parameterRef;
freq_highEnd.name = "freq highEnd";
freq_highEnd.value = 0;
freq_highEnd.step = 0.001;
parameterTuning_linear_sub(p1,freq_highEnd,0,2,[3,235,60,23],[65,235,300,23],[370,235,40,23],[415,235,43.5,23])

resolution = parameterRef;
resolution.name = "resolution";
resolution.value = 0;
resolution.step = 0.001;
parameterTuning_linear_sub(p1,resolution,0,2,[3,205,60,23],[65,205,300,23],[370,205,40,23],[415,205,43.5,23])

f1 = parameterRef;
f1.name = "f1";
f1.value = 1;
parameterTuning_sub(p1,f1,0.5,2,[3,170,15,23],[20,170,90,23],[112,179,38,14],[112,160,38,18])

f2 = parameterRef;
f2.name = "f2";
f2.value = 1;
parameterTuning_sub(p1,f2,0.5,2,[3,130,15,23],[20,130,90,23],[112,139,38,14],[112,120,38,18])

f3 = parameterRef;
f3.name = "f3";
f3.value = 1;
parameterTuning_sub(p1,f3,0.5,2,[3,90,15,23],[20,90,90,23],[112,99,38,14],[112,80,38,18])

f4 = parameterRef;
f4.name = "f4";
f4.value = 1;
parameterTuning_sub(p1,f4,0.5,2,[3,50,15,23],[20,50,90,23],[112,59,38,14],[112,40,38,18])

f5 = parameterRef;
f5.name = "f5";
f5.value = 1;
parameterTuning_sub(p1,f5,0.5,2,[3,10,15,23],[20,10,90,23],[112,19,38,14],[112,0,38,18])

a1 = parameterRef;
a1.name = "a1";
a1.value = 1;
parameterTuning_linear_sub(p1,a1,0,2,[157,170,15,23],[174,170,90,23],[266,179,38,14],[266,160,38,18])

a2 = parameterRef;
a2.name = "a2";
a2.value = 1;
parameterTuning_linear_sub(p1,a2,0,2,[157,130,15,23],[174,130,90,23],[266,139,38,14],[266,120,38,18])

a3 = parameterRef;
a3.name = "a3";
a3.value = 1;
parameterTuning_linear_sub(p1,a3,0,2,[157,90,15,23],[174,90,90,23],[266,99,38,14],[266,80,38,18])

a4 = parameterRef;
a4.name = "a4";
a4.value = 1;
parameterTuning_linear_sub(p1,a4,0,2,[157,50,15,23],[174,50,90,23],[266,59,38,14],[266,40,38,18])

a5 = parameterRef;
a5.name = "a5";
a5.value = 1;
parameterTuning_linear_sub(p1,a5,0,2,[157,10,15,23],[174,10,90,23],[266,19,38,14],[266,0,38,18])

t1 = parameterRef;
t1.name = "t1";
t1.value = 0;
parameterTuning_linear_sub(p1,t1,0,2,[311,170,15,23],[328,170,90,23],[420,179,38,14],[420,160,38,18])

t2 = parameterRef;
t2.name = "t2";
t2.value = 0;
parameterTuning_linear_sub(p1,t2,0,2,[311,130,15,23],[328,130,90,23],[420,139,38,14],[420,120,38,18])

t3 = parameterRef;
t3.name = "t3";
t3.value = 0;
parameterTuning_linear_sub(p1,t3,0,2,[311,90,15,23],[328,90,90,23],[420,99,38,14],[420,80,38,18])

t4 = parameterRef;
t4.name = "t4";
t4.value = 0;
parameterTuning_linear_sub(p1,t4,0,2,[311,50,15,23],[328,50,90,23],[420,59,38,14],[420,40,38,18])

t5 = parameterRef;
t5.name = "t5";
t5.value = 0;
parameterTuning_linear_sub(p1,t5,0,2,[311,10,15,23],[328,10,90,23],[420,19,38,14],[420,0,38,18])
end