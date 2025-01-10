function [pitch_ratio,amplitude_ratio,beat_ratio,vibrato_width,vibrato_rate,freq_lowEnd,...
    freq_highEnd,resolution,f1,f2,f3,f4,f5,a1,a2,a3,a4,a5,t1,t2,t3,t4,t5] = main_panel(p1)

pitch_ratio = parameterRef;
pitch_ratio.name = "Pitch ratio";
pitch_ratio.value = 1;
pitch_ratio.step = 0.0001;
parameterTuning(p1,pitch_ratio,0.5,2,[3,415,60,23],[65,415,300,23],[370,415,40,23],[415,415,50,23],[479,415,60,23],[540,415,48,23])

amplitude_ratio = parameterRef;
amplitude_ratio.name = "Ampltd ratio";
amplitude_ratio.value = 1;
amplitude_ratio.step = 0.001;
parameterTuning_linear(p1,amplitude_ratio,0,2,[3,385,60,23],[65,385,300,23],[370,385,40,23],[415,385,50,23],[479,385,60,23],[540,385,48,23])

beat_ratio = parameterRef;
beat_ratio.name = "Beat ratio";
beat_ratio.value = 0;
beat_ratio.step = 0.01;
parameterTuning_linear(p1,beat_ratio,0,2,[3,355,60,23],[65,355,300,23],[370,355,40,23],[415,355,50,23],[479,355,60,23],[540,355,48,23])

vibrato_width = parameterRef;
vibrato_width.name = "Vibrato wid";
vibrato_width.value = 0;
vibrato_width.step = 0.001;
parameterTuning_linear(p1,vibrato_width,0,0.5,[3,325,60,23],[65,325,300,23],[370,325,40,23],[415,325,50,23],[479,325,60,23],[540,325,48,23])

vibrato_rate = parameterRef;
vibrato_rate.name = "Vibrato rate";
vibrato_rate.value = 0;
vibrato_rate.step = 0.1;
parameterTuning_linear(p1,vibrato_rate,0,200,[3,295,60,23],[65,295,300,23],[370,295,40,23],[415,295,50,23],[479,295,60,23],[540,295,48,23])

freq_lowEnd = parameterRef;
freq_lowEnd.name = "Freq lowEnd";
freq_lowEnd.value = 0;
freq_lowEnd.step = 10;
parameterTuning_linear(p1,freq_lowEnd,0,1000,[3,265,60,23],[65,265,300,23],[370,265,40,23],[415,265,50,23],[479,265,60,23],[540,265,48,23])

freq_highEnd = parameterRef;
freq_highEnd.name = "Freq highEnd";
freq_highEnd.value = 10000;
freq_highEnd.step = 100;
parameterTuning_linear(p1,freq_highEnd,1000,11000,[3,235,64,23],[65,235,300,23],[370,235,40,23],[415,235,50,23],[479,235,60,23],[540,235,48,23])

resolution = parameterRef;
resolution.name = "Resolution";
resolution.value = 350;
resolution.step = 10;
parameterTuning_linear(p1,resolution,100,1000,[3,205,60,23],[65,205,300,23],[370,205,40,23],[415,205,50,23],[479,205,60,23],[540,205,48,23])

f1 = parameterRef;
f1.name = "f1";
f1.value = 1;
parameterTuning_sub(p1,f1,0.5,2,[3,170,20,23],[25,170,120,23],[150,179,40,14],[150,160,40,18])

f2 = parameterRef;
f2.name = "f2";
f2.value = 1;
parameterTuning_sub(p1,f2,0.5,2,[3,130,20,23],[25,130,120,23],[150,139,40,14],[150,120,40,18])

f3 = parameterRef;
f3.name = "f3";
f3.value = 1;
parameterTuning_sub(p1,f3,0.5,2,[3,90,20,23],[25,90,120,23],[150,99,40,14],[150,80,40,18])

f4 = parameterRef;
f4.name = "f4";
f4.value = 1;
parameterTuning_sub(p1,f4,0.5,2,[3,50,20,23],[25,50,120,23],[150,59,40,14],[150,40,40,18])

f5 = parameterRef;
f5.name = "f5";
f5.value = 1;
parameterTuning_sub(p1,f5,0.5,2,[3,10,20,23],[25,10,120,23],[150,19,40,14],[150,0,40,18])

a1 = parameterRef;
a1.name = "a1";
a1.value = 1;
parameterTuning_linear_sub(p1,a1,0,2,[203,170,20,23],[225,170,120,23],[350,179,40,14],[350,160,40,18])

a2 = parameterRef;
a2.name = "a2";
a2.value = 1;
parameterTuning_linear_sub(p1,a2,0,2,[203,130,20,23],[225,130,120,23],[350,139,40,14],[350,120,40,18])

a3 = parameterRef;
a3.name = "a3";
a3.value = 1;
parameterTuning_linear_sub(p1,a3,0,2,[203,90,20,23],[225,90,120,23],[350,99,40,14],[350,80,40,18])

a4 = parameterRef;
a4.name = "a4";
a4.value = 1;
parameterTuning_linear_sub(p1,a4,0,2,[203,50,20,23],[225,50,120,23],[350,59,40,14],[350,40,40,18])

a5 = parameterRef;
a5.name = "a5";
a5.value = 1;
parameterTuning_linear_sub(p1,a5,0,2,[203,10,20,23],[225,10,120,23],[350,19,40,14],[350,0,40,18])

t1 = parameterRef;
t1.name = "t1";
t1.value = 0;
parameterTuning_linear_sub(p1,t1,0,2,[403,170,20,23],[425,170,120,23],[549,179,40,14],[549,160,40,18])

t2 = parameterRef;
t2.name = "t2";
t2.value = 0;
parameterTuning_linear_sub(p1,t2,0,2,[403,130,20,23],[425,130,120,23],[549,139,40,14],[549,120,40,18])

t3 = parameterRef;
t3.name = "t3";
t3.value = 0;
parameterTuning_linear_sub(p1,t3,0,2,[403,90,20,23],[425,90,120,23],[549,99,40,14],[549,80,40,18])

t4 = parameterRef;
t4.name = "t4";
t4.value = 0;
parameterTuning_linear_sub(p1,t4,0,2,[403,50,20,23],[425,50,120,23],[549,59,40,14],[549,40,40,18])

t5 = parameterRef;
t5.name = "t5";
t5.value = 0;
parameterTuning_linear_sub(p1,t5,0,2,[403,10,20,23],[425,10,120,23],[549,19,40,14],[549,0,40,18])
end