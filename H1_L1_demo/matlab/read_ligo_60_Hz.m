% Compute maximum ampitude of 60*n hz tones in H1 data

%N point FFT scaling = N/2
% Let signal x[n]= I * cos(2*pi*f1.*n/fs) + Q * sin(2*pi*f1.*n/fs)
% Discrete Fourier Transform is given by X[k]= [sum over n=0:N-1] x[n]*exp( -j*2*pi*k*n/N)
% Let f1/fs = k1/N
% X[k]= [sum over n=0:N-1] [I * cos(2*pi*k1.*n/N) + Q * sin(2*pi*k1.*n/N)] * [cos(2*pi*k*n/N) - j* sin(2*pi*k*n/N) ]
% X[k1]= (I/2) * N -j*(Q/2)*N = [I-jQ] * (N/2)

%Set test_FFT_amplitude=1 if you want to test the FFT scaling of 60*n Hz tones
test_FFT_amplitude=0;


inFile=fopen('H-H1_LOSC_4_V1-1126259446-32.txt','r');

str_1= fgets(inFile);
str_1= fgets(inFile);
str_1= fgets(inFile);

N1=65536*2; 
for count=1:N1,
  dataArray_H(count)=fscanf(inFile,'%f',1);
end

fclose(inFile);

timeArray=[0:1:N1-1]/4096;
timeArray_1=timeArray-16; %32 seconds data; -16 to +16 seconds from tevent

Fs=4096; %sampling frequency

%FFT of H1 data
fft_sig=fft(dataArray_H);

%60*n hz tones; include frequencies +- 1 Hz surrounding each tone; 
f1=60;f2=120; f3=180;
freqArray=[f1-1:f1+1 f2-1:f2+1 f3-1:f3+1];
%ampArray = abs(fft_sig(freqArray*32+1) )/(N1/2); %frequency resolution=Fs/N1 = 1/32 Hz

%FFT scaling = N/2
ampArrayCos = real(fft_sig(freqArray*32+1) )/(N1/2); 
ampArraySin = -imag(fft_sig(freqArray*32+1) )/(N1/2); 

% reconstruct 60*n hz tones from H1 data
sig_2=zeros(1,N1);
for count=1:length(ampArrayCos),
%sig_2= sig_2 + ampArray(count).* cos(2*pi.*[0:1:N1-1]*freqArray(count)/Fs) ;
sig_2= sig_2 + ampArrayCos(count).* cos(2*pi.*[0:1:N1-1]*freqArray(count)/Fs) + ampArraySin(count).* sin(2*pi.*[0:1:N1-1]*freqArray(count)/Fs) ;
end

%compute power of 60*n hz tones and power ratio
disp('Maximum value of 60*n Hz tones')
max_tones= max(abs(sig_2))


disp('Maximum value of H1 signal')
max_gw=1e-21

disp('Maximum value of Raw H1 unprocessed signal')
max_total=max(abs(dataArray_H))


%max_total/max_gw
%max_total/max_tones

%Plot
figure(1)
subplot(2,1,1)
hold off
plot(timeArray_1(65536+(-800:800)), dataArray_H(65536+(-800:800)))
grid on
title('Raw H1 signal')
text(0, -3e-19, ' dominated by 10-15 Hz component') %gtext(' dominated by 10-15 Hz component') 

%figure(2)
subplot(2,1,2)
hold off
plot(timeArray_1(65536+(-800:800)), sig_2(65536+(-800:800)))
grid on
title('Extracted 60*n hz signal from Raw H1 signal')
xlabel('Time in seconds')


%Test FFT Tone amplitudes
if test_FFT_amplitude==1,
    
%60*n hz tones; include frequencies +- 1 Hz surrounding each tone; 
f1=60;f2=120; f3=180;
freqArray=[f1 f2 f3];

% Use any desired input amplitudes
ampArrayCosIn=[1 -2 3];
ampArraySinIn=[-4 5 -6];

sig_2=zeros(1,N1);
for count=1:length(ampArrayCosIn),
  sig_2= sig_2 +   ampArrayCosIn(count).* cos(2*pi.*[0:1:N1-1]*freqArray(count)/Fs) + ampArraySinIn(count).* sin(2*pi.*[0:1:N1-1]*freqArray(count)/Fs) ;
end
dataArray_H= sig_2;

%FFT of H1 data
fft_sig=fft(dataArray_H);

%FFT scaling = N/2
ampArrayCos = real(fft_sig(freqArray*32+1) )/(N1/2); 
ampArraySin = -imag(fft_sig(freqArray*32+1) )/(N1/2); 

% reconstruct 60*n hz tones from H1 data
sig_2=zeros(1,N1);
for count=1:length(ampArrayCos),
  sig_2= sig_2 + ampArrayCos(count).* cos(2*pi.*[0:1:N1-1]*freqArray(count)/Fs) + ampArraySin(count).* sin(2*pi.*[0:1:N1-1]*freqArray(count)/Fs) ;
end

disp('Input amplitudes of 60 Hz, 120Hz, 180Hz Cosine and Sine Tones.....')
ampArrayCosIn
ampArraySinIn

disp('Output amplitudes of 60 Hz, 120Hz, 180Hz Cosine and Sine Tones.....')
ampArrayCos
ampArraySin

end

