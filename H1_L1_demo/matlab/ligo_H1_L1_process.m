%Read H1/L1 data and process using whiteneing and 20-300Hz filter, as per
%LIGO python scripts
     
inFile=fopen('H-H1_LOSC_4_V1-1126259446-32.txt','r');


str_1= fgets(inFile);
str_1= fgets(inFile);
str_1= fgets(inFile);

N1=65536*2; 
for count=1:N1,
  dataArray_H(count)=fscanf(inFile,'%f',1);
end

fclose(inFile);

timeArray=[0:1:N1-1]/4096 - 65536/4096; %16384;

Fs=4096;




%%PSD
inFile=fopen('H1_Pxx.dat','r');

N1=2048; 
for count=1:N1,
  H1_Pxx(count)=fscanf(inFile,'%f',1);
end

fclose(inFile);



inFile=fopen('H1_psd_2.dat','r');

N1=65536; 
for count=1:N1,
  H1_psd_2(count)=fscanf(inFile,'%f',1);
end

fclose(inFile);


%Synthesize H1_Pxx
Fs = 4096;   
x = dataArray_H;
nfft = 2048*2;
Pxx = abs(fft(x,nfft)).^2/length(x)/Fs;

 
[H1_Pxx_my,freq_my] = pwelch(x, hanning(nfft),0, nfft, Fs);
 


new_freq=[0:2048/length(H1_psd_2):2048];

H1_psd_2_my=sqrt(interp1(freq_my,H1_Pxx_my,new_freq'));



%whitening
array_1=[H1_psd_2_my(1:65536)' H1_psd_2_my(65537:-1:2)'];
fft_sig=fft(dataArray_H);
fft_sig_2=fft_sig./array_1;
sig_2=real(ifft(fft_sig_2));


%BPF 20-300Hz filter
filterOrder=4; 
Wn=[20 300]*2/Fs;
ftype = 'bandpass';
[B,A] = butter(filterOrder,Wn,ftype); %

%Don't use 20-300Hz filter
sig_3_filtered = sig_2; %filtfilt(B,A,sig_2);
%sig_3_filtered = filtfilt(B,A,sig_2);


%scan 20-300Hz filtered version from file.
N1=66445;
N2=68902;

inFile=fopen('H1_filtered.dat','r'); %Ideal processed H1 from Python scripts. GW_tutorial.py

N3=131072-1; 
for count=1:N3,
  H1_filtered(count)=fscanf(inFile,'%f',1);
end

fclose(inFile);

%Plots
figure(11)
hold off
plot(timeArray(N1:N2),sig_3_filtered(N1:N2)*max(H1_filtered(N1:N2))/max(sig_3_filtered(N1:N2)))
hold on
plot(timeArray(N1:N2),H1_filtered(N1:N2),'r-')
grid on
%title('H1 strain filtered. 30-350Hz')
title('Blue: H1 strain without filtering. 0-2048 Hz. Red: H1 strain with filtering. 20-300 Hz')
xlabel('Time in seconds')

sig_0=sig_3_filtered(N1-10000:N1+2048-1-10000);
sig_1=sig_3_filtered(N1:N1+2048-1);

fft_sig_0=fft(sig_0);
fft_sig_1=fft(sig_1);

inFile=fopen('GW150914_4_NR_waveform.txt','r');

N1=2769;
for count=1:N1,
  timeArray(count)=fscanf(inFile,'%f',1);
  dataArray_NR(count)=fscanf(inFile,'%f',1);
end

fclose(inFile);
fft_sig_NR=fft(dataArray_NR(2769-2048+1:2769));

freqArray=[0:2:2048*2-1];

figure(21)
subplot(3,1,1)
hold off
plot(freqArray,abs(fft_sig_0(1:2048)))
grid on
title('FFT of H1 whitened signal during non-GW block')


subplot(3,1,2)
hold off
plot(freqArray,abs(fft_sig_1(1:2048)))
grid on
title('FFT of H1 whitened signal during GW block')

subplot(3,1,3)
hold off
plot(freqArray,abs(fft_sig_NR(1:2048)))
grid on
title('FFT of NR template')
xlabel('Frequency in Hz')


freqArray=[10:2:2048-1];

figure(20)
subplot(3,1,1)
hold off
loglog(freqArray,abs(fft_sig_0(5+1:2048/2)))
grid on
title('FFT of H1 whitened signal during non-GW block')

subplot(3,1,2)
hold off
loglog(freqArray,abs(fft_sig_1(5+1:2048/2)))
grid on
title('FFT of H1 signal  during GW block')

subplot(3,1,3)
hold off
loglog(freqArray,abs(fft_sig_NR(5+1:2048/2)))
grid on
title('FFT of NR template')
xlabel('Frequency in Hz')
