
%Read H1/L1 data and process using , impulse removal, interplation and 50-900Hz filter, as per
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

timeArray=[0:1:N1-1]/16384;
Fs=4096;

%FFT
fft_dataArray_H= fft(dataArray_H);




%Filter
filterOrder=4; 
Wn=[50 900]*2/Fs; K2=10; %K2=100;
[B,A] = cheby2(filterOrder,50,Wn,'bandpass');

dataArray_H_filtered = filtfilt(B,A,dataArray_H);



%Compute power, starting from center


N0=66680-1024*50; 
N1=1024; %4 Hz resolution



 sig=dataArray_H_filtered; %( N0+(count0-1)*N1+(1:N1) );
 f_step=Fs/length(sig); %f_step=Fs/(65536*2);
 fft_sig=fft(sig);
 
%Remove impulses
fft_sig_save=fft_sig;    
toneIndex=[60 120 180 300 303 332 502 505 508]*32;
K2=10; 
fftThrIndex=[];
for count1=1:length(toneIndex),
  fftThrIndex= [fftThrIndex toneIndex(count1)-32*1:toneIndex(count1)+32*1] ; %find( abs(fft_sig)> fftThr);
end
fft_sig( fftThrIndex) = 0;


%Interpolate linearly zero values
count=1;
while count < length(fftThrIndex), %/2,
    
  freqArray(1)= (fftThrIndex(count)-2)*f_step;
  fftRealArray(1)= real(fft_sig(fftThrIndex(count)-1));
  fftImagArray(1)= imag(fft_sig(fftThrIndex(count)-1));
  
  count2=2;
  while count<length(fftThrIndex) && fftThrIndex(count+1)- fftThrIndex(count)==1 ,
      count=count+1;
      count2=count2+1;
  end

  
  freqArray(2)= (fftThrIndex(count))*f_step;
  fftRealArray(2)= real(fft_sig(fftThrIndex(count)+1));
  fftImagArray(2)= imag(fft_sig(fftThrIndex(count)+1));
  
  new_freq=[freqArray(1):(freqArray(2)-freqArray(1))/count2:freqArray(2)];
  newValuesReal= interp1(freqArray,fftRealArray,new_freq');
  newValuesImag= interp1(freqArray,fftImagArray,new_freq');

  
  fft_sig(fftThrIndex(count-count2+2:count) )= newValuesReal(2:end-1)+ i.* newValuesImag(2:end-1);
  fft_sig(length(fft_sig)+2-(fftThrIndex(count-count2+2:count)) )= newValuesReal(2:end-1)- i.* newValuesImag(2:end-1);

  count=count+1;
  
end


sig_mod=real(ifft(fft_sig));


freqArray=[10:1/32:2048-1/32];

fft_sig_mod = fft(sig_mod);

%Plots
figure(2)
subplot(2,1,1)
hold off
loglog(freqArray, abs(fft_dataArray_H(10*32+1:65536)))
xlabel('Frequency in Hz')
title('FFT of H1 Raw signal')
grid on
subplot(2,1,2)
hold off
loglog(freqArray, abs(fft_sig_mod(10*32+1:65536)))
xlabel('Frequency in Hz')
title('FFT of H1 signal with impulsive tones removed, filtered with 50-900Hz')
grid on




inFile=fopen('H1_filtered.dat','r');

N3=131072-1; 
for count=1:N3,
  H1_filtered(count)=fscanf(inFile,'%f',1);
end

fclose(inFile);

N1=66445; 
N2=68902;


figure(1)
hold off
plot(sig_mod(N1:N2))
grid on
hold on
plot(H1_filtered(N1:N2)*max(sig_mod(N1:N2))/max(H1_filtered(N1:N2)),'r-')
title('Red: H1 whitened filtered with 20-300Hz filter. Blue: H1 in 50-900Hz, with impulses removed and interpolated')




sig_0=sig_mod(N1-30000:N1+2048-1-30000);
sig_1=sig_mod(N1:N1+2048-1);

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
    
figure(22)
subplot(3,1,1)
hold off
plot(freqArray,abs(fft_sig_0(1:2048)))
grid on
title('FFT of H1 signal 50-900Hz during non-GW block. Impulses removed')

subplot(3,1,2)
hold off
plot(freqArray,abs(fft_sig_1(1:2048)))
grid on
title('FFT of H1 signal 50-900Hz during GW block. Impulses removed')

subplot(3,1,3)
hold off
plot(freqArray,abs(fft_sig_NR(1:2048)))
grid on
title('FFT of NR template')
xlabel('Frequency in Hz')

freqArray=[10:2:2048-1];

figure(23)
subplot(3,1,1)
hold off
loglog(freqArray,abs(fft_sig_0(5+1:2048/2)))
grid on
title('FFT of H1 signal 50-900Hz during non-GW block. Impulses removed')

subplot(3,1,2)
hold off
loglog(freqArray,abs(fft_sig_1(5+1:2048/2)))
grid on
title('FFT of H1 signal 50-900Hz during GW block. Impulses removed')

subplot(3,1,3)
hold off
loglog(freqArray,abs(fft_sig_NR(5+1:2048/2)))
grid on
title('FFT of NR template')
xlabel('Frequency in Hz')


