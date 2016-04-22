%Compute Power in 50-300Hz region

%inFile=fopen('H-H1_LOSC_4_V1-1126257414-4096.txt','r');
inFile=fopen('L-L1_LOSC_4_V1-1126257414-4096.txt','r');

str_1= fgets(inFile);
str_1= fgets(inFile);
str_1= fgets(inFile);

N1=65536*2*128; 
for count=1:N1,

  dataArray_H(count)=fscanf(inFile,'%f',1);
end

fclose(inFile);


Fs=4096;
timeArray=[0:1:N1-1]/Fs;


%30-900Hz filter   
filterOrder=4; 
Wn=[30 900]*2/Fs;
[B,A] = cheby2(filterOrder,50,Wn,'bandpass');



%dataArray_H_filtered = dataArray_H; %
dataArray_H_filtered = filtfilt(B,A,dataArray_H);

sig=dataArray_H_filtered;





%Compute power, starting from center


N0=66680+1024*4*2032-1024*4*1000; 
N1=1024*4; %1 Hz resolution

%60*n Hz Tones;
toneFreqArray=[60 120 180 ]; 
toneIndex= toneFreqArray/1 + 1;
toneIndexOffset=1; 


for count0=1:1000*2,
    
 %FFT
 sig=dataArray_H_filtered( N0+(count0-1)*N1+(1:N1) );
 f_step=Fs/length(sig); %f_step=Fs/(65536*2);
 fft_sig=fft(sig);



 pwrArray_0(count0)=sqrt( sum( abs(sig).^2 )/length(sig) );

toneRegionArray=zeros(1,length(fft_sig)/2);
for count=1:length(toneIndex),
toneRegionArray(toneIndex(count)+(-toneIndexOffset:toneIndexOffset))= ones(1,toneIndexOffset*2+1)*max(abs(fft_sig));
end

%Low band

N3=round(300/f_step);
toneRegionLowArray=zeros(1,N3);
for count=1:3,
toneRegionLowArray(toneIndex(count)+(-toneIndexOffset:toneIndexOffset))= ones(1,toneIndexOffset*2+1);
end


toneRegionLowArrayIndex=find( toneRegionLowArray == 1);

toneRegionLowArray(1:50)= -1;

signalRegionLowArrayIndex = find( toneRegionLowArray == 0);

%60,120,180Haz tones
pwrArray_1(count0) =  sum( abs(fft_sig(toneRegionLowArrayIndex)).^2 )/length(sig) ;

%50-300Hz excluding 60,120,180Hz
pwrArray_3(count0) =  sum( abs(fft_sig(signalRegionLowArrayIndex)).^2 )/length(sig) ;

%50-300Hz including 60,120,180Hz
pwrArray_5(count0) = sum( abs(fft_sig([ toneRegionLowArrayIndex signalRegionLowArrayIndex] )).^2 )/length(sig) ;



end %end for


timeArray_1=[-1000:1:1000-1];

figure(1)
hold off
plot(timeArray_1,pwrArray_1);
%title('toneRegionLowArray')
title('Power in 60Hz,120Hz and 180Hz. Averaged over 1 second.')
xlabel('Time in seconds. Time=0 corresponds to GW150914 region.')
grid on

figure(3)
hold off
plot(timeArray_1,pwrArray_3);
%title('signalRegionLowArray')
title('Power in 50Hz-300Hz range, excluding 60*n Hz tones.')
xlabel('Time in seconds. Time=0 corresponds to GW150914 region.')
text(100,3.5e-41,'Averaged over 1 second.')
%gtext('Averaged over 1 second.')
grid on

figure(5)
hold off
plot(timeArray_1,pwrArray_5);
%title('tone + signalRegionLowArray')
title('Power in 50Hz-300Hz range, including 60*n Hz tones.')
xlabel('Time in seconds. Time=0 corresponds to GW150914 region.')
text(100,5.8e-41,'Averaged over 1 second.')
%gtext('Averaged over 1 second.')
grid on

disp('Maximum Power in 60*n Hz tones...')
max(pwrArray_1)

sig_1=[pwrArray_3(1:1000) pwrArray_3(1001+(1:1000-1))]; %except mid block containing GW event
sig_2= pwrArray_5(1001) - mean(sig_1);
disp('EM component power in 50-300Hz region....')
sig_2

disp('Maximum amplitude in 60*n Hz Tones...')
sqrt(max(pwrArray_1))

disp('EM component amplitude  in 50-300Hz region....')
sqrt(sig_2)

disp('Percentage of Excess amplitude in 50-300Hz region....')
abs(sqrt(sig_2)-sqrt(max(pwrArray_1)))*100/sqrt(max(pwrArray_1))



