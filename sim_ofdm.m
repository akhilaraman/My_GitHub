%Nov 27, 2015
%Simulate OFDM with 4-ary signalling on each FFT sample, add Additive White Gaussian Nose(AWGN)
%Compute Bit Error Rate(BER) and compare with ideal BER.

%OFDM with M-ary PAM;  M=4;

N_fft=512;

numLoop=10*200; %Increase this number for higher BER accuracy
N=N_fft*numLoop;

count3=1;

for SNR=16:-1:4, %SNR=Signal to Noise Ratio

    disp('processing SNR in dB...')
    SNR
    %SNR=10*log10( sigPwr/noisePwr)  % average noise power=1.
    
    noisePwr=1;
    sigPwr= 10.^(SNR/10)*noisePwr;

berCount=0;
for count1=1:numLoop,
    
r0=rand(1,N_fft);
r1= r0 > .5;  %+1 or -1


r0=rand(1,N_fft);
r2= r0 > .5;  %+1 or -1

%One bit from r1 stream and one bit from r2 stream are taken, combined
%together and coded into one of 4 symbols. M-ary PAM, M=4

%Real part of FFT
for count=1:N_fft/2,
  if r1(count)==0 && r2(count)==0,
       fft_sig_real(count)=-3;
  end
  if r1(count)==0 && r2(count)==1,
       fft_sig_real(count)=-1;
  end
    if r1(count)==1 && r2(count)==1,
       fft_sig_real(count)=1;
    end
    if r1(count)==1 && r2(count)==0,
       fft_sig_real(count)=3;
    end
end

%Imaginary part of FFT
for count=N_fft/2+(1:N_fft/2),
  if r1(count)==0 && r2(count)==0,
       fft_sig_imag(count-N_fft/2)=-3;
  end
  if r1(count)==0 && r2(count)==1,
       fft_sig_imag(count-N_fft/2)=-1;
  end
    if r1(count)==1 && r2(count)==1,
       fft_sig_imag(count-N_fft/2)=1;
    end
    if r1(count)==1 && r2(count)==0,
       fft_sig_imag(count-N_fft/2)=3;
    end
end


fft_sig1=(fft_sig_real + i.*fft_sig_imag) * sqrt(N_fft/2) * sqrt(sigPwr);

fft_sig_tx =[ fft_sig1 conj(fft_sig1(N_fft/2)) conj(fft_sig1(N_fft/2:-1:2)) ];

%IFFT
sig_1 = ifft(fft_sig_tx);  %Signal in time domain is real.

noise_1= randn(1,N_fft)*sqrt(noisePwr);  %Zero mean Unit standard deviation AWGN

rx_1= sig_1 + noise_1;  %add real noise AWGN.

%FFT
fft_sig_rx = fft(rx_1);

fft_sym = [ real(fft_sig_rx(1:N_fft/2))  imag(fft_sig_rx(1:N_fft/2)) ]/( sqrt(N_fft/2)*sqrt(sigPwr));



fft_sig_tx_2 = [fft_sig_real fft_sig_imag];

%Decode FFT symbols and generate bits and compute BER
sym_1=-3; sym_2= -1; sym_3=1; sym_4=3;
symArray=[sym_1 sym_2 sym_3 sym_4];
for count=1:N_fft,
    dist_1=abs(fft_sym(count)-sym_1);
    dist_2=abs(fft_sym(count)-sym_2);
    dist_3=abs(fft_sym(count)-sym_3);
    dist_4=abs(fft_sym(count)-sym_4);
    [minVal minIndex]=min([dist_1 dist_2 dist_3 dist_4]);
    decoded_1(count)=symArray(minIndex);
    
    if decoded_1(count)~=fft_sig_tx_2(count),
        if abs(decoded_1(count)- fft_sig_tx_2(count)) ==4,
           berCount=berCount+2;
        else
            berCount=berCount+1;
        end
    end
end


end

%simulated BER
berCount_sim(count3)=berCount/(N_fft*numLoop*2);



%Ideal BER
A=sqrt(sigPwr); sigma=1;
p0=0.5* erfc( (A*1/sigma) * sqrt(1/2) );
p1=0.5* erfc( ( A*(1+2)/sigma) * sqrt(1/2) );
p2=0.5* erfc( ( A*(1+4)/sigma) * sqrt(1/2) );
p3=0.5* erfc( ( A*(1+6)/sigma) * sqrt(1/2) );
berCount_ideal(count3)=( 0.5*( p0+ 2*p1 + (p0-p1) ) + 0.5* ( (p0-p1)+(p1-p2)*2+ p2 ) )/2; %exact bit error



berCount_ideal(count3)
berCount_sim(count3)

SNRarray(count3)=SNR;
BitRatePerBandwidth=4;
EbN0array(count3)=SNR-10*log10(BitRatePerBandwidth);

count3=count3+1;

end %end loop

figure(1)
hold off
plot(SNRarray,berCount_ideal,'r-')
hold on
plot(SNRarray,berCount_sim,'b-')
grid on
title('Plot of Signal to Noise ratio Vs Bit Error Probability. Red: Ideal. Blue: simulated');
xlabel('SNR in dB');
ylabel('BER');

figure(2)
hold off
plot(EbN0array,berCount_ideal,'r-')
hold on
plot(EbN0array,berCount_sim,'b-')
grid on
title('Plot of Eb/N0 ratio Vs Bit Error Probability. Red: Ideal. Blue: simulated');
xlabel('Eb/N0 in dB');
ylabel('BER');

