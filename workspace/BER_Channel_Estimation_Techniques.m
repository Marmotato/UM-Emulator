
% % % % % % % % % % % % Initialization
clc;
clear all;

snrr=[0:1:40];

nb=1593; %number of symbols symbols
N=10;
pilot=1+1i; % pilot symbol
% % % % % % % % % % % % % % % % % % % %
% % % % % % % Message symbol Generation

% % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % %mapping 16QAM
inphase = [0.9 0.9 0.9 0.9 0.3 0.3 0.3 0.3];
quadr = [-0.9 -0.3 0.3 0.9 -0.9 -0.3 0.3 0.9];
inphase = [inphase;-inphase]; inphase = inphase(:);
quadr = [quadr;quadr]; quadr = quadr(:);
const = inphase + j*quadr;

M = length(const);
b = randi([0 M-1],nb,1); %generación random de símbolos
% % % % % % % % % % % % % % % % % % % %
tx=genqammod(b,const); %general quadrature amplitude modulation
[txp]=addpilot(tx,nb,pilot,N);%adding pilot and calculating no. of symbols per frame
len=length(txp);
ct=rayleighfading(len); %genarating rayliegh fading channel
twn=txp.*ct; %Multiplying Rayleigh channel coeeficients
% % % % % % % % % % % % % % noise generation
t=awgn(twn,30,'measured','db' );%%%%% SNR
% % % % % % % % % % % % % % % % % % % % % % %
[rx rxp]=extractpilot(t,N,nb); %extracting pilot at receiver
%-----------------------------------------------------------------
%Channel Estimation CSI
csi=rxp/pilot; %calculating CSI of pilot
% % % % % % % % % % interpolation techniques
chnstinf=interpft(csi,nb);%%%%%%%% fft interpolation

t1=1:1:1770;
t2=1:10:1770;%%%%%% pilot symbols position
m=1;k=0;
for i=1:10:1770
 t3(m:m+8)=t1( i+1:i+9);%%%%%%% msg symbols position
 m=m+9;
 k=k+1;
end
chnstinf2=interp1(t2,csi,t3,'spline');%%%%%%% spline cubic interpolation
chnstinf3=interp1(t2,csi,t3,'linear');%%%%%%% linear interpolation 
chnstinf4=interp1(t2,csi,t3,'pchip');%%%%%%%% cubic interpolation technique

% % % % % % % % % % %Calculating received signal --> ECUALIZACIÓN ZERO
% FORCING
for i=1:nb
RX(i)=rx(i)/chnstinf(i); %%%%calculating rxd symbols using fft interpolation technique
RX2(i)=rx(i)/chnstinf2(i);%%% using spline cubic interpolation
RX3(i)=rx(i)/chnstinf3(i);%%%%using linear interpolation
RX4(i)=rx(i)/chnstinf4(i);%%%% using cubic interpolation
end

[ch chp]=extpltcoef(ct,N,nb); %Extracting actual pilot coefficients and channel coefficients

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% plots %%%%%%%%%%%%%%%%%%%%%%%

QAM_noise=awgn(txp,15,'measured','db' );%%%%% SNR
figure,plot(real(QAM_noise),imag(QAM_noise),'x');
title('QAM AFFECTED WITH NOISE');
xlabel('REAL(DATA)');
ylabel('IMG(DATA)');

figure,plot(real(twn),imag(twn),'r.');
title('QAM AFFECTED WITH RAYLEIGH FADING');
xlabel('REAL(DATA)');
ylabel('IMG(DATA)');

figure,plot(real(rx),imag(rx),'r.');
title('QAM PLOT WITH NOISE AND RAYLEIGH FADING');
xlabel('REAL(DATA)');
ylabel('IMG(DATA)');

figure,
plot(real(RX),imag(RX),'r.');
hold on;
plot(real(tx),imag(tx),'o');
grid on;
title('QAM PLOT');
xlabel('REAL(DATA)');
ylabel('IMG(DATA)');
legend('PLOT AT RX con Ecualización','PLOT AT TX' );

%--------------------------------------------------------------------------------

%Función de transferencia del canal  y estimación del canal
cdb=10*log(abs(ch));
figure,plot(0:1:nb-1,cdb);
title('Typical Rayleigh fading channel');
xlabel('t in samples');
ylabel('Power in dB');

figure,plot(0:1:nb-1,(abs(ch)),'r'); %Plotting actual value
hold on;
plot(0:1:nb-1,(abs(chnstinf)),'b'); %Plotting estimated value of fft
legend ('Actual value','Estimated value fft');
hold off;
title('ACTUAL AND ESTIMATED CHANNEL COEFFICIENTS with FFT');
xlabel('Time in samples');
ylabel('Magnitude of coefficients');


figure,plot(0:1:nb-1,(abs(ch)),'r'); %Plotting actual value
hold on;
plot(0:1:nb-1,(abs(chnstinf2)),'b'); %Plotting estimated value by spline
legend ('Actual value','Estimated value spline');
hold off;
title('ACTUAL AND ESTIMATED CHANNEL COEFFICIENTS with Spline');
xlabel('Time in samples');
ylabel('Magnitude of coefficients');

figure,plot(0:1:nb-1,(abs(ch)),'r'); %Plotting actual value
hold on;
plot(0:1:nb-1,(abs(chnstinf3)),'b'); %Plotting estimated value by linear
legend ('Actual value','Estimated value linear');
hold off;
title('ACTUAL AND ESTIMATED CHANNEL COEFFICIENTS with Linear');
xlabel('Time in samples');
ylabel('Magnitude of coefficients');

figure,plot(0:1:nb-1,(abs(ch)),'r'); %Plotting actual value
hold on;
plot(0:1:nb-1,(abs(chnstinf4)),'b'); %Plotting estimated value by cubic
legend ('Actual value','Estimated value fft');
hold off;
title('ACTUAL AND ESTIMATED CHANNEL COEFFICIENTS with Cubic');
xlabel('Time in samples');
ylabel('Magnitude of coefficients');
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%% for BER curves%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% generating curves for snr ranging from 0dB to 40dB
for l=1:length(snrr)
tt=awgn(twn,snrr(l),'measured','db' );
% % % % % % % % % % % % % % % % % % % % % % %
[rxt rxpt]=extractpilot(tt,N,nb); %extracting pilot at receiver
csit=rxpt/pilot; %calculating CSI of pilot
% % % % % % % % % % interpolation techniques

chnstinft=interpft(csit,nb);
t1t=1:1:1770;
t2t=1:10:1770;
mt=1;kt=0;
for i=1:10:1770
 t3t(mt:mt+8)=t1t( i+1:i+9);
 mt=mt+9;
 kt=kt+1;
end
chnstinf2t=interp1(t2t,csit,t3t,'spline');
chnstinf3t=interp1(t2t,csit,t3t,'linear');
chnstinf4t=interp1(t2t,csit,t3t,'pchip');
% % % % % % % % % % %Calculating received signal ECUALIZACIÓN ZF
for i=1:nb
RXt(i)=rxt(i)/chnstinft(i);
RX2t(i)=rxt(i)/chnstinf2t(i);
RX3t(i)=rxt(i)/chnstinf3t(i);
RX4t(i)=rxt(i)/chnstinf4t(i);
end

%Demodulación de la señal
rt=genqamdemod(RXt,const);
r2t=genqamdemod(RX2t,const);
r3t=genqamdemod(RX3t,const);
r4t=genqamdemod(RX4t,const);
[no_of_error1(l),rate1(l)]=biterr(b.',rt) ; % error rate calculation for fft
[no_of_error2(l),rate2(l)]=biterr(b.',r2t) ; % error rate calculation for spline
[no_of_error3(l),rate3(l)]=biterr(b.',r3t) ; % error rate calculation for linear
[no_of_error4(l),rate4(l)]=biterr(b.',r4t) ; % error rate calculation for cubic
end
% % % % % % % % % % % % % BER plot % % % % % % % % % % %
figure,semilogy(snrr,rate1,'b-',snrr,rate2,'r-',snrr,rate3,'k-',snrr,rate4,'g-');
legend('fft','cubic spline','linear','cubic');
title('BER curves for different interpolation techniques');
xlabel('SNR in dB');
ylabel('BER');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%used functions
%function [txp]=addpilot(tx,nb,pilot,N) %Number of symbols per frame (Interpolation Factor)
%m=1;k=0;
%for i=1:N-1:nb
%txp(m)=pilot;
%txp(m+1:m+N-1)=tx(i:i+N-2);
%m=m+N;
%k=k+1;
%end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Remove pilot coefficients -- We actually do not use this function.
%function [ch chp]=extpltcoef(ct,N,nb)
%m=1;
%ch=[];
%for i=1:N:length(ct)
%chp(m)=ct(i); %extracting only pilot symbols' coeeficient
%ch=[ch ct(i+1:i+N-1)];
%m=m+1;
%end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Extraxting and removing pilot from transmitted signal
%function [rx rxp]=extractpilot(t,N,nb)
%m=1;
%rx=[];
%for i=1:N:length(t)
%rxp(m)=t(i); %extracting only pilot symbols
%rx=[rx t(i+1:i+N-1)]; %extracting data symbols
%m=m+1;
%end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Rayleigh fading
%function [c]=rayleighfading(m)
%clc;
%N=50; % Number of reflections
%fmax=100; %Max doppler shift
%A=1; %amplitude
%f=10000; %sampling frequency
%t=0:1/f:((m/10000)-(1/f)); %sampling time
%ct=zeros(1,m);
%ph=2*pi* rand(1,32);
%theta=2*pi*rand(1,32);
%fd=fmax*cos(theta); %doppler shift
%for n=1:m
%for i=1:32
%ct(n)=ct(n)+(A*exp(j*(2*pi*fd(i)*t(n)+ph(i))));
%end
%end
%c=ct/sqrt(N); %channel coefficient
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%