
% Script for plotting the eye diagram of the system in which the transmitted signal is
% filtered by the RC and LCP.

clc
N  = 10^5; % number of symbols, solo se tomaran las primeras 1000000 muestras
am = 2*(rand(1,N)>0.5)-1 + j*(2*(rand(1,N)>0.5)-1); %random binary sequence BPSK
fs = 10; % sampling frequency in Hz
alpha = 0.35; %roll-off factor
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%RAISED COSINE FILTER

% defining the sinc filter
sincNum = sin(pi*[-fs:1/fs:fs]); % numerator of the sinc function
sincDen = (pi*[-fs:1/fs:fs]); % denominator of the sinc function
sincDenZero = find(abs(sincDen) < 10^-10);
sincOp = sincNum./sincDen;
sincOp(sincDenZero) = 1; % sin(pix/(pix) =1 for x =0

%Raised Cosine Filter
alpha = 0.35;
cosNum = cos(alpha*pi*[-fs:1/fs:fs]);
cosDen = (1-(2*alpha*[-fs:1/fs:fs]).^2);
cosDenZero = find(abs(cosDen)<10^-10);
cosOp = cosNum./cosDen;
cosOp(cosDenZero) = pi/4;

gt_alpha35 = sincOp.*cosOp; %RC impulse responde
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%Linear Combination Pulse
% constante*PLP(n=1) +(1-constante)RC
constante=1.7;  %valor de la constante --> referencia

Exp1_Num=((1-constante)*sin(pi*[-fs:1/fs:fs]).*cos(pi*alpha*[-fs:1/fs:fs]));
Exp1_Den=pi*[-fs:1/fs:fs].*(1-(2*alpha*[-fs:1/fs:fs]).^2);
Exp1DenZero=find(abs(Exp1_Den) < 10^-10);
Exp1_Total=Exp1_Num./Exp1_Den;
Exp1_Total(Exp1DenZero)=0.5;

Exp2_Num=constante*sin(pi*[-fs:1/fs:fs]).*sin(pi*alpha*[-fs:1/fs:fs]);
Exp2_Den=pi^2*alpha*[-fs:1/fs:fs].^2;
Exp2DenZero=find(abs(Exp2_Den) < 10^-10);
Exp2_Total=Exp2_Num./Exp2_Den;
Exp2_Total(Exp2DenZero)=0.5;

Linear_Combination_Pulse_PLP=(Exp1_Total+Exp2_Total);
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% upsampling the transmit sequence 
amUpSampled = [am;zeros(fs-1,length(am))];
amU = amUpSampled(:).';
%--------------------------------------------------------------------------
% filtered sequence using convolution fucntion
st_alpha35 = conv(amU,gt_alpha35);
Linear_Combination_Pulse_PLP=conv(amU,Linear_Combination_Pulse_PLP);
%--------------------------------------------------------------------------
%One only takes the first 1000000 samples for 
st_alpha35 = st_alpha35([1:1000000]);
Linear_Combination_Pulse_PLP=Linear_Combination_Pulse_PLP([1:1000000]);

st_alpha35_reshape = reshape(st_alpha35,fs*2,N*fs/20).';
Linear_Combination_Pulse_PLP_reshape=reshape(Linear_Combination_Pulse_PLP,fs*2,N*fs/20).';
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%PLOT of Eye Diagrams
close all
figure(1);
plot([-0.99:1/fs:0.99],real(st_alpha35_reshape).','k');   
title('Eye Diagram of the RC pulse for an excess bandwidth alpha=0.35','LineWidth',.1);
xlabel('t/T')
ylabel('Amplitude') 
axis([-1 1 -1.75 1.75])
grid on
hold on


figure(2);
plot([-0.99:1/fs:0.99],real(Linear_Combination_Pulse_PLP_reshape).','k','LineWidth',.1);   
title('Eye Diagram of the LCP pulse for an excess bandwidth alpha=0.35');
xlabel('t/T')
ylabel('Amplitude') 
axis([-1 1 -1.75 1.75])
grid on
hold on






