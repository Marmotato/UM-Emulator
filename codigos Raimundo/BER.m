clear all;
clc;
close all;
% SNR Performance For VLC SYSTEM %
%
%
% BASIC PARAMETER REQUIRED %
Incidence = 70*pi/180;
TX_FOV = 70; % Transmitter Field Of View
RX_FOV = 90; % Receivers Field Of View
Tx = [2,2,2]; % Transmitter Location
%Rxp = [2,2]; % Receiver Location
W_Room = 4; % Width of Room
L_Room = 4; % Length of Room
H_Room = 2; % Height Between Transmitter and Receiver
R = 1; % Responsivity of Photodiode
Apd = 1e-4; % Area of PhotoDetector
Rb = 1e6; % Data rate of system
Iamp = 5e-12; % Amplifier Current
q = 1.6e-19; % Electron Charge
Bn = 50e6; % Noise Bandwidth
I2 = 0.562; % Noise Bandwidth Factor
PLED = 1; % Power Emitted by LED
index =1;
HLED = 1;
[W,L] = meshgrid(-(W_Room/2) : 0.50 : (W_Room/2)); % Consideer Length of BLock for Room
xydist = sqrt((W).^2 + (L).^2);
hdist = sqrt(xydist.^2 + HLED.^2);
%D = Tx – Rx;
%d = norm(D);
%Incidence = acos()
A_Irradiance = ((Tx(3)-HLED)./hdist);
%I(index) = Irradiance*180/pi;
%if abs(Incidence <= RX_FOV)
p = TX_FOV ;
Tx_FOV = (TX_FOV*pi)/180;
% BASIC CALCULATION IN VLC SYSTEM %
% Lambertian Pattern
m = real(-log(2)/log(cos(Tx_FOV)));
% Radiation Intensity at particular point
Ro = real(((m+1)/(2*pi)).*A_Irradiance^m);
% Transmitted power By LED
Ptx = PLED .* Ro;
% Channel Gain ( Channel Coefficient Of LOS Channel )
%Theta=atand(sqrt(sum((Tx-Rx).^2))/H_Room);
HLOS = (Apd./hdist.^2).*cos(Incidence).*Ro;
% Received Power By PhotoDetector
Prx = HLOS.*Ptx;
% Calculate Noise in System
Bs = Rb*I2;
Pn = Iamp/Rb;
Ptotal = Prx+Pn;
new_shot = 2*q*Ptotal*Bs;
new_amp = Iamp^2*Bn;
% Calculate SNR
new_total = new_shot + new_amp;
SNRl = (R.*Prx).^2./ new_total;
SNRdb = 10*log10(SNRl);
% else
% SNRl = 0;
% SNRdb = 0;
% end
index = index + 1;
% Plot Graph %
ber=qfunc(sqrt(SNRl));
figure;
mesh(W,L,SNRdb);
%mesh(SNRdb);
%ylim([0 30]);
title('SNR Distribution in Room');
xlabel('Length of Room');
ylabel('Width of Room');
zlabel('SNR in dB');

figure;
mesh(W,L,ber);
%mesh(SNRdb);
%ylim([0 30]);
title('BER Distribution in Room');
xlabel('Length of Room');
ylabel('Width of Room');
zlabel('BER');