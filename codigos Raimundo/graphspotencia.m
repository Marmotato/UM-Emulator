q=1.6e-19; % carga del electron 
R=0.72; % responsividad del fotodiodo
B=100e6;%Ancho de banda del filtro eléctrico que sigue el fotodetector [Hz] ¿¿¿¿¿¿
Ib=10e-6;%Corriente de fondo debida a la luz ambiente
I2=0.562;%Factor de ancho de banda de transductancia


k=1.38e-23;%Boltzman constant [m^2 kg s^-2 K^-1]
tk=295;%Absolute temperature
I3=0.0868;%Factor de ancho de banda
gol=1;%open loop voltage gain  
cpd=112e-8; %Capacitancia del fotodetector por unidad de área
Ap=0.0001;
FET=1.5;%FET channel noise 
gm=0.03;% FET transconductancia
%Señal eléctrica recibida

P_los=m_HLoS1;
P_nlos=ret+ret1;
P_sca=Hsca;
P_recibida=P_los+P_nlos+P_sca;

for ii=1:length(P_recibida)
for jj=1:length(P_recibida)
rdsht(ii,jj)=(2*q*R*P_recibida(ii,jj)*B)+(2*q*Ib*I2*B);%Ruido shot 
rdthrml=((8*pi*k*tk*cpd*Ap*I2*B^2)/gol)+(((16*pi^2*k*tk*FET)/gm)*cpd^2*Ap^2*I3*B^3);%Ruido termal
ruido(ii,jj)=rdsht(ii,jj)+rdthrml;
end
end

P_total=P_recibida+ruido;
P_totalv=reshape(10.*P_total,1,[]);
[f,tim] = ecdf(P_totalv);
figure
ecdf(P_totalv)

figure
surf(x,y, 10.*P_total);
axis([0 lx 0 ly min(min(10.*P_total)) max(max(10.*P_total))]);

P_rec_dBm=(10*log10(10.*P_total));
P_dbmv=reshape(P_rec_dBm,1,[]);
[f,tim] = ecdf(P_dbmv);
figure
ecdf(P_dbmv)

figure
surf(x,y, P_rec_dBm);
axis([0 lx 0 ly min(min(P_rec_dBm)) max(max(P_rec_dBm))]);







