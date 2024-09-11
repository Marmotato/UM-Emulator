function ber = get_ber(modulacion, escenario, metodo, distancia, snr, show_plot)
% parametros
n_ref = escenario(2);
vel = escenario(3);
freq = escenario(4);
rayleigh = escenario(5);
N=distancia; %intervalo pilotos
num_bits = 100000; % n° bits a considerar


% generación de modulación
if modulacion == "16QAM"
    mod_order = 4;
    nb = num_bits/mod_order; % n° de simbolos
    if metodo == "fft" || metodo == "cubic"
        nb = floor(nb/(N-1))*(N-1);
    end
    A1 = sqrt(1/10);
    A2 = 3*sqrt(1/10);
    pilot=A1+1j*A1; % pilot symbol
    inphase = [A2 A2 A2 A2 A1 A1 A1 A1];
    quadr = [-A2 -A1 A1 A2 -A2 -A1 A1 A2];
    inphase = [inphase;-inphase]; inphase = inphase(:);
    quadr = [quadr;quadr]; quadr = quadr(:);
    const = inphase + j*quadr;
    M = length(const);
    b = randi([0 M-1],nb,1); %generación random de símbolos
    tx=genqammod(b,const); %general quadrature amplitude modulation
elseif modulacion == "QPSK"
    mod_order = 2;
    nb = num_bits/mod_order;
    if metodo == "fft" || metodo == "cubic"
        nb = floor(nb/(N-1))*(N-1);
    end
    A = sqrt(1/2);
    pilot = A+1j*A;
    inphase = [A A -A -A];
    quadr = [A -A A -A];
    const = inphase + j*quadr;
    M = length(const);
    b = randi([0 M-1],nb,1); %generación random de símbolos
    tx=genqammod(b,const); %general quadrature amplitude modulation
else
    disp("Error en argumento modulacion")
end

% adicion de pilotos
if metodo == "fft" || metodo == "cubic"
    [txp]=addpilot(tx,nb,pilot,N);
else
    for i=1:1:nb
        txp(i) = tx(i);
    end
end
len=length(txp);

% generación del canal
if rayleigh == 1
    ct = ( randn( 1, len ) + j*randn( 1, len ) )*sqrt(1/2);
elseif rayleigh == 0
    ct = doppler_fading(len, n_ref, vel, freq);
else
    disp("Error en argumento rayleigh")
end

% simulación de transmición
twn=txp.*ct; %Multiplying Rayleigh channel coeficients
Tx=awgn(twn,snr,'measured','db' ); % noise generation

% extracción de pilotos
if metodo == "fft" || metodo == "cubic"
    [rx rxp]=extractpilot(Tx,N,nb); %extracting pilot at receiver
    csi=rxp/pilot; %calculating CSI of pilot
elseif metodo == "perfect"
    rx = Tx;
else
    disp("Error en argumento metodo")
end

% interpolation para estimar canal
if metodo == "fft"
    channel=interpft(csi,nb);
elseif metodo == "cubic"
    aux = floor(nb/(N-1))+nb;
    t1=1:1:aux;
    t2=1:N:aux; %pilot symbols position
    m=1;k=0;
    for i=1:N:aux
        t3(m:m+N-2)=t1( i+1:i+N-1); %msg symbols position
        m=m+N-1;
        k=k+1;
    end
    channel=interp1(t2,csi,t3,'spline');
else
    channel = ct;
end

% ecualizador
for i=1:nb
RX(i)=rx(i)/channel(i);
end

% demodulación
rt=genqamdemod(RX,const);
[no_of_error1,rate1]=biterr(b.',rt) ; % error rate calculation for fft
ber = rate1;

% graficos
if show_plot
    if metodo ~= "perfect"
        [ch chp]=extpltcoef(ct,N,nb); %Extracting actual pilot coefficients and channel coefficients
    else
        ch = ct;
    end

    % SNR = -5dB
    snr_value = -5;
    % solo ruido
    txp_noise=awgn(txp,snr_value,'measured','db' );%%%%% SNR
    figure,plot(real(txp_noise),imag(txp_noise),'x', Color=[0.6 0 0.4]);
    title("Constelación recibida afectada con AWGN, SNR = "+string(snr_value)+"dB");
    xlabel('REAL(DATA)');
    ylabel('IMG(DATA)');
    % ruido mas canal antes
    txp_ch = txp.*ct;
    txp_noise=awgn(txp_ch,snr_value,'measured','db' );
    figure,plot(real(txp_noise),imag(txp_noise),'x',  Color = [0 0.7 0.3]);
    title("Constelación recibida afectada por el canal y AWGN, SNR = "+string(snr_value)+"dB");
    xlabel('REAL(DATA)');
    ylabel('IMG(DATA)');
    % ruido mas canal despues
    [rx rxp]=extractpilot(txp_noise,N,nb); %extracting pilot at receiver
    csi=rxp/pilot; %calculating CSI of pilot
    channel=interpft(csi,nb);
    for i=1:nb
    RX(i)=rx(i)/channel(i);
    end
    figure,plot(real(RX),imag(RX),'x',  Color = [0 0.3 0.9]);
    title("Constelación recibida afectada por el canal y AWGN despues de ecualizar, SNR = "+string(snr_value)+"dB");
    xlabel('REAL(DATA)');
    ylabel('IMG(DATA)');

    % SNR = 0dB
    snr_value = 0;
    % solo ruido
    txp_noise=awgn(txp,snr_value,'measured','db' );%%%%% SNR
    figure,plot(real(txp_noise),imag(txp_noise),'x',  Color=[0.6 0 0.4]);
    title("Constelación recibida afectada con AWGN, SNR = "+string(snr_value)+"dB");
    xlabel('REAL(DATA)');
    ylabel('IMG(DATA)');
    % ruido mas canal antes
    txp_ch = txp.*ct;
    txp_noise=awgn(txp_ch,snr_value,'measured','db' );
    figure,plot(real(txp_noise),imag(txp_noise),'x', Color = [0 0.7 0.3]);
    title("Constelación recibida afectada por el canal y AWGN, SNR = "+string(snr_value)+"dB");
    xlabel('REAL(DATA)');
    ylabel('IMG(DATA)');
    % ruido mas canal despues
    [rx rxp]=extractpilot(txp_noise,N,nb); %extracting pilot at receiver
    csi=rxp/pilot; %calculating CSI of pilot
    channel=interpft(csi,nb);
    for i=1:nb
    RX(i)=rx(i)/channel(i);
    end
    figure,plot(real(RX),imag(RX),'x', Color = [0 0.3 0.9]);
    title("Constelación recibida afectada por el canal y AWGN despues de ecualizar, SNR = "+string(snr_value)+"dB");
    xlabel('REAL(DATA)');
    ylabel('IMG(DATA)');

    % SNR = 15dB
    snr_value = 15;
    % solo ruido
    txp_noise=awgn(txp,snr_value,'measured','db' );%%%%% SNR
    figure,plot(real(txp_noise),imag(txp_noise),'x', Color=[0.6 0 0.4]);
    title("Constelación recibida afectada con AWGN, SNR = "+string(snr_value)+"dB");
    xlabel('REAL(DATA)');
    ylabel('IMG(DATA)');
    % ruido mas canal antes
    txp_ch = txp.*ct;
    txp_noise=awgn(txp_ch,snr_value,'measured','db' );
    figure,plot(real(txp_noise),imag(txp_noise),'x', Color = [0 0.7 0.3]);
    title("Constelación recibida afectada por el canal y AWGN, SNR = "+string(snr_value)+"dB");
    xlabel('REAL(DATA)');
    ylabel('IMG(DATA)');
    % ruido mas canal despues
    [rx rxp]=extractpilot(txp_noise,N,nb); %extracting pilot at receiver
    csi=rxp/pilot; %calculating CSI of pilot
    channel=interpft(csi,nb);
    for i=1:nb
    RX(i)=rx(i)/channel(i);
    end
    figure,plot(real(RX),imag(RX),'x', Color = [0 0.3 0.9]);
    title("Constelación recibida afectada por el canal y AWGN despues de ecualizar, SNR = "+string(snr_value)+"dB");
    xlabel('REAL(DATA)');
    ylabel('IMG(DATA)');

    % SNR = 30dB
    snr_value = 30;
    % solo ruido
    txp_noise=awgn(txp,snr_value,'measured','db' );%%%%% SNR
    figure,plot(real(txp_noise),imag(txp_noise),'x',  Color=[0.6 0 0.4]);
    title("Constelación recibida afectada con AWGN, SNR = "+string(snr_value)+"dB");
    xlabel('REAL(DATA)');
    ylabel('IMG(DATA)');
    % ruido mas canal antes
    txp_ch = txp.*ct;
    txp_noise=awgn(txp_ch,snr_value,'measured','db' );
    figure,plot(real(txp_noise),imag(txp_noise),'x',  Color = [0 0.7 0.3]);
    title("Constelación recibida afectada por el canal y AWGN, SNR = "+string(snr_value)+"dB");
    xlabel('REAL(DATA)');
    ylabel('IMG(DATA)');
    % ruido mas canal despues
    [rx rxp]=extractpilot(txp_noise,N,nb); %extracting pilot at receiver
    csi=rxp/pilot; %calculating CSI of pilot
    channel=interpft(csi,nb);
    for i=1:nb
    RX(i)=rx(i)/channel(i);
    end
    figure,plot(real(RX),imag(RX),'x', Color = [0 0.3 0.9]);
    title("Constelación recibida afectada por el canal y AWGN despues de ecualizar, SNR = "+string(snr_value)+"dB");
    xlabel('REAL(DATA)');
    ylabel('IMG(DATA)');
    
%     figure,plot(real(twn),imag(twn),'r.');
%     title('Símbolos afectados por el canal');
%     xlabel('REAL(DATA)');
%     ylabel('IMG(DATA)');
%     
%     figure,plot(real(rx),imag(rx),'r.');
%     title('Símbolos afectados por el canal y ruido AWGN');
%     xlabel('REAL(DATA)');
%     ylabel('IMG(DATA)');
%     
%     figure,
%     plot(real(RX),imag(RX),'r');
%     hold on;
%     plot(real(tx),imag(tx),'o');
%     grid on;
%     title('QAM PLOT');
%     xlabel('REAL(DATA)');
%     ylabel('IMG(DATA)');
%     legend('PLOT AT RX con Ecualización','PLOT AT TX' );
% 
    %Función de transferencia del canal  y estimación del canal
    cdb=10*log(abs(ch));
    figure,plot(0:1:nb-1,cdb, LineWidth=1.3, Color = [0 0.1 0.6]);
    title('Canal en función del tiempo');
    subtitle("n_{ref}="+string(n_ref)+", "+"V="+string(vel)+" km/h"+", "+"f_c="+string(freq/1e6)+" MHz");
    %subtitle("Infinitos rebotes Rayleigh multiplicativo")
    xlabel('t in samples');
    ylabel('Power in dB');
%     
%     figure,plot(0:1:nb-1,(abs(ch)),color ='r'); %Plotting actual value
%     hold on;
%     plot(0:1:nb-1,(abs(channel)),'b'); %Plotting estimated value of fft
%     legend ('Actual value','Estimated value fft');
%     hold off;
%     title('ACTUAL AND ESTIMATED CHANNEL COEFFICIENTS with FFT');
%     xlabel('Time in samples');
%     ylabel('Magnitude of coefficients');
end

end