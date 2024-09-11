% escenarios a evaluar
escenario1 = [1 5 30 700*1e6 0];
escenario2 = [2 5 30 3.5*1e9 0];
escenario3 = [3 5 120 700*1e6 0];
escenario4 = [4 5 120 3.5*1e9 0];
escenario5 = [5 40 30 700*1e6 0];
escenario6 = [6 40 30 3.5*1e9 0];
escenario7 = [7 40 120 700*1e6 0];
escenario8 = [8 40 120 3.5*1e9 0];
escenario9 = [9 0 0 0 1];

% extraccion de data
qam_1 = load("results/qam_1.mat").qam_1;
qam_2 = load("results/qam_2.mat").qam_2;
qam_3 = load("results/qam_3.mat").qam_3;
qam_4 = load("results/qam_4.mat").qam_4;
qam_5 = load("results/qam_5.mat").qam_5;
qam_6 = load("results/qam_6.mat").qam_6;
qam_7 = load("results/qam_7.mat").qam_7;
qam_8 = load("results/qam_8.mat").qam_8;
qam_9 = load("results/qam_9.mat").qam_9;

qpsk_1 = load("results/qpsk_1.mat").qpsk_1;
qpsk_2 = load("results/qpsk_2.mat").qpsk_2;
qpsk_3 = load("results/qpsk_3.mat").qpsk_3;  
qpsk_4 = load("results/qpsk_4.mat").qpsk_4;
qpsk_5 = load("results/qpsk_5.mat").qpsk_5;
qpsk_6 = load("results/qpsk_6.mat").qpsk_6;
qpsk_7 = load("results/qpsk_7.mat").qpsk_7;
qpsk_8 = load("results/qpsk_8.mat").qpsk_8;
qpsk_9 = load("results/qpsk_9.mat").qpsk_9;

%graficos
plot_ber_by_modulation(qam_1, "16QAM", qpsk_1, "QPSK", escenario1);
plot_ber_by_modulation(qam_2, "16QAM", qpsk_2, "QPSK", escenario2);
plot_ber_by_modulation(qam_3, "16QAM", qpsk_3, "QPSK", escenario3);
plot_ber_by_modulation(qam_4, "16QAM", qpsk_4, "QPSK", escenario4);
plot_ber_by_modulation(qam_5, "16QAM", qpsk_5, "QPSK", escenario5);
plot_ber_by_modulation(qam_6, "16QAM", qpsk_6, "QPSK", escenario6);
plot_ber_by_modulation(qam_7, "16QAM", qpsk_7, "QPSK", escenario7);
plot_ber_by_modulation(qam_8, "16QAM", qpsk_8, "QPSK", escenario8);
plot_ber_by_modulation(qam_9, "16QAM", qpsk_9, "QPSK", escenario9);

function plot_ber_curves(ber_list, name)
% plot_ber_curves(qam_1, "qam_1")
    snrs = -2:1:30;
    c1 = ber_list{1};
    c2 = ber_list{2};
    c3 = ber_list{3};
    c4 = ber_list{4};
    c5 = ber_list{5};
    c6 = ber_list{6};
    c7 = ber_list{7};
    figure,semilogy(snrs,c1, snrs,c2,snrs,c3,snrs,c4,snrs,c5,snrs,c6,snrs,c7);
    legend('fft-5','fft-10','fft-20','cubic-5','cubic-10','cubic-20', 'perfect');
    title('BER for '+name);
    xlabel('SNR in dB');
    ylabel('BER');
end

function plot_ber_by_modulation(mod1, name1, mod2, name2, escenario)
n_ref = escenario(2);
vel = escenario(3);
freq = escenario(4);
rayleigh = escenario(5);
snrs = -2:1:30;

m1c1 = mod1{1};
m1c2 = mod1{2};
m1c3 = mod1{3};
m1c4 = mod1{4};
m1c5 = mod1{5};
m1c6 = mod1{6};
m1c7 = mod1{7};
m2c1 = mod2{1};
m2c2 = mod2{2};
m2c3 = mod2{3};
m2c4 = mod2{4};
m2c5 = mod2{5};
m2c6 = mod2{6};
m2c7 = mod2{7};

for i=2:1:length(snrs)
    snr = 10^(snrs(i)/10);
    ber_teo(i) = 0.5*(1-sqrt(abs(snr/(1+snr))));
end

figure;
semilogy(snrs,m1c1,snrs,m1c2,snrs,m1c3,snrs,m1c4,snrs,m1c5,snrs,m1c6,snrs,m1c7, LineWidth=1.5)
hold on
semilogy(snrs,m2c1,'--',snrs,m2c2,'--',snrs,m2c3,'--',snrs,m2c4,'--',snrs,m2c5,'--',snrs,m2c6,'--',snrs,m2c7,'--', LineWidth=1.5)
hold on
semilogy(snrs, ber_teo, 'k', LineWidth=1.5);
hold on
legend(name1+'-fft-5',name1+'-fft-10',name1+'-fft-20',name1+'-cubic-5',name1+'-cubic-10',name1+'-cubic-20', name1+'-perfect', name2+'-fft-5',name2+'-fft-10',name2+'-fft-20',name2+'-cubic-5',name2+'-cubic-10',name2+'-cubic-20', name2+'-perfect', 'BER Teórico');
title("BER para escenario N° "+escenario(1));
if rayleigh == 0
    subtitle("n_{ref}="+string(n_ref)+", "+"V="+string(vel)+" km/h"+", "+"f_c="+string(freq/1e6)+" MHz");
else
    subtitle("Caso Rayleigh Fading con infinitos rebotes");
end
xlabel('SNR in dB');
ylabel('BER');

end