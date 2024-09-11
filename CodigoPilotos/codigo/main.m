clear; clc;
% parametros generales
nro_ejecuciones = 21;
snrs = -2:1:30;
modulaciones = ["16QAM", "QPSK"];

% escenarios: [n n_ref, vel, carrier Rayleigh(bool)]
escenario1 = [1 5 30 700*1e6 0];
escenario2 = [2 5 30 3.5*1e9 0];
escenario3 = [3 5 120 700*1e6 0];
escenario4 = [4 5 120 3.5*1e9 0];
escenario5 = [5 40 30 700*1e6 0];
escenario6 = [6 40 30 3.5*1e9 0];
escenario7 = [7 40 120 700*1e6 0];
escenario8 = [8 40 120 3.5*1e9 0];
escenario9 = [9 0 0 0 1];

% interpolacion y distancia pilotos
metodos_interp = ["fft" "cubic" "perfect"];
distancias = [5 10 20];

% configuracion
conf = {metodos_interp distancias snrs nro_ejecuciones};

% resultados para 16QAM
qam_1 = simulate("16QAM", escenario1, conf); save("results/qam_1.mat", 'qam_1');
qam_2 = simulate("16QAM", escenario2, conf); save("results/qam_2.mat", 'qam_2');
qam_3 = simulate("16QAM", escenario3, conf); save("results/qam_3.mat", 'qam_3');
qam_4 = simulate("16QAM", escenario4, conf); save("results/qam_4.mat", 'qam_4');
qam_5 = simulate("16QAM", escenario5, conf); save("results/qam_5.mat", 'qam_5');
qam_6 = simulate("16QAM", escenario6, conf); save("results/qam_6.mat", 'qam_6');
qam_7 = simulate("16QAM", escenario7, conf); save("results/qam_7.mat", 'qam_7');
qam_8 = simulate("16QAM", escenario8, conf); save("results/qam_8.mat", 'qam_8');
qam_9 = simulate("16QAM", escenario9, conf); save("results/qam_9.mat", 'qam_9');

% resultados QPSK
qpsk_1 = simulate("QPSK", escenario1, conf); save("results/qpsk_1.mat", 'qpsk_1');
qpsk_2 = simulate("QPSK", escenario2, conf); save("results/qpsk_2.mat", 'qpsk_2');
qpsk_3 = simulate("QPSK", escenario3, conf); save("results/qpsk_3.mat", 'qpsk_3');
qpsk_4 = simulate("QPSK", escenario4, conf); save("results/qpsk_4.mat", 'qpsk_4');
qpsk_5 = simulate("QPSK", escenario5, conf); save("results/qpsk_5.mat", 'qpsk_5');
qpsk_6 = simulate("QPSK", escenario6, conf); save("results/qpsk_6.mat", 'qpsk_6');
qpsk_7 = simulate("QPSK", escenario7, conf); save("results/qpsk_7.mat", 'qpsk_7');
qpsk_8 = simulate("QPSK", escenario8, conf); save("results/qpsk_8.mat", 'qpsk_8');
qpsk_9 = simulate("QPSK", escenario9, conf); save("results/qpsk_9.mat", 'qpsk_9');

function BER = simulate(modulacion, escenario, conf)
    BER = {};
    metodos_interp = conf{1};
    distancias = conf{2};
    snrs = conf{3};
    nro_ejecuciones = conf{4}(1);
    i = escenario(1);
    disp("Working on: "+modulacion+"_"+string(i))
    for metodo = metodos_interp
        if metodo ~= "perfect"
            ND = length(distancias);
        else
            ND = 1;
        end
        for i=1:1:ND
            distancia = distancias(i);
            ber_list = [];
            for snr = snrs
                ber_snr = [];
                for n = 1:1:nro_ejecuciones
                    ber = get_ber(modulacion, escenario, metodo, distancia, snr, false);
                    ber_snr = [ber_snr ber];
                end
                ber_snr = mean(ber_snr);
                ber_list = [ber_list ber_snr];
            end
            BER{end+1} = ber_list;
        end
    end
end