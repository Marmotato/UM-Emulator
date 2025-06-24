

%% Example: Filling in the Measurements
% In a real scenario, you'd replace this loop with actual data collection
% or file reading. Below is just an example of assigning random values.
i = d/10 +1;% POSICION ACTUAL
m = m;% MEDICION ACTUAL

fprintf('Salvando medidas para Distancia = %d, Medicion #%d\n', d, m);

% Assign some dummy values (replace with your measured values)
Measurements(i).BER(m)      = BER;  % e.g., a random BER
Measurements(i).BER_noise(m)  = BER_noise;  %
Measurements(i).BER_PCSI(m) = BER_PCSI;  % e.g., a random BER_PCSI

Measurements(i).CIRog(m,:) = Final_response;
Measurements(i).CIRbpsk(m, :) = extractedCIRNorm';
Measurements(i).CIRlab(m, :)  = extractedCIRNorm_noise';
Measurements(i).CIRnoH       = extractedCIRNorm_noise_noH;  % 21 x N


%% Compute Mean and Standard Deviation for Each Position
meanBER      = zeros(nPositions, 1);
stdBER       = zeros(nPositions, 1);
meanBER_PCSI = zeros(nPositions, 1);
stdBER_PCSI  = zeros(nPositions, 1);

for i = 1:nPositions
    meanBER(i)      = mean(Measurements(i).BER);
    stdBER(i)       = std(Measurements(i).BER);
    meanBER_PCSI(i) = mean(Measurements(i).BER_PCSI);
    stdBER_PCSI(i)  = std(Measurements(i).BER_PCSI);
end

fprintf('=======================================================%d\n', d, m);
fprintf('Guardado para Distancia = %d, Medicion #%d\n', d, m);
fprintf('=======================================================%d\n', d, m);


%% Display or Plot Results
% For example, display results in the Command Window:
% disp(table(positions', meanBER, stdBER, meanBER_PCSI, stdBER_PCSI, ...
%     'VariableNames', {'Position','Mean_BER','Std_BER','Mean_BER_PCSI','Std_BER_PCSI'}));

% Alternatively, you can create plots or error bars:
% figure;
% errorbar(positions, meanBER, stdBER, 'o-');
% hold on;
% errorbar(positions, meanBER_PCSI, stdBER_PCSI, 'x-');
% legend('BER','BER\_PCSI','Location','Best');
% xlabel('Position');
% ylabel('Measurement Value');
% title('Mean and Standard Deviation of BER & BER\_PCSI by Position');
% grid on;


