
%%__main
%toma manual de datos

%% Define Measurement Parameters
positions = 0:10:100;          % Positions: 0, 10, 20, ..., 90
nPositions = numel(positions);
nMeasures = 21;               % Number of measurements per position

N = 281;

%% Initialize Struct to Store Measurements
Measurements = struct;

for i = 1:nPositions
    Measurements(i).Position      = positions(i);        % e.g., 0, 10, 20, ...
    Measurements(i).MeasurementID = (1:nMeasures)';      % Column vector [1..21]
    Measurements(i).BER           = zeros(nMeasures, 1); % Preallocate for BER
    Measurements(i).BER_noise       = zeros(nMeasures, 1); % Preallocate for BER_Noise
    Measurements(i).BER_PCSI      = zeros(nMeasures, 1); % Preallocate for BER_PCSI


    % CIRog is measured once per position (single vector).
    % Initialize to zeros for demonstration; replace with real data as needed.
    Measurements(i).CIRog         = zeros(nMeasures, N);          % 1 x N

    % CIRbpsk and CIRlab each have 21 measurements (same size as nMeasures).
    % One approach is to store them in a 2D numeric array: 21 rows (one per measurement),
    % each row is a CIR vector of length N. If your CIR lengths vary from measurement
    % to measurement, you might instead store them in cell arrays.
    Measurements(i).CIRbpsk       = zeros(nMeasures, N);  % 21 x N
    Measurements(i).CIRlab        = zeros(nMeasures, N);  % 21 x N
    Measurements(i).CIRnoH       = zeros(nMeasures, N);  % 21 x N
end

%% Loading the Measurements struct from the file
 loadedData = load('0507_MeasurementsDataLOS20dB.mat', 'Measurements');
 Measurements = loadedData.Measurements;


%% 
d = 100;
m = 14;
run("txscript.m")

%% 
run("rxscript.m")
run("measuresscript.m")

% %% Saving the Measurements struct to a file
save('0509_MeasurementsDataLOS0dB.mat', 'Measurements'); %cambiar por dB

d = d;
m = m+1;

if m >24
    d = d;
    m = m;
else
    run("txscript.m")
end
 


%% Loading the Measurements struct from the file
loadedData = load('0519_MeasurementsDataLOS20dB.mat', 'Measurements');
Measurements = loadedData.Measurements;

%% 

% Extract positions into a vector
positions = arrayfun(@(x) x.Position, Measurements);

% Create a logical index to exclude one or more specific distances
% Example: exclude distance = 0
excludedPositions = 0;  % you can also specify multiple, e.g., [0, 10, 20]
includeIdx = ~ismember(positions, excludedPositions);

% Subset only the Measurements you want
MeasurementsFiltered = Measurements(includeIdx);

for i = 1:10
    for j = 1:21
        if MeasurementsFiltered(i).BER(j) == 0
            MeasurementsFiltered(i).BER(j) = 10E-6;
        end

    end
end

% Recompute positions from filtered measurements
positionsFiltered = arrayfun(@(x) x.Position, MeasurementsFiltered);

% Compute mean BER for the filtered measurements
meanBERFiltered = arrayfun(@(x) mean(x.BER), MeasurementsFiltered);

% Plot the filtered data
figure;
semilogy(positionsFiltered, meanBERFiltered, 'o-', 'LineWidth', 1.5, 'MarkerSize', 6);
grid on;
xlabel('Distance (cm)');
ylabel('Mean BER');
title('Mean BER vs. Distance');

% Create a logical index to exclude one or more specific distances
% Example: exclude distance = 0
excludedPositions = 0;  % you can also specify multiple, e.g., [0, 10, 20]
includeIdx = ~ismember(positions, excludedPositions);

% Subset only the Measurements you want
MeasurementsFiltered = Measurements(includeIdx);

% Recompute positions from filtered measurements
positionsFiltered = arrayfun(@(x) x.Position, MeasurementsFiltered);

for i = 1:10
    for j = 1:21
        if MeasurementsFiltered(i).BER_PCSI(j) == 0
            MeasurementsFiltered(i).BER_PCSI(j) = 10E-6;
        end

    end
end

% Compute mean BER for the filtered measurements
meanBERFiltered = arrayfun(@(x) mean(x.BER_PCSI), MeasurementsFiltered);

% Plot the filtered data
hold on 
semilogy(positionsFiltered, meanBERFiltered, 'diamond-', 'LineWidth', 1.5, 'MarkerSize', 6);

% Create a logical index to exclude one or more specific distances
% Example: exclude distance = 0
excludedPositions = 0;  % you can also specify multiple, e.g., [0, 10, 20]
includeIdx = ~ismember(positions, excludedPositions);

% Subset only the Measurements you want
MeasurementsFiltered = Measurements(includeIdx);

for i = 1:10
    for j = 1:21
        if MeasurementsFiltered(i).BER_noise(j) == 0
            MeasurementsFiltered(i).BER_noise(j) = 10E-6;
        end

    end
end

% Recompute positions from filtered measurements
positionsFiltered = arrayfun(@(x) x.Position, MeasurementsFiltered);

% Compute mean BER for the filtered measurements
meanBERFiltered = arrayfun(@(x) mean(x.BER_noise), MeasurementsFiltered);

% Plot the filtered data
hold on 
semilogy(positionsFiltered, meanBERFiltered, 'square-', 'LineWidth', 1.5, 'MarkerSize', 6);
legend( 'BPSK Pilot', 'Perfect CSI', 'Pseudo-random Noise')


%% 


% Plot the filtered data
figure
semilogy(positionsFiltered, meanBERFiltered, 'o-', 'LineWidth', 1.5, 'MarkerSize', 6);
grid on;
xlabel('Distance');
ylabel('Mean BER');
title('Mean BER (Imperfect CSI) vs. Distance');
hold on
errorbar(positionsFiltered, meanBERFiltered, stdBER, 'LineStyle', 'none', 'Color', 'blue');
set(gca, 'YScale', 'log');


%% 

% Plot the filtered data
figure
semilogy(positionsFiltered, meanBERFiltered_PCSI, 'diamond-', 'LineWidth', 1.5, 'MarkerSize', 6);
grid on;
xlabel('Distance');
ylabel('Mean BER');
title('Mean BER (Perfect CSI) vs. Distance');
hold on
errorbar(positionsFiltered, meanBERFiltered_PCSI, stdBER_PCSI, 'LineStyle', 'none', 'Color', 'blue');
set(gca, 'YScale', 'log');

% legend('Imperfect CSI', 'Perfect CSI')

%% 

%% Assume 'Measurements' struct already in the workspace
%  Measurements(i).CIRbpsk is 21 x N for each position i

nPositions = numel(MeasurementsFiltered);
N = size(MeasurementsFiltered(1).CIRbpsk, 2);  % Length of each CIR vector
MeanCIRbpsk = zeros(nPositions, N);    % Preallocate for mean vectors

for i = 1:nPositions
    % Compute mean across the 21 rows (measurements) for position i
    % This yields a 1xN row vector
    MeanCIRbpsk(i, :) = mean(MeasurementsFiltered(i).CIRbpsk, 1);
end

% Now 'MeanCIRbpsk' is an nPositions x N array
%  where each row corresponds to the mean CIR vector for that position.

% (Optional) If you also have Measurements(i).CIRlab, do the same:
MeanCIRlab = zeros(nPositions, N);
for i = 1:nPositions
    MeanCIRlab(i, :) = mean(Measurements(i).CIRlab, 1);
end

MeanCIRog = zeros(nPositions, N);
for i = 1:nPositions
    MeanCIRog(i, :) = mean(Measurements(i).CIRog, 1);
end

i = 5;
figure
plot([1:281],MeanCIRog(i,:), '*-', 'LineWidth', 1.5, 'MarkerIndices', 1:5:length(MeanCIRog));
hold on
plot([1:281],MeanCIRlab(i,:), 'square-', 'LineWidth', 1.5, 'MarkerIndices', 1:5:length(MeanCIRog));
plot([1:281],MeanCIRbpsk(i,:), 'diamond-', 'LineWidth', 1.5, 'MarkerIndices', 1:5:length(MeanCIRog));
hold off
grid on;
title(['CIR Plot for Distance ',num2str(i*10),' cm']);
legend('Original CIR','Noise Recovered CIR','BPSK symbol Recovered CIR')

%% 



%% Assume 'Measurements' struct already in the workspace
%  Measurements(i).CIRbpsk is 21 x N for each position i

nPositions = numel(MeasurementsFiltered);
N = size(MeasurementsFiltered(1).CIRnoH, 2);  % Length of each CIR vector
MeanCIRnoH = zeros(nPositions, N);    % Preallocate for mean vectors

for i = 1:nPositions
    % Compute mean across the 21 rows (measurements) for position i
    % This yields a 1xN row vector
    MeanCIRnoH(i, :) = mean(MeasurementsFiltered(i).CIRnoH, 1);
end

i = 4;
figure
plot([1:281],MeanCIRnoH(i,:), '*-', 'LineWidth', 1.5, 'MarkerIndices', 1:5:length(MeanCIRog));
grid on;
title(['Laboratory CIR Plot for Distance ',num2str(i*10),' cm']);
%legend('Original CIR','Noise Recovered CIR','BPSK symbol Recovered CIR')





