%% Close and Clear
close all;
clear;

%% Define the bit stream

% Number of waveforms and bits
WAVEFORMS = 500;
BITS = 100;
EXTEND_VALUE = 7;

% Generate random bit stream
bitStream = randi([0, 1], WAVEFORMS, BITS);

% Extend each bit 7 times
extendedBitStream = repelem(bitStream, 1, EXTEND_VALUE);

% Voltage Level
VOLTAGE_LEVEL = 4;

% Unipolar encoding
unipolar = extendedBitStream * VOLTAGE_LEVEL;

% Polar Non Return to Zero (NRZ) encoding
polarNRZ = (2 * extendedBitStream - 1) * VOLTAGE_LEVEL;

% Polar Return to Zero (RZ) 4-3 encoding
polarRZ = polarNRZ; % Copy the polar NRZ signal
polarRZ(:, 5:7:end) = 0; % Clears the 5th bit of each 7 bits
polarRZ(:, 6:7:end) = 0; % Clears the 6th bit of each 7 bits
polarRZ(:, 7:7:end) = 0; % Clears the 7th bit of each 7 bits

%% Plot the signals

figure('Position', [100, 100, 1200, 800]);
plotSignal(5, 1, bitStream(1, 1:11), 'Original Bit Stream (1, 1:11)', [-2 2], 1);
plotSignal(5, 1, extendedBitStream(1, 1:71), 'Extended Bit Stream (1, 1:11)', [-2 2], 2);
plotSignal(5, 1, unipolar(1, 1:71), 'Unipolar Signal (1, 1:11)', [-6 6], 3);
plotSignal(5, 1, polarNRZ(1, 1:71), 'Polar NRZ Signal (1, 1:11)', [-6 6], 4);
plotSignal(5, 1, polarRZ(1, 1:71), 'Polar RZ Signal (1, 1:11)', [-6 6], 5);

%% Apply random shifts to the waveforms of each encoding

shiftedUnipolar = applyRandomShifts(unipolar);
shiftedPolarNRZ = applyRandomShifts(polarNRZ);
shiftedPolarRZ = applyRandomShifts(polarRZ);

%% Plot the shifted signals

figure('Position', [100, 100, 1200, 800]);
plotSignal(4, 1, unipolar(1, 1:71), 'Unipolar Signal (1, 1:11)', [-6 6], 1);
plotSignal(4, 1, shiftedUnipolar(1, 1:71), 'Shifted Unipolar Signal (1, 1:11)', [-6 6], 2);
plotSignal(4, 1, unipolar(2, 1:71), 'Unipolar Signal (2, 1:11)', [-6 6], 3);
plotSignal(4, 1, shiftedUnipolar(2, 1:71), 'Shifted Unipolar Signal (2, 1:11)', [-6 6], 4);

figure('Position', [100, 100, 1200, 800]);
plotSignal(4, 1, polarNRZ(1, 1:71), 'PolarNRZ Signal (1, 1:11)', [-6 6], 1);
plotSignal(4, 1, shiftedPolarNRZ(1, 1:71), 'Shifted PolarNRZ Signal (1, 1:11)', [-6 6], 2);
plotSignal(4, 1, polarNRZ(2, 1:71), 'PolarNRZ Signal (2, 1:11)', [-6 6], 3);
plotSignal(4, 1, shiftedPolarNRZ(2, 1:71), 'Shifted PolarNRZ Signal (2, 1:11)', [-6 6], 4);

figure('Position', [100, 100, 1200, 800]);
plotSignal(4, 1, polarRZ(1, 1:71), 'PolarRZ Signal (1, 1:11)', [-6 6], 1);
plotSignal(4, 1, shiftedPolarRZ(1, 1:71), 'Shifted PolarRZ Signal (1, 1:11)', [-6 6], 2);
plotSignal(4, 1, polarRZ(2, 1:71), 'PolarRZ Signal (2, 1:11)', [-6 6], 3);
plotSignal(4, 1, shiftedPolarRZ(2, 1:71), 'Shifted PolarRZ Signal (2, 1:11)', [-6 6], 4);

%% Calculate Outputs

outputUnipolar = shiftedUnipolar;
outputPolarNRZ = shiftedPolarNRZ;
outputPolarRZ = shiftedPolarRZ;

%% Statistical mean

% Calculate the probability of each output for UniPolar encoding cross all waveforms
meanUnipolar = sum(outputUnipolar) / WAVEFORMS;
meanPolarNRZ = sum(outputPolarNRZ) / WAVEFORMS;
meanPolarRZ = sum(outputPolarRZ) / WAVEFORMS;

% Plot the mean of the Unipolar signal
figure('Position', [100, 100, 1200, 800]);
plotSignal(3, 1, meanUnipolar, 'Mean of Unipolar Signal', [-1 5], 1, false);
plotSignal(3, 1, meanPolarNRZ, 'Mean of Polar NRZ Signal', [-5 5], 2, false);
plotSignal(3, 1, meanPolarRZ, 'Mean of Polar RZ Signal', [-5 5], 3, false);

%% Stationarity

%% Functions

% Function to apply random shifts to the bit stream
function shiftedBitStream = applyRandomShifts(extendedBitStream)
    WAVEFORMS = size(extendedBitStream, 1);
    shiftAmounts = randi([0, 6], WAVEFORMS, 1);
    shiftedBitStream = zeros(size(extendedBitStream));

    for i = 1:WAVEFORMS
        shiftedBitStream(i, :) = circshift(extendedBitStream(i, :), -shiftAmounts(i), 2);
    end

end

% Function to plot signals
function plotSignal(n, m, signal, titleText, yLimits, subplotIndex, xticksFlag)
    if nargin < 7
        xticksFlag = true;
    end
    
    subplot(n, m, subplotIndex);
    stairs(signal, 'LineWidth', 1.5);
    title(titleText);
    ylim(yLimits);
    grid on;
    if xticksFlag
        xticks(0:1:length(signal));
    end
    xlim([1 length(signal)]);
    yline(0, '--');
end
