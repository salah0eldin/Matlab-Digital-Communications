clear;     % Clear all variables from the workspace
close all; % Close all open figures

%% Control Flags
N = 7;                          % Number of repeatetion of bits
A = 4;                          % Define the amplitude level
num_of_bits = 100;              % Number of bits per waveform
num_of_realizations = 500;      % Number of random realizations (ensamble size)
time_samples = N * num_of_bits; % Total number of samples in the waveform

%% Initialize Matrices to Store Signal Realizations
% Creating zero matrices to store the generated waveforms for each line coding scheme
unipolar_ensamble  = zeros(num_of_realizations, time_samples);
NRZ_polar_ensamble = zeros(num_of_realizations, time_samples);
RZ_polar_ensamble  = zeros(num_of_realizations, time_samples);

%% Generate ensambles with random delay for different line coding schemes
for i = 1:num_of_realizations
    
    x = randi([0 1], 1, num_of_bits + 1); % Generate a random binary sequence of 101 bits
    x = repelem(x, N);                    % Repeat each bit 7 times to simulate pulse width of 70ms
    delay = randi([0 N-1]);               % Generate a random delay between 0 and 6 samples
    
    % Unipolar Encoding: 0 -> 0, 1 -> A, and Take 700 samples starting from the delay
    unipolar_ensamble(i, :) = x(N + 1 - delay : time_samples + N - delay) * A;
    
    % NRZ Polar Encoding: 0 -> -A, 1 -> A
    NRZ_polar_ensamble(i, :) = x(N + 1 - delay : time_samples + N - delay) * 2 * A - A * ones(1, time_samples);
    
    % RZ Polar Encoding: Modify signal to introduce return-to-zero behavior
    for j = N:N:(time_samples + N)
        x(j)     = 0.5; % Set last sample of each bit to 0.5
        x(j - 1) = 0.5; % Set second last sample of each bit to 0.5
        x(j - 2) = 0.5; % Set third last sample of each bit to 0.5
    end
    
    RZ_polar_ensamble(i, :) = x(N + 1 - delay : time_samples + N - delay) * 2 * A - A * ones(1, time_samples);
end

%% Plot the ensamble realizations for all three line coding schemes
figure;

% Plot Unipolar realization
subplot(3, 1, 1);
plot(unipolar_ensamble(num_of_realizations - 1, :), 'LineWidth', 0.5, 'Color', [0.2, 0.4, 0.6]);
grid on;
xlabel('Time');
ylabel('Value');
title('Unipolar Ensamble');
ylim([-A - 1, A + 1]); % Set Y-axis limits

% Plot NRZ Polar realization
subplot(3, 1, 2);
plot(NRZ_polar_ensamble(num_of_realizations - 1, :), 'LineWidth', 0.5, 'Color', [0.6, 0.4, 0.2]);
grid on;
xlabel('Time');
ylabel('Value');
title('NRZ Polar Ensamble');
ylim([-A - 1, A + 1]);

% Plot RZ Polar realization
subplot(3, 1, 3);
plot(RZ_polar_ensamble(num_of_realizations - 1, :), 'LineWidth', 0.5, 'Color', [0.4, 0.6, 0.2]);
grid on;
xlabel('Time');
ylabel('Value');
title('RZ Polar Ensamble');
ylim([-A - 1, A + 1]);

%% 1- Compute Statistical Mean for each line coding scheme
unipolar_statistical_mean  = zeros(1, time_samples); % Initialize vector for mean values
NRZ_polar_statistical_mean = zeros(1, time_samples);
RZ_polar_statistical_mean  = zeros(1, time_samples);

% Compute mean by summing all realizations and dividing by num_of_realizations
for i = 1:num_of_realizations
    unipolar_statistical_mean  = unipolar_ensamble(i, :)  + unipolar_statistical_mean;
    NRZ_polar_statistical_mean = NRZ_polar_ensamble(i, :) + NRZ_polar_statistical_mean;
    RZ_polar_statistical_mean  = RZ_polar_ensamble(i, :)  + RZ_polar_statistical_mean;
end

unipolar_statistical_mean  = unipolar_statistical_mean  / num_of_realizations;
NRZ_polar_statistical_mean = NRZ_polar_statistical_mean / num_of_realizations;
RZ_polar_statistical_mean  = RZ_polar_statistical_mean  / num_of_realizations;

%% Plot Statistical Mean for all three line coding schemes
figure;

% Unipolar Statistical Mean
subplot(3, 1, 1);
plot(unipolar_statistical_mean, 'LineWidth', 1, 'Color', [0.6, 0.2, 0.2]);
hold on;
yline(A / 2, '--', 'LineWidth', 1, 'Color', [0.2, 0.2, 0.2]); % Theoretical mean line at y = A/2
hold off;
grid on;
xlabel('Time');
ylabel('Value');
title('Statistical Mean of Unipolar Ensamble');
ylim([-A - 1, A + 1]);

% NRZ Polar Statistical Mean
subplot(3, 1, 2);
plot(NRZ_polar_statistical_mean, 'LineWidth', 1, 'Color', [0.2, 0.6, 0.2]);
hold on;
yline(0, '--', 'LineWidth', 1, 'Color', [0.2, 0.2, 0.2]); % Theoretical mean line at y = 0
hold off;
grid on;
xlabel('Time');
ylabel('Value');
title('Statistical Mean of NRZ Polar Ensamble');
ylim([-2 * A - 1, 2 * A + 1]);

% RZ Polar Statistical Mean
subplot(3, 1, 3);
plot(RZ_polar_statistical_mean, 'LineWidth', 1, 'Color', [0.2, 0.2, 0.6]);
hold on;
yline(0, '--', 'LineWidth', 1, 'Color', [0.2, 0.2, 0.2]); % Theoretical mean line at y = 0
hold off;
grid on;
xlabel('Time');
ylabel('Value');
title('Statistical Mean of RZ Polar Ensamble');
ylim([-2 * A - 1, 2 * A + 1]);

%% 2- Compute Statistical Autocorrelation
% Initialize autocorrelation matrices
unipolar_statistical_autocorrelation  = zeros(time_samples / 2, time_samples / 2);
nrz_polar_statistical_autocorrelation = zeros(time_samples / 2, time_samples / 2);
rz_polar_statistical_autocorrelation  = zeros(time_samples / 2, time_samples / 2);

% Compute autocorrelation for different time delays
for time = 1:(time_samples / 2)
    for time_diff = 0:(time_samples / 2 - 1)
        unipolar_statistical_autocorrelation(time,  time_diff + 1) = dot(unipolar_ensamble(:, time),  unipolar_ensamble(:,  time + time_diff));
        nrz_polar_statistical_autocorrelation(time, time_diff + 1) = dot(NRZ_polar_ensamble(:, time), NRZ_polar_ensamble(:, time + time_diff));
        rz_polar_statistical_autocorrelation(time,  time_diff + 1) = dot(RZ_polar_ensamble(:, time),  RZ_polar_ensamble(:,  time + time_diff));
    end
end

% Normalize autocorrelation by the number of realizations
unipolar_statistical_autocorrelation  = unipolar_statistical_autocorrelation  / num_of_realizations;
nrz_polar_statistical_autocorrelation = nrz_polar_statistical_autocorrelation / num_of_realizations;
rz_polar_statistical_autocorrelation  = rz_polar_statistical_autocorrelation  / num_of_realizations;

% Mirror the second dimension in negative
unipolar_statistical_autocorrelation_conc  = [fliplr(unipolar_statistical_autocorrelation(:, 2:end)),  unipolar_statistical_autocorrelation];
nrz_polar_statistical_autocorrelation_conc = [fliplr(nrz_polar_statistical_autocorrelation(:, 2:end)), nrz_polar_statistical_autocorrelation];
rz_polar_statistical_autocorrelation_conc  = [fliplr(rz_polar_statistical_autocorrelation(:, 2:end)),  rz_polar_statistical_autocorrelation];

% Calculate autocorrelation at avg time for all tau
unipolar_autocorrelation_avg  = sum(unipolar_statistical_autocorrelation_conc, 1)  / (time_samples / 2);
nrz_polar_autocorrelation_avg = sum(nrz_polar_statistical_autocorrelation_conc, 1) / (time_samples / 2);
rz_polar_autocorrelation_avg  = sum(rz_polar_statistical_autocorrelation_conc, 1)  / (time_samples / 2);

% Theoretical autocorrelations
unipolar_autocorrelation_theoretical  = A^2 / 4 * ones(1, time_samples/2);
nrz_polar_autocorrelation_theoretical = zeros(1, time_samples/2);
rz_polar_autocorrelation_theoretical  = zeros(1, time_samples/2);
for i = 1:N
    unipolar_autocorrelation_theoretical(i) = A^2 / 2 * (1 - (i - 1) / (2*N));
end
for i = 1:N
    nrz_polar_autocorrelation_theoretical(i) = A^2 * (1 - (i - 1) / N);
end
for i = 1:N-3
    rz_polar_autocorrelation_theoretical(i) = 4 * A^2 / 7 * (1 - 7 * (i - 1) / (4*N));
end

% Mirror the second dimension in negative
unipolar_autocorrelation_theoretical_conc  = [fliplr(unipolar_autocorrelation_theoretical(2:end)),  unipolar_autocorrelation_theoretical];
nrz_polar_autocorrelation_theoretical_conc = [fliplr(nrz_polar_autocorrelation_theoretical(2:end)), nrz_polar_autocorrelation_theoretical];
rz_polar_autocorrelation_theoretical_conc  = [fliplr(rz_polar_autocorrelation_theoretical(2:end)),  rz_polar_autocorrelation_theoretical];

%% 3D Plot of Statistical Autocorrelation for all three line coding schemes

% Create meshgrid for 3D plot
% [X, Y] = meshgrid(-time_samples+1:time_samples-1, 1:time_samples);
[X, Y] = meshgrid(-time_samples/2+1:time_samples/2-1, 1:time_samples/2);

figure;

% Unipolar Autocorrelation
subplot(3, 1, 1);
surf(X, Y, unipolar_statistical_autocorrelation_conc, 'EdgeColor', 'none');
colorbar;
xlabel('Time Diff');
ylabel('Time');
zlabel('Autocorrelation');
title('3D Plot of Unipolar Autocorrelation');
view(3);

% NRZ Polar Autocorrelation
subplot(3, 1, 2);
surf(X, Y, nrz_polar_statistical_autocorrelation_conc, 'EdgeColor', 'none');
colorbar;
xlabel('Time Diff');
ylabel('Time');
zlabel('Autocorrelation');
title('3D Plot of NRZ Polar Autocorrelation');
view(3);

% RZ Polar Autocorrelation
subplot(3, 1, 3);
surf(X, Y, rz_polar_statistical_autocorrelation_conc, 'EdgeColor', 'none');
colorbar;
xlabel('Time Diff');
ylabel('Time');
zlabel('Autocorrelation');
title('3D Plot of RZ Polar Autocorrelation');
view(3);

%% Plot Statistical Autocorrelation for all three line coding schemes
figure;
subplot(3, 1, 1);
plot(X, unipolar_autocorrelation_avg, 'LineWidth', 1, 'Color', [0.2, 0.4, 0.6]);
hold on;
plot(X, unipolar_autocorrelation_theoretical_conc, '--', 'LineWidth', 1, 'Color', [0.8, 0.2, 0.2]); 
hold off;
grid on;
xlabel('Time Diff'); ylabel('Value');
title('Statistical Autocorrelation of Unipolar Ensamble');
ylim([A/2, A^2/2 + 1]);
legend('Statistical', 'Theoretical');

subplot(3, 1, 2);
plot(X, nrz_polar_autocorrelation_avg, 'LineWidth', 1, 'Color', [0.6, 0.4, 0.2]);
hold on;
plot(X, nrz_polar_autocorrelation_theoretical_conc, '--', 'LineWidth', 1, 'Color', [0.8, 0.2, 0.2]); 
hold off;
grid on;
xlabel('Time Diff'); ylabel('Value');
title('Statistical Autocorrelation of NRZ Polar Ensamble');
ylim([-2, A^2 + 4]);
legend('Statistical', 'Theoretical');

subplot(3, 1, 3);
plot(X, rz_polar_autocorrelation_avg, 'LineWidth', 1, 'Color', [0.2, 0.6, 0.4]);
hold on;
plot(X, rz_polar_autocorrelation_theoretical_conc, '--', 'LineWidth', 1, 'Color', [0.8, 0.2, 0.2]); 
hold off;
grid on;
xlabel('Time Diff'); ylabel('Value');
title('Statistical Autocorrelation of RZ Polar Ensamble');
ylim([-2, A^2 / 2 + 4]);
legend('Statistical', 'Theoretical');

