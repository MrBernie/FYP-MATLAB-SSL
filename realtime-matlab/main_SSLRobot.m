%2024-06-25
% Pi Rendong
% 4 Mic SSL robot
% Microphone array Realtime

clc
clear all
close all


%%%%%%%%%%%%System Setup%%%%%%%%%%%%
%%%%%%%%%%%%Theorectical Localization Result%%%%%%%


% Define the cylinder dimensions
cylinder_radius = 0.3;
cylinder_length = 3;

% Create the cylinder vertices
[cylinder_y, cylinder_z, cylinder_x] = cylinder(cylinder_radius, 20);
cylinder_x = cylinder_x * cylinder_length - 0.5;
cylinder_y = cylinder_y + 0;
cylinder_z = cylinder_z + 0;

% Define the sensor point coordinates
% 4 mics (x1, x2, x3, x4) (y1, y2, y3, y4), (z1, z2, z3, z4)
sensor_x = [0.11, 0.02, 0.02,   0.06];
sensor_y = [0.10, 0.11, 0.05, 0];
sensor_z = [0 0 0 0];

% Plot the cylinder and sensor points
figure;
plot3(cylinder_x, cylinder_y, cylinder_z, '-', 'LineWidth', 1,'color', [0.5 0.5 0.5]);
hold on;
plot3(sensor_x, sensor_y, sensor_z, 'r*', 'MarkerSize', 10);
axis equal;
xlabel('X');
ylabel('Y');
zlabel('Z');
title('Cylinder, Sensor and Source Locations');
grid on;
view(3); % Set the view to 3D


% Define the sound source position
sound_source_x =0.41;
sound_source_y = 0.07;
sound_source_z = 0;
plot3(sound_source_x, sound_source_y, sound_source_z, 'bo', 'MarkerSize', 10, 'MarkerFaceColor', 'b');


% Calculate the distances from the sound source to the receivers
distances = zeros(1, 4);
for i = 1:4
    distances(i) = sqrt((sound_source_x - sensor_x(i))^2 + (sound_source_y - sensor_y(i))^2 + (sound_source_z - sensor_z(i))^2);
end

% Calculate the arrival time differences
sound_speed = 343; % speed of sound in m/s
arrival_time = distances / sound_speed;
disp(['Arrival Time to Sensor 1 : (', num2str(arrival_time(1)*1000) ') ms']);
disp(['Arrival Time to Sensor 2 : (', num2str(arrival_time(2)*1000) ') ms']);
disp(['Arrival Time to Sensor 3 : (', num2str(arrival_time(3)*1000) ') ms']);
disp(['Arrival Time to Sensor 4 : (', num2str(arrival_time(4)*1000) ') ms']);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%2D plane Triangulation SSL%%%%%%%%%

%#1， 3 sensor triangulation localization in x and y plane, neglect z height
xs_true=sound_source_x;
ys_true=sound_source_y;

% Speed of sound
v = 343; % m/s, replace with the actual value

% Sensor positions
sensor1 = [sensor_x(1), sensor_y(1)];
sensor2 = [sensor_x(2), sensor_y(2)];
sensor3 = [sensor_x(3), sensor_y(3)];
% Sound source position
source = [xs_true, ys_true];

% Calculate distances
distance1 = norm(source - sensor1);
distance2 = norm(source - sensor2);
distance3 = norm(source - sensor3);

% Calculate time delays
delay1 = distance1 / v;
delay2 = distance2 / v;
delay3 = distance3 / v;

% Display time delays
% Nelget z-axis height difference
disp(['Time Delay between Sensor 1 and Sound Source: ', num2str(delay1*1000), ' milliseconds']);
disp(['Time Delay between Sensor 2 and Sound Source: ', num2str(delay2*1000), ' milliseconds']);
disp(['Time Delay between Sensor 3 and Sound Source: ', num2str(delay3*1000), ' milliseconds']);


disp(['Time Delay between Sensor 1 and Sensor 2: ', num2str((delay1-delay2)*1000), ' milliseconds']);
disp(['Time Delay between Sensor 2 and Sensor 3: ', num2str((delay2-delay3)*1000), ' milliseconds']);

min_samp_rate_req=1/min(delay1-delay2, delay1-delay2);
disp(['Minium sampling rate required ', num2str(min_samp_rate_req), 'Hz']);



% Known values
x1 = sensor_x(1); % Replace with the actual value
y1 = sensor_y(1); % Replace with the actual value
x2 = sensor_x(2); % Replace with the actual value
y2 = sensor_y(2); % Replace with the actual value
x3 = sensor_x(3); % Replace with the actual value
y3 = sensor_y(3); % Replace with the actual value

%z1=delay1-delay2;
%z2=delay2-delay3;

z1 = (sqrt((xs_true - x1)^2 + (ys_true - y1)^2) - sqrt((xs_true - x2)^2 + (ys_true - y2)^2))/v ;% Replace with the actual value
z2 = (sqrt((xs_true - x2)^2 + (ys_true - y2)^2) - sqrt((xs_true - x3)^2 + (ys_true - y3)^2))/v ;% Replace with the actual value

% Define the system of equations
eqn1 = @(xs, ys) sqrt((xs - x1)^2 + (ys - y1)^2) - sqrt((xs - x2)^2 + (ys - y2)^2) - z1 * v;
eqn2 = @(xs, ys) sqrt((xs - x2)^2 + (ys - y2)^2) - sqrt((xs - x3)^2 + (ys - y3)^2) - z2 * v;

% Solve the system of equations
equations = @(vars) [eqn1(vars(1), vars(2)); eqn2(vars(1), vars(2))];
initialGuess = [0; 0]; % Initial guess for (xs, ys)
estimatedCoordinates = fsolve(equations, initialGuess);

% Extract xs and ys
xs_3sensorXY = estimatedCoordinates(1);
ys_3sensorXY = estimatedCoordinates(2);

% Display the estimated coordinates
disp(['Method 1: 3 sensor triangulation localization in x and y plane, neglect z height']);
disp(['Estimated Coordinates (xs_3sensorXY, ys_3sensorXY): (', num2str(xs_3sensorXY), ', ', num2str(ys_3sensorXY), ')']);

disp(['Actual Coordinates (xs_true, ys_true): (', num2str(xs_true), ', ', num2str(ys_true), ')']);



%Method #2， 5 sensor triangulation localization in x and y plane, neglect z height

%%%Change to array sensors%%%
% Known values
x1 = sensor_x(1); % Replace with the actual value
y1 = sensor_y(1); % Replace with the actual value
x2 = sensor_x(2); % Replace with the actual value
y2 = sensor_y(2); % Replace with the actual value
x3 = sensor_x(3); % Replace with the actual value
y3 = sensor_y(3); % Replace with the actual value
x4 = sensor_x(4); % Replace with the actual value
y4 = sensor_y(4); % Replace with the actual value



z1 = (sqrt((xs_true - x1)^2 + (ys_true - y1)^2) - sqrt((xs_true - x2)^2 + (ys_true - y2)^2))/v; % Replace with the actual value
z2 = (sqrt((xs_true - x2)^2 + (ys_true - y2)^2) - sqrt((xs_true - x3)^2 + (ys_true - y3)^2))/v; % Replace with the actual value
z3 = (sqrt((xs_true - x3)^2 + (ys_true - y3)^2) - sqrt((xs_true - x4)^2 + (ys_true - y4)^2))/v; % Replace with the actual value


% Define the system of equations
equations = @(vars) [
    sqrt((vars(1) - x1)^2 + (vars(2) - y1)^2) - sqrt((vars(1) - x2)^2 + (vars(2) - y2)^2) - z1 * v;
    sqrt((vars(1) - x2)^2 + (vars(2) - y2)^2) - sqrt((vars(1) - x3)^2 + (vars(2) - y3)^2) - z2 * v;
    sqrt((vars(1) - x3)^2 + (vars(2) - y3)^2) - sqrt((vars(1) - x4)^2 + (vars(2) - y4)^2) - z3 * v;
    ];

% Solve the system of equations using least squares
initialGuess = [0; 0]; % Initial guess for (xs, ys)
estimatedCoordinates = lsqnonlin(equations, initialGuess);

% Extract xs and ys
xs_4sensorXY = estimatedCoordinates(1);
ys_4sensorXY = estimatedCoordinates(2);

% Display the estimated coordinates
disp(['Estimated Coordinates (xs_4sensorXY, ys_4sensorXY): (', num2str(xs_4sensorXY), ', ', num2str(ys_4sensorXY), ')']);
disp(['Actual Coordinates (xs_true, ys_true): (', num2str(xs_true), ', ', num2str(ys_true), ')']);

%Plot localization result
figure
plot3(cylinder_x, cylinder_y, cylinder_z, '-', 'LineWidth', 2,'color', [0.5 0.5 0.5]);
hold on;
plot3(sensor_x, sensor_y, sensor_z, 'r*', 'MarkerSize', 10);
axis equal;
xlabel('X');
ylabel('Y');
zlabel('Z');
title('Theorectical SSL result: Blue-Actual, Black-3 sensor, Red-4 sensor,');
grid on;
view(3); % Set the view to 3D

plot3(sound_source_x, sound_source_y, sound_source_z, 'bx', 'MarkerSize', 6, 'MarkerFaceColor', 'b');
plot3(xs_3sensorXY, ys_3sensorXY, 0 ,'ko', 'MarkerSize', 4, 'MarkerFaceColor', 'k');
plot3(xs_4sensorXY, ys_4sensorXY, 0 ,'ro', 'MarkerSize', 4, 'MarkerFaceColor', 'r');






%% - <> PARAMETERS <> -
% operating system
% 1 - macOS
% 2 - windows
inxOS = 2;

% switch real-time recording. [boolean]
booRT = 0;

% channel number
nMic = 4 ;
% bit-depth
nBit = 16;
% sample rate [Hz]
fs = 48000 ;
% time resolution [sec]
td = 1/fs;

% [NOTE] One second period is recommended.
% individual frame size in number
nFrame = 1;
szFrame = nFrame*fs;
% total frame size in seconds
tt = (szFrame-1)*td;

% [MONITOR]
fprintf('\n# PARAMETER #\n')
fprintf('> Sample frequency: %.2f kHz\n',fs/1e3)
fprintf('> Frame/Buffer size: 2^%i = %i\n',log2(szFrame),szFrame)
fprintf('> Sample length: %.2f sec\n',tt)

%% - <> RECORDING & SAVE DATA <> -
duration = 3; % [sec]
nSample = 1;

[aryData,colTime] = ...
    get_audioData_store(inxOS,fs,szFrame,nMic,nBit,nSample,duration);


%s1-s4 record time domain signals
s1 = aryData(:,1);
s2 = aryData(:,2);
s3 = aryData(:,3);%
s4 = aryData(:,4);%

p(1,:) = s1'; % Store sound pressure for ith microphone
p(2,:) = s2'; % Store sound pressure for ith microphone
p(3,:) = s3'; % Store sound pressure for ith microphone
p(4,:) = s4'; % Store sound pressure for ith microphone

t=colTime;
delta_t=colTime(2)-colTime(1);%time step
fs=1/delta_t %virtual sampling rate
N=length(t);
%frequencies = (0:(N/2))*(fs/N); % Frequency vector

N_channel=4; %Number of evaluation points

for ii = 1:N_channel
    Data_fft(:,ii) = abs(fft(aryData(:,ii)));
end


N_mic=4;%Number of mics
% Plot sound pressure signals
figure;
t_axis = t;
for ii = 1:N_mic
    subplot(N_mic, 1, ii);
    plot(t_axis, p(ii,:)); xlim([0 max(t)]);
    %    ii_pos=-(N_mic-1)*d_mic/2+(ii-1)*d_mic;%position of i-th mic
    %    title(sprintf('Microphone %d , Position %s m', ii, num2str(ii_pos)))
    ylabel('Pressure')
end
xlabel('Time (s)');
sgtitle('Sound pressure signals at microphone array');

figure;
% Plot frequency spectra
for i = 1:N_mic
    P_i = fft(p(i,:), [], 2); % Fourier transform of sound pressure signal
    f_axis = linspace(0, fs/2, length(P_i)/2+1); % Frequency axis
    P_i = P_i(1:length(P_i)/2+1); % Keep positive frequencies
    subplot(4, 1, i);
    plot(f_axis, 10*log10(abs(P_i))); xlim([0 fs/2]);

    %    ii_pos=-(N_mic-1)*d_mic/2+(i-1)*d_mic;%position of i-th mic
    %    title(sprintf('Microphone %d , Position %s m', i, num2str(ii_pos)))
    ylabel('Sound Pressure Level')

end
xlabel('Frequency (Hz)');
sgtitle('Frequency spectra from microphone array');


%%%%100-10000 Hz
figure;
% Plot frequency spectra
for i = 1:N_mic
    P_i = fft(p(i,:), [], 2); % Fourier transform of sound pressure signal
    f_axis = linspace(0, fs/2, length(P_i)/2+1); % Frequency axis
    P_i = P_i(1:length(P_i)/2+1); % Keep positive frequencies
    subplot(4, 1, i);
    plot(f_axis, 10*log10(abs(P_i))); xlim([100 10000]);

    %    ii_pos=-(N_mic-1)*d_mic/2+(i-1)*d_mic;%position of i-th mic
    %    title(sprintf('Microphone %d , Position %s m', i, num2str(ii_pos)))
    ylabel('Sound Pressure Level')

end
xlabel('Frequency (Hz)');
sgtitle('Frequency spectra from microphone array');





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%Real-time collected data%%%%%%
%%%%%%Mic Array Coordinate%%%%%%%%%%

% Get audio data
channel1 = s1;
channel2 = s2;
channel3 = s3;
channel4 = s4;


% Perform FFT analysis
nfft = N; % Number of points for FFT analysis
f = fs*(0:(nfft/2))/nfft; % Frequency vector


Y1 = fft(channel1, nfft); % FFT of microphone 1
Y2 = fft(channel2, nfft); % FFT of microphone 2
Y3 = fft(channel3, nfft); % FFT of microphone 3
Y4 = fft(channel4, nfft); % FFT of microphone 4

P1 = abs(Y1/nfft); % Single-sided amplitude spectrum of microphone 1
P2 = abs(Y2/nfft); % Single-sided amplitude spectrum of microphone 2
P3 = abs(Y3/nfft); % Single-sided amplitude spectrum of microphone 1
P4 = abs(Y4/nfft); % Single-sided amplitude spectrum of microphone 2



% Compute cross-correlation
[correlation1, lags1] = xcorr(channel1, channel2);
[correlation2, lags2] = xcorr(channel2, channel3);
[correlation3, lags3] = xcorr(channel3, channel4);

% Find max correlation and corresponding lag
[maxCorr1, idx1] = max(correlation1);
[maxCorr2, idx2] = max(correlation2);
[maxCorr3, idx3] = max(correlation3);
lag1 = lags1(idx1);
lag2 = lags2(idx2);
lag3 = lags3(idx3);

% Calculate TDOA (in seconds)
tdoa1 = lag1 / fs;
tdoa2 = lag2 / fs;
tdoa3 = lag3 / fs;

z1=tdoa1;
z2=tdoa2;
z3=tdoa3;





% Display Measured time delays
% Nelget z-axis height difference

disp(['Time Delay between Sensor 1 and Sensor 2: ', num2str((tdoa1)*1000), ' milliseconds']);
disp(['Time Delay between Sensor 2 and Sensor 3: ', num2str((tdoa2)*1000), ' milliseconds']);
disp(['Time Delay between Sensor 3 and Sensor 4: ', num2str((tdoa3)*1000), ' milliseconds']);



%%Three sensor, 2 TDOA%%%

% Define the system of equations
eqn1 = @(xs, ys) sqrt((xs - x1)^2 + (ys - y1)^2) - sqrt((xs - x2)^2 + (ys - y2)^2) - z1 * v;
eqn2 = @(xs, ys) sqrt((xs - x2)^2 + (ys - y2)^2) - sqrt((xs - x3)^2 + (ys - y3)^2) - z2 * v;

% Solve the system of equations
equations = @(vars) [eqn1(vars(1), vars(2)); eqn2(vars(1), vars(2))];
initialGuess = [0; 0]; % Initial guess for (xs, ys)
estimatedCoordinates = fsolve(equations, initialGuess);

% Extract xs and ys
xs_3sensorXY = estimatedCoordinates(1);
ys_3sensorXY = estimatedCoordinates(2);

% Display the estimated coordinates
disp(['Me fdsathod 1: 3 sensor triangulation localization in x and y plane, neglect z height']);
disp(['Experimental Estimated Coordinates (xs_3sensorXY, ys_3sensorXY): (', num2str(xs_3sensorXY), ', ', num2str(ys_3sensorXY), ')']);

disp(['Actual Coordinates (xs_true, ys_true): (', num2str(xs_true), ', ', num2str(ys_true), ')']);



%Method #2， 5 sensor triangulation localization in x and y plane, neglect z height


% Define the system of equations
equations = @(vars) [
    sqrt((vars(1) - x1)^2 + (vars(2) - y1)^2) - sqrt((vars(1) - x2)^2 + (vars(2) - y2)^2) - z1 * v;
    sqrt((vars(1) - x2)^2 + (vars(2) - y2)^2) - sqrt((vars(1) - x3)^2 + (vars(2) - y3)^2) - z2 * v;
    sqrt((vars(1) - x3)^2 + (vars(2) - y3)^2) - sqrt((vars(1) - x4)^2 + (vars(2) - y4)^2) - z3 * v;
    ];

% Solve the system of equations using least squares
initialGuess = [0; 0]; % Initial guess for (xs, ys)
estimatedCoordinates = lsqnonlin(equations, initialGuess);

% Extract xs and ys
xs_4sensorXY = estimatedCoordinates(1);
ys_4sensorXY = estimatedCoordinates(2);

% Display the estimated coordinates
disp(['Experimental Estimated Coordinates (xs_4sensorXY, ys_4sensorXY): (', num2str(xs_4sensorXY), ', ', num2str(ys_4sensorXY), ')']);
disp(['Actual Coordinates (xs_true, ys_true): (', num2str(xs_true), ', ', num2str(ys_true), ')']);

%Plot localization result
figure
plot3(cylinder_x, cylinder_y, cylinder_z, '-', 'LineWidth', 2,'color', [0.5 0.5 0.5]);
hold on;
plot3(sensor_x, sensor_y, sensor_z, 'r*', 'MarkerSize', 10);
axis equal;
xlabel('X');
ylabel('Y');
zlabel('Z');
title('Experimental SSL result: Blue-Actual, Black-3 sensor, Red-4 sensor,');
grid on;
view(3); % Set the view to 3D

plot3(sound_source_x, sound_source_y, sound_source_z, 'bx', 'MarkerSize', 6, 'MarkerFaceColor', 'b');
plot3(xs_3sensorXY, ys_3sensorXY, 0 ,'ko', 'MarkerSize', 4, 'MarkerFaceColor', 'k');
plot3(xs_4sensorXY, ys_4sensorXY, 0 ,'ro', 'MarkerSize', 4, 'MarkerFaceColor', 'r');








