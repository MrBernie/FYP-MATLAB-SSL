%Pi Rendong
% 4 Mic SSL robot
clc
clear all
close all

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
sensor_x = [0, 0.01, 0.02, 0.03];
sensor_y = [0, 0.05, 0.1, 0.15];
sensor_z = [0 0 0 0];

% Plot the cylinder and sensor points
figure;
plot3(cylinder_x, cylinder_y, cylinder_z, '-', 'LineWidth', 2,'color', [0.5 0.5 0.5]);
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
sound_source_x = 1;
sound_source_y = 0;
sound_source_z = 0.01;
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


%%Three sensor, 2 TDOA%%%

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
x5 = 0; % Replace with the actual value  - Put one more sensor in the middile
y5 = 0; % Replace with the actual value  - Put one more sensor in the middile


z1 = (sqrt((xs_true - x1)^2 + (ys_true - y1)^2) - sqrt((xs_true - x2)^2 + (ys_true - y2)^2))/v; % Replace with the actual value
z2 = (sqrt((xs_true - x2)^2 + (ys_true - y2)^2) - sqrt((xs_true - x3)^2 + (ys_true - y3)^2))/v; % Replace with the actual value
z3 = (sqrt((xs_true - x3)^2 + (ys_true - y3)^2) - sqrt((xs_true - x4)^2 + (ys_true - y4)^2))/v; % Replace with the actual value
z4 = (sqrt((xs_true - x4)^2 + (ys_true - y4)^2) - sqrt((xs_true - x5)^2 + (ys_true - y5)^2))/v; % Replace with the actual value


% Define the system of equations
equations = @(vars) [
    sqrt((vars(1) - x1)^2 + (vars(2) - y1)^2) - sqrt((vars(1) - x2)^2 + (vars(2) - y2)^2) - z1 * v;
    sqrt((vars(1) - x2)^2 + (vars(2) - y2)^2) - sqrt((vars(1) - x3)^2 + (vars(2) - y3)^2) - z2 * v;
    sqrt((vars(1) - x3)^2 + (vars(2) - y3)^2) - sqrt((vars(1) - x4)^2 + (vars(2) - y4)^2) - z3 * v;
    sqrt((vars(1) - x4)^2 + (vars(2) - y4)^2) - sqrt((vars(1) - x5)^2 + (vars(2) - y5)^2) - z4 * v;
    ];

% Solve the system of equations using least squares
initialGuess = [0; 0]; % Initial guess for (xs, ys)
estimatedCoordinates = lsqnonlin(equations, initialGuess);

% Extract xs and ys
xs_5sensorXY = estimatedCoordinates(1);
ys_5sensorXY = estimatedCoordinates(2);

% Display the estimated coordinates
disp(['Estimated Coordinates (xs_5sensorXY, ys_5sensorXY): (', num2str(xs_5sensorXY), ', ', num2str(ys_5sensorXY), ')']);
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
title('SSL result: Blue-Actual, Black-3 sensor, Red-5 sensor,');
grid on;
view(3); % Set the view to 3D

plot3(sound_source_x, sound_source_y, sound_source_z, 'bx', 'MarkerSize', 6, 'MarkerFaceColor', 'b');
plot3(xs_3sensorXY, ys_3sensorXY, 0 ,'ko', 'MarkerSize', 4, 'MarkerFaceColor', 'k');
plot3(xs_5sensorXY, ys_5sensorXY, 0 ,'ro', 'MarkerSize', 4, 'MarkerFaceColor', 'r');









