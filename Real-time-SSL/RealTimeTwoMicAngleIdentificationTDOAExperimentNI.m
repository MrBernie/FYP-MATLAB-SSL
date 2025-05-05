%%2023-07-28
%%Two microphone time delay of arrival program
%%For angle identification

clc
clear all
close all

micDist = 0.1; % Distance between microphones (m)
numAngles = 180; % Number of angles to evaluate
speedOfSound = 343; % Speed of sound (m/s)


% Load the Data Acquisition Toolbox
if ~exist('daq','file')
    disp('Data Acquisition Toolbox is not installed or not available.');
    return;
end

% Create a DAQ session
s = daq.createSession('ni');

% Define acquisition parameters
fs = 51200; % Sampling rate (Hz)
duration = 0.2; % Duration of acquisition (s)

% Add an analog input channel for microphone
ch1 = s.addAnalogInputChannel('cDAQ1Mod1', 0,  'IEPE'); % Modify 'Dev1' to match your DAQ device ID and '0' to match the channel number of your microphone
ch2 = s.addAnalogInputChannel('cDAQ1Mod1', 1,  'IEPE'); % Modify 'Dev1' to match your DAQ device ID and '0' to match the channel number of your microphone

% Set acquisition parameters
s.Rate = fs;
s.DurationInSeconds = duration;




while true

    [data, time] = s.startForeground();




    % Get audio data
    %audioData = getaudiodata(recObj);
    channel1 = data(:,1);
    channel2 = data(:,2);


    channel1_fft = fft(channel1);
    channel2_fft = fft(channel2);

    N = size(channel1, 1); % Length of the audio data
    frequencies = (0:(N/2))*(fs/N); % Frequency vector


    xdft1=channel1_fft(1:N/2+1);
    xdft2=channel1_fft(1:N/2+1);
    psdx1 = (1/(fs*N)) * abs(xdft1).^2;
    psdx1(2:end-1) = 2*psdx1(2:end-1);
    psdx2 = (1/(fs*N)) * abs(xdft2).^2;
    psdx2(2:end-1) = 2*psdx2(2:end-1);
    Chanel_Sig1=psdx1;%1靠近absorber
    Chanel_Sig1_s=smooth(Chanel_Sig1,50);
    Chanel_Sig2=psdx2;%2靠近声源
    Chanel_Sig2_s=smooth(Chanel_Sig2,50);



    % Compute cross-correlation
    [correlation, lags] = xcorr(channel1, channel2);

    % Find max correlation and corresponding lag
    [maxCorr, idx] = max(correlation);
    lag = lags(idx);

    % Calculate TDOA (in seconds)
    tdoa = lag / fs;

    % Estimate direction
    angles = linspace(-90, 90, numAngles);
    tau = (micDist / speedOfSound) .* sind(angles);
    [val, angleIdx] = min(abs(tau - tdoa));
    sourceAngle = angles(angleIdx);

    % Calculate directivity pattern
    directivity = zeros(1, numAngles);
    for i = 1:numAngles
        directivity(i) = abs(sum(channel1 .* circshift(channel2, round(fs * tau(i)))));
    end
    

    % Plot directivity pattern
    %figure;
    %subplot(1,2,1)
    polarplot(deg2rad(angles), directivity/max(directivity), 'LineWidth', 2);
    hold on;
    polarplot(deg2rad(sourceAngle), directivity(angleIdx)/max(directivity), 'ro', 'MarkerSize', 10, 'LineWidth', 2);
    title(['Estimated Source Angle: ' num2str(sourceAngle) '°']);
    legend('Directivity Pattern', 'Estimated Source');
    set(gca, 'ThetaZeroLocation', 'top', 'ThetaDir', 'clockwise');
    hold off;

    drawnow;

    
end
