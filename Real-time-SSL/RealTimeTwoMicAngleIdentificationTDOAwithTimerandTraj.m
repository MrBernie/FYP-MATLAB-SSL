%%2024-02-15
%%Two microphone time delay of arrival program
%%For angle identification
%%Real Time-Audio Interface with two cheap microphones

clc
clear all
close all

% Parameters
fs = 44100; % Sample rate (Hz)
maxTime=150; % Maximum recording time
%duration = 1; % Recording duration (s)
micDist = 0.16; % Distance between microphones (m)
numAngles = 180; % Number of angles to evaluate
speedOfSound = 343; % Speed of sound (m/s)

N = 1024; % FFT length

win = hann(N); % Hanning window
win = ones(N,1); % Hanning window


% % Record audio signals
% recObj = audiorecorder(fs, 16, 2);

% Set up the audio device reader
reader = audioDeviceReader('SampleRate', fs, 'NumChannels', 2);

% Set up the frequency vector
f = (0:N/2-1)/N*fs;

t0=tic;

icount=1;
while toc(t0)<maxTime

    % Read a frame from the audio device
    x = reader();

    % Apply the window function
    x = x .* win;

    % Get audio data
    %audioData = getaudiodata(recObj);
    channel1 = x(:, 1);
    channel2 = x(:, 2);


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





    subplot(2, 2, 1);%Plot Mic 1 FFT data and smoothed
    plot(frequencies, 10*log10(Chanel_Sig1));
    %plot(frequencies,10*log10(Chanel_Sig1_s),'r')
    legend('Frequency Spectrum')
    xlabel('Frequency (Hz)');ylabel('Channel 1');
    title(['Channel 1 Frequency Domain']);
    axis([100 10000 -160 -80])
    grid on

    subplot(2, 2, 3);%Plot Mic 1 FFT data and smoothed
    plot(frequencies, 10*log10(Chanel_Sig2));
    %plot(frequencies,10*log10(Chanel_Sig2_s),'r')
    legend('Frequency Spectrum')
    xlabel('Frequency (Hz)');ylabel('Channel 2');
    title(['Channel 2 Frequency Domain']);
    axis([100 10000 -160 -80])
    grid on


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
    subplot(2, 2, 2);

    polarplot(deg2rad(angles), directivity/max(directivity), 'LineWidth', 2);
    hold on;
    polarplot(deg2rad(sourceAngle), directivity(angleIdx)/max(directivity), 'ro', 'MarkerSize', 10, 'LineWidth', 2);
    title(['Estimated Source Angle: ' num2str(sourceAngle) '°']);
    %legend('Directivity Pattern', 'Estimated Source');
    set(gca, 'ThetaZeroLocation', 'top', 'ThetaDir', 'clockwise');
    hold off;




    angle_traj(icount)=sourceAngle;

    current_time=toc(t0)

    current_angle=angle_traj(icount);

    subplot(2, 2, 4);
    time_vector(icount)=current_time;
    angle_vector(icount)=current_angle;
    plot(time_vector, angle_vector,'b--','Linewidth',1);
    hold on
    plot(time_vector, smooth(angle_vector,5),'r','Linewidth',2);
    legend('Real Time','Smoothed')
    xlabel('Time (s)');
    ylabel('Estiamted angle');
    title('Trajactory');
    xlim([0 maxTime])
    ylim([-90 90])
    grid on


    drawnow;
    icount=icount+1;
    %time2frame_ratio(icount)=icount/toc(t0)
end



figure
plot((1:length(angle_traj))/6,angle_traj,'b--','Linewidth',1)
hold on
plot((1:length(angle_traj))/6,smooth(angle_traj,10),'r','Linewidth',2)
ylim([-90,90])
grid on
ylabel('Estimated source angle')
xlabel('Time frame (s)')
legend('Direct output','Averaged')