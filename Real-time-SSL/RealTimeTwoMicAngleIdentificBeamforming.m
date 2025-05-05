%% 2025-03-29
%% Two microphone beamforming program
%% For angle identification
%% Real Time-Audio Interface with two cheap microphones

clc
clear all
close all

% Parameters
fs = 44100; % Sample rate (Hz)
maxTime = 60; % Maximum recording time
%duration = 1; % Recording duration (s)
micDist = 0.16; % Distance between microphones (m)
numAngles = 180; % Number of angles to evaluate
speedOfSound = 343; % Speed of sound (m/s)

N = 1024; % FFT length

win = hann(N); % Hanning window
%win = ones(N,1); % Rectangular window

% Set up the audio device reader
reader = audioDeviceReader('SampleRate', fs, 'NumChannels', 2);

% Set up the frequency vector
f = (0:N/2-1)/N*fs;

t0 = tic;

icount = 1;
while toc(t0) < maxTime

    % Read a frame from the audio device
    x = reader();

    % Apply the window function
    x = x .* win;

    % Get audio data
    channel1 = x(:, 1);
    channel2 = x(:, 2);

    % 获取频域数据
    channel1_fft = fft(channel1, N);
    channel2_fft = fft(channel2, N);

    % 构建完整频率向量
    f_full = (0:N-1)*(fs/N);
    f_full(f_full >= fs/2) = f_full(f_full >= fs/2) - fs;  % 转为正负频率

    % 波束形成计算
    angles = linspace(-90, 90, numAngles);
    directivity = zeros(1, numAngles);
    
    for i = 1:numAngles
        theta = angles(i);
        tau = (micDist/speedOfSound)*sind(theta);
        
        % 生成相位补偿因子
        phase_shift = exp(-1i*2*pi*f_full*tau).';
        
        % 波束形成计算
        beam_signal = channel1_fft + channel2_fft .* phase_shift;
        
        % 计算能量（考虑全频段）
        directivity(i) = sum(abs(beam_signal).^2);
    end

    % 找到最大能量对应的角度
    [~, angleIdx] = max(directivity);
    sourceAngle = angles(angleIdx);

    % Plot frequency spectrum of each channel
    frequencies = (0:N/2)*(fs/N);
    channel1_fft_data = channel1_fft(1:N/2+1);
    channel2_fft_data = channel2_fft(1:N/2+1);
    psdx1 = (1/(fs*N)) * abs(channel1_fft_data).^2;
    psdx1(2:end-1) = 2*psdx1(2:end-1);
    psdx2 = (1/(fs*N)) * abs(channel2_fft_data).^2;
    psdx2(2:end-1) = 2*psdx2(2:end-1);
    Chanel_Sig1 = psdx1;
    Chanel_Sig1_s = smooth(Chanel_Sig1, 50);
    Chanel_Sig2 = psdx2;
    Chanel_Sig2_s = smooth(Chanel_Sig2, 50);

    subplot(2, 2, 1); % Plot Mic 1 FFT data and smoothed
    plot(frequencies, 10*log10(Chanel_Sig1));
    legend('Frequency Spectrum')
    xlabel('Frequency (Hz)');
    ylabel('Channel 1');
    title(['Channel 1 Frequency Domain']);
    axis([100 10000 -160 -80])
    grid on

    subplot(2, 2, 3); % Plot Mic 1 FFT data and smoothed
    plot(frequencies, 10*log10(Chanel_Sig2));
    legend('Frequency Spectrum')
    xlabel('Frequency (Hz)');
    ylabel('Channel 2');
    title(['Channel 2 Frequency Domain']);
    axis([100 10000 -160 -80])
    grid on

    % Plot directivity pattern
    subplot(2, 2, 2);
    polarplot(deg2rad(angles), directivity/max(directivity), 'LineWidth', 2);
    hold on;
    polarplot(deg2rad(sourceAngle), directivity(angleIdx)/max(directivity), 'ro', 'MarkerSize', 10, 'LineWidth', 2);
    title(['Estimated Source Angle: ' num2str(sourceAngle) '°']);
    set(gca, 'ThetaZeroLocation', 'top', 'ThetaDir', 'clockwise');
    hold off;

    % Record angle trajectory
    angle_traj(icount) = sourceAngle;

    current_time = toc(t0);
    current_angle = angle_traj(icount);

    subplot(2, 2, 4);
    time_vector(icount) = current_time;
    angle_vector(icount) = current_angle;
    plot(time_vector, angle_vector, 'b--', 'LineWidth', 1);
    hold on;
    plot(time_vector, smooth(angle_vector, 5), 'r', 'LineWidth', 2);
    legend('Real Time', 'Smoothed')
    xlabel('Time (s)');
    ylabel('Estimated Angle');
    title('Trajectory');
    xlim([0 maxTime])
    ylim([-90 90])
    grid on

    drawnow;
    icount = icount + 1;
end

figure
plot((1:length(angle_traj))/6, angle_traj, 'b--', 'LineWidth', 1)
hold on
plot((1:length(angle_traj))/6, smooth(angle_traj, 10), 'r', 'LineWidth', 2)
ylim([-90, 90])
grid on
ylabel('Estimated Source Angle')
xlabel('Time Frame (s)')
legend('Direct Output', 'Averaged')
