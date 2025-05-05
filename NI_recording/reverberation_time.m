function RT60_values = calculate_RT60_per_band(signal, fs, bands_low, bands_high)
    % 计算输入信号在不同频段的RT60
    % signal: 输入的声音信号
    % fs: 采样频率
    % bands_low: 频段下频率（Hz）的数组
    % bands_high: 频段上频率（Hz）的数组
    
    % 确保信号为列向量
    signal = signal(:);
    
    % 初始化RT60值数组
    RT60_values = zeros(length(bands_low), 1);
    
    for i = 1:length(bands_low)
        % 使用带通滤波器提取特定频段
        band_signal = bandpass(signal, [bands_low(i), bands_high(i)], fs);
        
        % 计算该频段的RT60
        RT60_values(i) = calculate_RT60(band_signal, fs, bands_low(i), bands_high(i));
    end
end

function RT60 = calculate_RT60(signal, fs, band_low, band_high)
    % 计算输入信号的RT60
    signal = signal(:);
    
    start_decay = -5; % 指定开始衰减的信噪比(dB)
    stop_decay = -20; % 指定结束衰减的信噪比:(dB)
    window_size = round(0.01*fs); % 窗口
    hop_size = round(0.01*fs); % 步长

    function RMS = rms(signal)
        RMS = sqrt(mean(signal.^2));
    end
    signal_pressure_begin = rms(signal(1:window_size));
    reference_pressure = 20e-6; % 20 µPa
    spl_begin = 20 * log10(signal_pressure_begin/reference_pressure);

    window_count = 0;

    for i = 1:hop_size:length(signal)-window_size
        windowed_signal = signal(i:i+window_size-1);
        window_signal_pressure = rms(windowed_signal);
        spl = 20*log10(window_signal_pressure/reference_pressure);

        if spl > spl_begin
            spl_begin = spl;
            window_count = 0;
        end

        if spl - spl_begin < start_decay
            window_count = window_count+1;
        end

        if spl - spl_begin < stop_decay
            break;
        end
    end
    decay_time_ratio = -60/(stop_decay-start_decay); % 计算哀减分贝比
    RT60 = window_count * window_size / fs * decay_time_ratio;

    spl_values = [];
    for i = 1:hop_size:length(signal)-window_size
        window = signal(i:i+window_size-1);
        window_signal_power = rms(window);
        spl = 20*log10(window_signal_power/reference_pressure);
        spl_values(end+1) = spl;
    end

    % 绘制SPL随时间变化的图形
    time_vector = (0:length(spl_values)-1) * (hop_size / fs); % 时间向量（秒）

    figure;
    plot(time_vector, spl_values, 'LineWidth', 0.8);
    title_string = sprintf('SPL over Time (%.0fHz ~ %.0fHz)', band_low, band_high);
    title(title_string);
    xlabel('Time (seconds)');
    ylabel('SPL (dB)');
    grid on;
end

close all;
clear all;

% 读取音频文件
[fileName, filePath] = uigetfile('*.wav', '选择音频文件');
audioFile = fullfile(filePath, fileName);
[signal, fs] = audioread(audioFile); % 读取多通道音频文件
signal = signal(:,1);
signal = signal((2.9*fs):end-1); % 截取需要是时间

% 定义要分析的频段（Hz）
% bands_low = [125, 250, 500, 1000, 2000, 4000]; 
% bands_high = [250, 500, 1000, 2000, 4000, 8000];
% bands_low = 1000:1000:7000;
% bands_high = 2000:1000:8000;
bands_mid = [125, 250, 500, 1000, 2000, 4000, 8000, 16000];
bands_low = bands_mid / (2^(1/6));
bands_high = bands_mid * (2^(1/6));

% 确保信号为列向量
signal = signal(:);

% 初始化RT60值数组
RT60_values = zeros(length(bands_low), 1);

for i = 1:length(bands_low)
    % 使用带通滤波器提取特定频段
    band_signal = bandpass(signal, [bands_low(i), bands_high(i)], fs);

    % 计算该频段的RT60
    RT60_values(i) = calculate_RT60(band_signal, fs, bands_low(i), bands_high(i));
end

% 输出每个频段的RT60值
for i = 1:length(bands_low)
    fprintf('中心频率为 %.0f Hz 的RT60: %.4f \n', bands_mid(i), RT60_values(i));
end



% % 使用带通滤波器提取特定频段
% band_signal = bandpass(signal, [bands_low(4), bands_high(4)], fs);
% 
% % 计算该频段的RT60
% RT60_values = calculate_RT60(band_signal, fs, bands_low(4), bands_high(4));
% fprintf('频率 %.0f Hz 到频率 %.0f HZ的 RT60: %.3f 秒\n', bands_low(4), bands_high(4), RT60_values);
