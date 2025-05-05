close all;
clear all;
%% --
% 读取音频文件
[fileName, filePath] = uigetfile('*.wav', '选择音频文件');
audioFile = fullfile(filePath, fileName);
[signal, fs] = audioread(audioFile); % 读取双通道音频文件

ch_1 = signal(:,1);
ch_2 = signal(:,2);

% time_windows = 1/fs; % 时间窗对应的时间
% p_ref = 20e-6;
% spl_ch1 = 20 * log10(spl_ch1) / p_ref;
% spl_ch2 = 20 * log10(spl_ch2) / p_ref;

%% -- SPL

% 参数设置
window_duration = 0.005; % 时间窗长度 (秒)
window_samples = round(window_duration * fs); % 每个时间窗的采样点数
num_windows = floor(length(ch_1) / window_samples); % 时间窗数量

% 初始化声压级数组
spl_ch1 = zeros(num_windows, 1);
spl_ch2 = zeros(num_windows, 1);
time_windows = (0:num_windows-1) * window_duration; % 时间窗对应的时间

% 计算每个时间窗的声压级
for i = 1:num_windows
    % 获取当前时间窗的数据
    start_idx = (i-1) * window_samples + 1;
    end_idx = i * window_samples;
    window_data_ch1 = ch_1(start_idx:end_idx);
    window_data_ch2 = ch_2(start_idx:end_idx);

    % 声压级计算公式：SPL = 20 * log10(rms(p) / p_ref)
    p_ref = 20e-6; % 声压参考值 (20 µPa)
    spl_ch1(i) = 20 * log10(rms(window_data_ch1) / p_ref);
    spl_ch2(i) = 20 * log10(rms(window_data_ch2) / p_ref);
end

%% ---- T60 Dumb
decay_start_time = 2.44;
decay_end_time = 2.45;
spl_diff = spl_ch1(round(decay_start_time / window_duration)) - spl_ch1(round(decay_end_time / window_duration));
time_diff = decay_end_time - decay_start_time;
t60_diff = 60/spl_diff;
t60 = time_diff*t60_diff;
fprintf('Decay start time: %.2f seconds\n', decay_start_time);
fprintf('Decay end time: %.2f seconds\n', decay_end_time);
fprintf('T60: %.2f seconds\n', t60);


%% ----  T60
% 
% % 找到衰减开始时间
% % decay_start_time = find_decay_start(spl_ch1, window_duration, 500);
% decay_start_time = 2.44;
% 
% if isnan(decay_start_time)
%     error('未找到衰减开始时间');
% end
% 
% % 找到初始声压级（第一个窗口）
% initial_spl = spl_ch1(round(decay_start_time / window_duration));
% 
% % 找到 T10
% decay_threshold = initial_spl - 10; % 衰减10 dB
% decay_indices = find(spl_ch1 <= decay_threshold, 1, 'first');
% decay_end_time = time_windows(decay_indices);
% 
% if ~isempty(decay_indices)
%     t10 = decay_end_time - decay_start_time; % T10 对应的时间
% else
%     t10 = NaN; % 如果没有找到衰减到10 dB，则设为NaN
% end
% 
% % 计算 T60
% t60 = t10 * 6;
% 
% % 输出结果
% fprintf('Decay start time: %.2f seconds\n', decay_start_time);
% fprintf('Decay end time: %.2f seconds\n', decay_end_time);
% fprintf('T10: %.2f seconds\n', t10);
% fprintf('T60: %.2f seconds\n', t60);

%% -- Plot
% 绘制声压级随时间变化图
figure;

% 通道1声压级图
subplot(2, 1, 1);
plot(time_windows, spl_ch1, 'b', 'LineWidth', 0.8);
hold on;
scatter(decay_start_time, spl_ch1(round(decay_start_time / window_duration)), 80, 'r', 'filled', 'DisplayName', 'Decay Start');
scatter(decay_end_time, spl_ch1(round(decay_end_time / window_duration)), 80, 'g', 'filled', 'DisplayName', 'T20 End');
title('Sound Pressure Level vs Time (Channel 1)');
xlabel('Time (s)');
ylabel('SPL (dB)');
grid on;

% 通道2声压级图
subplot(2, 1, 2);
plot(time_windows, spl_ch2, 'r', 'LineWidth', 0.8);
title('Sound Pressure Level vs Time (Channel 2)');
xlabel('Time (s)');
ylabel('SPL (dB)');
grid on;

% 总标题
sgtitle('Sound Pressure Level vs Time for Both Channels');

%% -- Functions


% 定义衰减开始的阈值（例如，信号能量下降到初始能量的 90%）
function decay_start_time = find_decay_start(spl, window_duration, slope_threshold)
    % 计算声压级的变化率（斜率）
    spl_diff = diff(spl); % 计算相邻声压级之间的差值
    time_diff = window_duration; % 每个时间窗的长度

    % 计算斜率
    slopes = spl_diff / time_diff; % 计算每个窗的斜率

    % 找到斜率超过阈值的第一个点
    for i = 1:length(slopes)
        if slopes(i) < -slope_threshold % 斜率小于负阈值，表示下降
            decay_start_time = (i + 1) * window_duration; % 返回对应的时间
            return;
        end
    end

    decay_start_time = NaN; % 如果没有找到衰减点，返回 NaN
end
