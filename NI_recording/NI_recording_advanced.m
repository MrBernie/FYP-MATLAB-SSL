%% 清理环境
clear;
clc;

%% 创建NI数据采集会话
s = daq('ni');

%% 添加模拟输入通道配置
% 根据实际硬件修改设备ID和通道号
channel_1 = addinput(s, 'cDAQ1Mod1', 'ai0', 'IEPE');
channel_2 = addinput(s, 'cDAQ1Mod1', 'ai1', 'IEPE');

%% 配置采集参数
Fs = 51200;              % 采样率
duration = 5;            % 录制时长(秒)
samplesPerFrame = 2048;  % 每帧采样数
s.Rate = Fs;             % 设置采样率
totalSamples = Fs * duration;

%% 预分配存储空间
audioData = zeros(totalSamples, 2);
currentIndex = 1;

%% 创建进度条
h = waitbar(0, '准备开始采集...', 'Name','NI采集进度');
s.IsContinuous = true;  % 启用连续采集模式

try
    start(s);  % 启动后台采集
    
    %% 分帧采集循环
    while currentIndex <= totalSamples
        % 计算剩余采样数
        remaining = totalSamples - currentIndex + 1;
        readSamples = min(samplesPerFrame, remaining);
        
        % 读取当前帧数据
        [frameData, ~] = read(s, readSamples, 'OutputFormat','Matrix');
        
        % 存储数据
        audioData(currentIndex:currentIndex+readSamples-1, :) = frameData;
        currentIndex = currentIndex + readSamples;
        
        % 更新进度条
        progress = currentIndex / totalSamples;
        waitbar(progress, h, ...
            sprintf('采集进度: %.1f%% (已采集 %.1fs/%.1fs)',...
            progress*100, currentIndex/Fs, duration));
    end
    
catch ME
    % 异常处理
    stop(s);
    release(s);
    delete(h);
    rethrow(ME);
end

%% 停止采集并清理资源
stop(s);
release(s);
delete(h);

%% 数据后处理与保存
% 归一化处理
audioData = audioData / max(abs(audioData(:)));

% 保存为WAV文件
filename = 'NI_recording.wav';
audiowrite(filename, audioData, Fs);
fprintf('文件已保存: %s\n', filename);

%% 时域波形显示
figure;
t = (0:totalSamples-1)/Fs;

subplot(2,1,1);
plot(t, audioData(:,1));
grid on;
title('NI通道1 - 时域波形');
xlabel('时间 (s)'), ylabel('幅值');
xlim([0 duration]);

subplot(2,1,2);
plot(t, audioData(:,2), 'r'), grid on;
title('NI通道2 - 时域波形');
xlabel('时间 (s)'), ylabel('幅值');
xlim([0 duration]);