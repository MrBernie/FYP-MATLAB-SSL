%% 数据采集参数设置
Fs = 51200;                 % 采样率
samplesPerFrame = 2048;     % 每帧采样数
numChannels = 2;            % 双通道
bufferLength = 5;           % 时域显示缓冲时长(秒)

%% 创建NI数据采集会话
s = daq('ni');
s.Rate = Fs;
% s.IsContinuous = true;      % 连续采集模式

% 添加输入通道（根据实际硬件配置修改）
channel_1 = addinput(s, 'cDAQ1Mod1', 'ai0', 'IEPE');
% channel_1.IEPEEnable = true;    % 启用IEPE供电
% channel_1.TerminalConfig = 'SingleEnded';  % 根据硬件连接方式设置
% channel_1.Range = [-10, 10];    % 根据传感器规格调整

channel_2 = addinput(s, 'cDAQ1Mod1', 'ai1', 'IEPE');
% channel_2.IEPEEnable = true;
% channel_2.TerminalConfig = 'SingleEnded';
% channel_2.Range = [-10, 10];

%% 初始化显示参数
% 时域显示设置
timeAxis = (-bufferLength:1/Fs:0)';    % 时间轴
timeBuffer = zeros(length(timeAxis), numChannels);

% 频域显示设置
NFFT = 4096;                          % FFT点数
freqAxis = Fs/2 * linspace(0,1,NFFT/2+1); % 频率轴
win = hann(samplesPerFrame);          % 窗函数

%% 创建图形窗口（保持不变）
figure('Name','实时音频分析','NumberTitle','off');
subplot(2,1,1);
hTimePlot(1) = plot(timeAxis, timeBuffer(:,1), 'b');
hold on;
hTimePlot(2) = plot(timeAxis, timeBuffer(:,2), 'r');
xlim([-bufferLength 0]);
ylim([-0.05 0.05]); % 可能需要根据实际信号范围调整
xlabel('时间 (秒)');
ylabel('幅值');
title('实时时域波形');
legend('通道1','通道2');

subplot(2,1,2);
hFreqPlot(1) = semilogx(freqAxis, zeros(1,length(freqAxis)), 'b');
hold on;
hFreqPlot(2) = semilogx(freqAxis, zeros(1,length(freqAxis)), 'r');
xlim([20 20000]);
ylim([-120 0]);
xlabel('频率 (Hz)');
ylabel('幅值 (dB)');
title('实时频谱分析');
legend('通道1','通道2');
grid on;
set(gca,'XTick',[20 50 100 200 500 1000 2000 5000 10000 20000]);

%% 实时处理循环
start(s); % 开始连续采集
disp('开始实时采集 (按Ctrl+C停止)...');
while ishandle(gcf)
    try
        % 读取最新数据块
        data = read(s, samplesPerFrame, 'OutputFormat', 'Matrix');
        
        % 更新时域缓冲区
        timeBuffer = [timeBuffer(samplesPerFrame+1:end,:); data];
        
        % 更新时域波形
        set(hTimePlot(1), 'YData', timeBuffer(:,1));
        set(hTimePlot(2), 'YData', timeBuffer(:,2));
        
        % 计算频谱
        for ch = 1:numChannels
            frameWin = data(:,ch) .* win;
            Y = fft(frameWin, NFFT);
            P2 = abs(Y/NFFT);
            P1 = 20*log10(P2(1:NFFT/2+1));
            set(hFreqPlot(ch), 'YData', P1);
        end
        
        drawnow limitrate;
    catch ME
        disp(['错误: ' ME.message]);
        break;
    end
end

%% 清理资源
stop(s);
release(s);
disp('采集已停止');