%% 音频采集参数设置
Fs = 51200;                 % 采样率
samplesPerFrame = 2048;     % 每帧采样数(平衡实时性和延迟)
numChannels = 2;            % 双通道
bufferLength = 5;           % 时域显示缓冲时长(秒)

%% 创建音频采集对象
% 修改为实际设备名称
deviceReader = audioDeviceReader(...
    'Device', '麦克风 (USB Audio CODEC )', ...  
    'SampleRate', Fs,...
    'NumChannels', numChannels,...
    'SamplesPerFrame', samplesPerFrame);

%% 初始化显示参数
% 时域显示设置
timeAxis = (-bufferLength:1/Fs:0)';    % 时间轴(滚动显示)
timeBuffer = zeros(length(timeAxis), numChannels);

% 频域显示设置
NFFT = 4096;                          % FFT点数
freqAxis = Fs/2 * linspace(0,1,NFFT/2+1); % 频率轴(单边谱)
win = hann(samplesPerFrame);          % 窗函数

%% 创建图形窗口
figure('Name','实时音频分析','NumberTitle','off');
subplot(2,1,1);
hTimePlot(1) = plot(timeAxis, timeBuffer(:,1), 'b');
hold on;
hTimePlot(2) = plot(timeAxis, timeBuffer(:,2), 'r');
xlim([-bufferLength 0]);
ylim([-0.25 0.25]);
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
disp('开始实时采集 (按Ctrl+C停止)...');
while ishandle(gcf) % 窗口存在时持续运行
    % 读取音频数据
    audioFrame = deviceReader();
    
    % 更新时域缓冲区
    timeBuffer = [timeBuffer(samplesPerFrame+1:end,:); audioFrame];
    
    % 更新时域波形
    for ch = 1:numChannels
        set(hTimePlot(ch), 'YData', timeBuffer(:,ch));
    end
    
    % 计算频谱
    for ch = 1:numChannels
        frameWin = audioFrame(:,ch) .* win;
        Y = fft(frameWin, NFFT);
        P2 = abs(Y/NFFT);
        P1 = 20*log10(P2(1:NFFT/2+1)); % 转换为dB
        
        % 更新频谱图
        set(hFreqPlot(ch), 'YData', P1);
    end
    
    % 强制刷新图形
    drawnow limitrate;
end

%% 清理资源
release(deviceReader);
disp('采集已停止');