% 读取音频文件
[fileName, filePath] = uigetfile('*.wav', '选择音频文件');
audioFile = fullfile(filePath, fileName);
[multiChannelSignal, fs] = audioread(audioFile); % 读取多通道音频文件

% 初始化参数
numChannels = size(multiChannelSignal, 2); % 获取通道数
signalPowers = zeros(numChannels, 1); % 初始化信号功率数组
noisePowers = zeros(numChannels, 1); % 初始化噪声功率数组
SNR_dB = zeros(numChannels, 1); % 初始化信噪比数组

% 带通滤波器设计（400-600 Hz）
lowCutoff = 100; % 低截止频率
highCutoff = 1500; % 高截止频率
[b, a] = butter(4, [lowCutoff, highCutoff] / (fs / 2), 'bandpass'); % 4阶巴特沃斯带通滤波器

% 计算每个通道的信号和噪声功率
for i = 1:numChannels
    currentSignal = multiChannelSignal(:, i);
    
    % 提取 400-600 Hz 的信号
    filteredSignal = filter(b, a, currentSignal);
    
    % 计算信号功率
    signalPower = var(filteredSignal);
    signalPowers(i) = signalPower;
    
    % 噪声估计：从原始信号中减去提取的信号
    noiseEstimate = currentSignal - filteredSignal;
    
    % 计算噪声功率
    noisePower = var(noiseEstimate);
    noisePowers(i) = noisePower;
    
    % 计算信噪比（分贝单位）
    if noisePower > 0
        SNR_dB(i) = 10 * log10(signalPower / noisePower);
    else
        SNR_dB(i) = Inf; % 如果噪声功率为零，设定为无穷大
    end
end

% 显示结果
disp('各通道的信号功率:');
disp(signalPowers);
disp('各通道的噪声功率:');
disp(noisePowers);
disp('各通道的信噪比 (dB):');
disp(SNR_dB);