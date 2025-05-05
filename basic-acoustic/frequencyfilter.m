close all;
clear all;

% 读取音频文件
[fileName, filePath] = uigetfile('*.wav', '选择音频文件');
audioFile = fullfile(filePath, fileName);
[p, fs] = audioread(audioFile); % 读取多通道音频文件
p = double(p);
dt = 1/fs;
t = (0:length(p)-1)/fs;
t = transpose(t);

% 设计带通滤波器
low_cutoff = 950; % 下限频率
high_cutoff = 1050; % 上限频率

% 使用巴特沃斯滤波器设计函数
[b, a] = butter(4, [low_cutoff high_cutoff] / (fs / 2), 'bandpass');

% 应用滤波器
pf = filter(b, a, p);

%% ------------------------------Plot original Signal-------------------------------%
figure;
tiledlayout(6,1)
for i = 1 : 6
    nexttile();
    plot(t, p(:,i));  
    title(['Original Signal Channel ' num2str(i)]);
end

%% ------------------------------Plot filtered Signal-------------------------------%
figure;
tiledlayout(6,1)
for i = 1 : 6
    nexttile();
    plot(t, pf(:,i));  
    title(['Filtered Signal Channel ' num2str(i)]);
end

%% ------------------------------Write filtered signal to audio file-------------------------------%
output_path = 'filtered_audio.wav'; % 指定输出路径和文件名
audiowrite(output_path, pf, fs); % 写入音频文件

disp(['Filtered audio saved to: ' output_path]);