close all;
clear all;

% 设置输入和输出文件夹
inputFolder = 'D:\OneDrive\OneDrive - The Hong Kong Polytechnic University\PolyU Folder\FYP\code\basic-acoustic-plot\dataset\data_nov_4'; % 替换为你的输入文件夹路径
outputFolder = 'D:\OneDrive\OneDrive - The Hong Kong Polytechnic University\PolyU Folder\FYP\code\basic-acoustic-plot\dataset\data_nov_4_modified'; % 替换为你的输出文件夹路径

% 创建输出文件夹（如果不存在）
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

% 获取输入文件夹中所有的.wav文件
audioFiles = dir(fullfile(inputFolder, '*.wav'));

% 遍历每个音频文件
for k = 1:length(audioFiles)
    % 读取音频文件
    [y, Fs] = audioread(fullfile(inputFolder, audioFiles(k).name));
    
    % 找到每个通道的最大绝对值
    maxVal = max(y, [], 'all'); % 获取所有通道的最大绝对值
    minVal = min(y, [], 'all');
    
    y_normalized = 2 * ((y - minVal) / (maxVal - minVal)) - 1;

    % 确保数据类型为整数16位（int16）
    % y_normalized = int16(y_normalized);
    
    % 构造输出文件名
    outputFileName = fullfile(outputFolder, audioFiles(k).name);
    
    % 写入处理后的音频数据到新文件
    audiowrite(outputFileName, y_normalized, Fs);

    disp(['处理完成: ', audioFiles(k).name]);
end

disp('所有音频文件处理完成！');