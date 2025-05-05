clear all;
close all;

%% --------------------------------------------------------------
% 设置要遍历的文件夹路径
folderPath = 'D:\OneDrive\OneDrive - The Hong Kong Polytechnic University\PolyU Folder\FYP\code\basic-acoustic-plot\dataset\nov_11_pred_result'; % 替换为你的文件夹路径

% 获取该文件夹下所有的 .txt 文件
files = dir(fullfile(folderPath, '*.txt'));

% 初始化一个 cell 数组来存储文件名和数据
data = cell(length(files), 2);

% 遍历每个文件并读取数据
for k = 1:length(files)
    % 获取当前文件的完整路径
    fileName = fullfile(folderPath, files(k).name);
    
    % 读取当前文件的数据
    % 假设每个文件中只有一列数字
    fileData = readmatrix(fileName); % 使用 readmatrix 读取数据
    
    % 将文件名存储在第一列
    data{k, 1} = files(k).name;
    
    % 将读取的数据存储在第二列
    data{k, 2} = fileData; % 存储整个数据列
end

% 显示结果
disp('所有文件的数据已读取。');
disp(data);

%% --------------------------------------------------------------
for i=1:size(data,1)

    % 获取第一个 cell 的第二列数据
    yData = data{i, 2}; % 读取第一个文件的数据
    
    % 确保 yData 是一个列向量
    if iscolumn(yData)
        yData = yData; % 保持原样
    else
        yData = yData'; % 转置为列向量
    end
    
    % 获取 n 的大小
    n = length(yData);
    
    % 创建 x 数据（横坐标）
    xData = 1:n; % 从 1 到 n 的数组
    
    % 绘制折线图
    figure; % 创建新图形窗口
    plot(xData, yData, '-o'); % '-o' 表示带有圆点的线条
    xlabel('Predict Point (n)'); % 横坐标标签
    ylabel('Angle (degrees)'); % 纵坐标标签
    fileName = data{i, 1}; % 获取第一个文件的文件名
    titleStr = sprintf('Result from the coordinate: %s',fileName(1:end-4));
    title(titleStr);
    grid on; % 显示网格
end



