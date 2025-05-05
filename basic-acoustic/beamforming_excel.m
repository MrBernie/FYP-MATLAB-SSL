clear all;
close all;

% --------------------------------Read Files-------------------------------%
% 选择音频文件夹
folderPath = uigetdir('', '选择音频文件夹');
audioFiles = dir(fullfile(folderPath, '*.wav')); % 获取所有.wav文件

% 创建一个Excel文件用于保存结果
outputFile = fullfile(folderPath, 'results.xlsx');
headers = {'FileName', 'low_cutoff','high_cutoff','frequency', 'LocationX', 'LocationY', 'DegreeOfArrival'};
xlswrite(outputFile, headers); % 写入表头

% 初始化行索引
rowIndex = 2; % 从第二行开始写入数据

for fileIdx = 1:length(audioFiles)
    % 读取音频文件
    audioFile = fullfile(folderPath, audioFiles(fileIdx).name);
    disp(['Processing data ',audioFiles(fileIdx).name])
    [p, fs] = audioread(audioFile); % 读取多通道音频文件
    p = double(p);
    dt = 1/fs;
    t = (0:length(p)-1)/fs;
    t = transpose(t);

    % ------------------------------Frequency filter----------------------------------%
    low_cutoff = 1900; % 下限频率
    high_cutoff = 2100; % 上限频率
    f = 2000; % 关心的频率

    % 使用巴特沃斯滤波器设计函数
    [b, a] = butter(4, [low_cutoff high_cutoff] / (fs / 2), 'bandpass');

    % 应用巴特沃斯滤波器
    p = filter(b, a, p);

    % % ------------------------------Normalize original signal-------------------------------%
    % pn_max = zeros(size(p, 1), 1);
    % for i = 1:size(p,1)
    %     pn_max(i) = max(p(i,:)); % y的第i行的最大元素的值
    %     pn(i,:) = p(i,:)/pn_max(i);
    % end

    %% ------------------------------Beamforming-------------------------------%
    T = length(p)/fs;
    R = p*p'/T; % 接收数据的自协方差矩阵  
    w = 2*pi*f;  % 角频率
    c = 334; % 声速

    zi = zeros(6,1); % 生成一个M*1维的零矩阵
    yi = [0; 0.052; 0.052; 0; -0.052; -0.052]*1;
    xi = [-0.06; -0.03; 0.03; 0.06; 0.03; -0.03]*1;

    step_x = 0.05;  
    step_y = 0.05;
    d_x = 1;
    d_y = 1;
    
    y = (-1*d_x:step_y:d_x);
    x = (-1*d_x:step_x:d_x);  
    z = 0;

    for k1=1:length(y)
        for k2=1:length(x)
            Ri = sqrt((x(k2)-xi).^2+(y(k1)-yi).^2+(z-zi).^2);  
            Ri2 = sqrt((x(k2)-0).^2+(y(k1)-0).^2+(z-0).^2);
            Rn = Ri-Ri2;   
            b = exp(-j*w*Rn/c); 
            Pcbf(k1,k2) = abs(b'*R*b); 
        end
    end

    %% -------------------------------------Normalization-------------------------------------%
    for k1 = 1:length(y)
        pp(k1) = max(Pcbf(k1,:)); 
    end

    Pcbf = Pcbf/max(pp);
    
    [maxValue, linearIndex] = max(Pcbf, [], 'all');
    [row, column] = ind2sub(size(Pcbf), linearIndex);
    
    y_max = row*step_x-d_x;
    x_max = column*step_y-d_y;

    disp(['Processing file: ', audioFiles(fileIdx).name]);
    
    % 计算角度（弧度）
    theta_rad = atan2(y_max, x_max);
    
    % 转换为度
    theta_deg = rad2deg(theta_rad);

    disp(['Maximum number ', num2str(maxValue), ', location:(', num2str(x_max), ',', num2str(y_max), ')']);
    disp(['Index of the maximum ', num2str(row), ',', num2str(column)]);
    disp(['Degree of Arrival ', num2str(theta_deg)]);

    % 将结果写入Excel文件
    resultsRow = {audioFiles(fileIdx).name, low_cutoff,high_cutoff,f, x_max, y_max, theta_deg};
    
    xlswrite(outputFile, resultsRow, 'Sheet1', ['A' num2str(rowIndex)]); % 写入数据
    rowIndex = rowIndex + 1; % 更新行索引
end

disp('所有音频文件处理完成，结果已写入Excel。');