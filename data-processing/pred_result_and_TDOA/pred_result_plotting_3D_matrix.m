%% Plot 3D diagram for a single audio prediction result

% close all;
clear all;

% 读取txt文件
[fileName, filePath] = uigetfile('*.txt', '选择预测结果文件');

% 检查用户是否选择了文件
if isequal(fileName, 0)
    disp('取消了文件选择');
    return; % 退出程序
else
    disp(['选择了文件: ', fullfile(filePath, fileName)]);
end

try
    % 尝试读取文件
    data = readmatrix(fullfile(filePath, fileName));

    % 获取数据维度
    [num_frames, num_angles] = size(data);

    % 创建坐标网格
    time_frames = 1:num_frames;
    angles = 1:num_angles;
    [X, Y] = meshgrid(time_frames * 0.2, angles);

    % 绘制3D表面图
    figure;
    surf(X', Y', data);

    % 添加坐标轴标签
    xlabel('Time Sections');
    ylabel('Angle (degree)');
    zlabel('Prediction Result');

    % 添加标题
    title_string = strjoin({fileName, 'Time-Angle-Prediction 3D'}, ' - ');
    title(title_string);

    % 可以根据需要调整视角
    % view(3); % 默认3D视角
    view([0 90]); % 自定义视角（方位角，仰角）

    % 可以添加颜色条
    colorbar;
    
    % % 绘制2D图像
    % figure;
    % imagesc(data'); % 使用imagesc函数绘制2D图像
    % 
    % % 添加坐标轴标签
    % xlabel('时间帧');
    % ylabel('角度');
    % 
    % % 添加标题
    % title_string = strjoin({fileName, '时间帧-角度 2D 图'}, ' - ');
    % title(title_string);
    % 
    % % 设置 y 轴的刻度，显示角度值
    % % yticks(1:num_angles); % 确保有足够的刻度
    % % yticklabels(string(angles)); % 假设你的角度是 1 到 180
    % 
    % % 添加颜色条
    % colorbar;
    % 
    % % 调整布局以防止重叠
    % sgtitle(['文件: ', fileName]); % 总标题，显示文件名
    drawnow; % 确保图形更新
    
catch ME
    % 捕获错误
    disp(['读取文件时发生错误: ', ME.message]);
end

%% Plot 3D Diagram for entire folder and saved as files

% 清除工作区
clear all;
close all;

% 可视化路径选择界面
try
    input_folder = uigetdir('', '请选择输入文件夹');
    if input_folder == 0
        error('用户取消输入文件夹选择');
    end
    
    output_folder = uigetdir('', '请选择输出文件夹');
    if output_folder == 0
        error('用户取消输出文件夹选择');
    end
catch ME
    errordlg(['路径选择错误: ' ME.message], '系统错误');
    return;
end

% 自动创建输出目录
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

% 每个音频时间长度（手动修改）(秒)
total_time_length = 15;

%input_folder底下子文件夹
pred_matrix_dir = fullfile(input_folder, 'pred_matrix');

% 获取所有txt文件列表
file_list = dir(fullfile(pred_matrix_dir, '*.txt'));

% 遍历处理每个文件
for n = 1:length(file_list)
    current_file = fullfile(pred_matrix_dir, file_list(n).name);
    
    try
        % 读取数据矩阵
        data = readmatrix(current_file);
        [num_frames, num_angles] = size(data);
        seconds_per_frames = total_time_length/num_frames;
        
        % 创建新图形窗口
        fig = figure('Visible', 'off'); % 不显示图形窗口
        set(fig, 'Position', [0 0 1280 720]); % 设置分辨率画布
        
        % 绘制3D表面图
        X = linspace(0,total_time_length,num_frames);
        Y = 1:num_angles;
        surf(X, Y, data');
        
        % 图形标注设置
        xlabel('Time Frames');
        ylabel('Angle (degree)');
        zlabel('Prediction Value');
        title_str = strrep(file_list(n).name, '_', '\_'); % 转义下划线
        title([title_str ' - 3D Plot']);
        % view(3); % 标准3D视角
        view([0 90]);
        colorbar;
        grid on;
        
        % 设置保存路径
        [~, name] = fileparts(file_list(n).name);
        output_path = fullfile(output_folder, [name '_3Dplot.png']);
        
        % 保存为高清图片（可修改为.fig保存原始图形）
        exportgraphics(fig, output_path, 'Resolution', 300);
        close(fig); % 关闭图形释放内存
        
        fprintf('成功处理: %s\n', file_list(n).name);
        
    catch ME
        fprintf('文件 %s 处理失败: %s\n', file_list(n).name, ME.message);
        if exist('fig', 'var') && ishandle(fig)
            close(fig); % 确保异常时关闭图形
        end
    end
end

disp('批量处理完成！');

%% 语音活动检测可视化增强版（带动态时间窗和亮度调节）

clear all;
close all;
clc;

try
    input_folder = uigetdir('', '请选择输入文件夹');
    if input_folder == 0
        error('用户取消输入文件夹选择');
    end
    
    output_folder = uigetdir('', '请选择输出文件夹');
    if output_folder == 0
        error('用户取消输出文件夹选择');
    end
catch ME
    errordlg(['路径选择错误: ' ME.message], '系统错误');
    return;
end

if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

% 参数配置
fs = 16000;           % 采样率
min_saturation = 0.4; % 最小饱和度
min_brightness = 0.4; % 最小亮度
total_time_length = 15; % 单个音频时长

pred_matrix_dir = fullfile(input_folder, 'pred_matrix');
vad_out_dir = fullfile(input_folder, 'vad_out');
original_audio_dir = fullfile(input_folder, 'original_signal');

file_list = dir(fullfile(pred_matrix_dir, '*.txt'));

for n = 1:length(file_list)
    current_file = fullfile(pred_matrix_dir, file_list(n).name);
    
    try
        % ========== 数据读取 ==========
        data = readmatrix(current_file);
        [num_frames, num_angles] = size(data);
        
        % ========== VAD文件处理 ==========
        [~, pred_name] = fileparts(file_list(n).name);
        vad_name = strrep(pred_name, '-pred-matrix', '-vad-out');
        vad_file = fullfile(vad_out_dir, [vad_name '.txt']);
        
        if ~exist(vad_file, 'file')
            error('VAD文件不存在: %s', vad_file);
        end
        vad_data = readmatrix(vad_file);
        
        % ========== 动态时间窗计算 ==========
        total_samples = size(vad_data, 1);
        samples_per_frame = floor(total_samples / num_frames);
        remainder_samples = total_samples - num_frames*samples_per_frame;
        
        % ========== 计算VAD均值 ==========
        avg_values = zeros(1, num_frames);
        for i = 1:num_frames
            % 处理余数采样点
            if i == num_frames && remainder_samples > 0
                start_idx = (i-1)*samples_per_frame + 1;
                end_idx = start_idx + samples_per_frame + remainder_samples - 1;
            else
                start_idx = (i-1)*samples_per_frame + 1;
                end_idx = i*samples_per_frame;
            end
            
            % 提取VAD数据并计算均值
            window_vad = vad_data(start_idx:end_idx, :);
            avg_values(i) = mean(window_vad(:), 'double');
        end
        
        % ========== 颜色映射处理 ==========
        % 转置数据以匹配维度
        cdata = data';  % 现在维度为 [num_angles × num_frames]
        
        % 亮度调整矩阵
        brightness_scale = avg_values * (1 - min_brightness) + min_brightness;
        brightness_matrix = repmat(brightness_scale, num_angles, 1);
        
        % 按列归一化预测值
        cdata = normalize(cdata,'range');
        % 画图归一化预测值
        cmin = min(cdata(:));
        cmax = max(cdata(:));
        normalized_data = (cdata - cmin) / (cmax - cmin);
        
        % 生成基础颜色
        cmap = parula(256);
        rgb_data = ind2rgb(round(normalized_data*255), cmap);
        
        % 调整亮度（HSV空间操作）
        hsv_data = rgb2hsv(rgb_data);
        % hsv_data(:,:,2) = hsv_data(:,:,2) .* brightness_matrix;
        hsv_data(:,:,3) = hsv_data(:,:,3) .* brightness_matrix;
        adjusted_rgb = hsv2rgb(hsv_data);
        
        % ========== 图形绘制 ==========
        fig = figure('Visible', 'off');
        set(fig, 'Position', [0 0 1280 720], 'Color', 'white');
        
        % 绘制表面图
        X = linspace(0,total_time_length,num_frames);
        Y = 1:num_angles;
        % ax1 = axes;
        surf(X, Y, cdata,...
             'FaceColor', 'texturemap',...
             'CData', adjusted_rgb);
        
        % 视图设置
        view([0 90]);
        axis tight;
        xlabel(sprintf('Time (Second)'));
        ylabel('Angle (Degree)');
        title_str = strrep(file_list(n).name, '_', '\_');
        title([title_str ' - Prediction Result with VAD']);
        colorbar;
        caxis([cmin cmax]);
        grid on;
        
        % ========== 保存输出 ==========
        output_path = fullfile(output_folder, [pred_name '_enhanced.png']);
        exportgraphics(fig, output_path, 'Resolution', 300);
        close(fig);
        
        fprintf('成功处理: %s\n', file_list(n).name);

        % ====== 原始声音信号处理 =======
        % [~, pred_name] = fileparts(file_list(n).name);
        % vad_name = strrep(pred_name, '-pred-matrix', '-vad-out');
        % vad_file = fullfile(vad_out_dir, [vad_name '.txt']);
        % 
        % if ~exist(vad_file, 'file')
        %     error('VAD文件不存在: %s', vad_file);
        % end
        % vad_data = readmatrix(vad_file);

        
    catch ME
        fprintf('文件 %s 处理失败: %s\n', file_list(n).name, ME.message);
        if exist('fig', 'var') && ishandle(fig)
            close(fig);
        end
    end
end

disp('批量处理完成！');


