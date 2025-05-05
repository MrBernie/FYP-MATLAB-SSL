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

% 带通滤波器参数设置
fc_low = 100;        % 低截止频率(Hz)
fc_high = 8000;      % 高截止频率(Hz)
filter_order = 12;     % 滤波器阶数（双极点）
file_ext = '*.wav';   % 支持的文件扩展名

% 获取文件列表
file_list = dir(fullfile(input_folder, file_ext));

% 批量处理循环
for n = 1:length(file_list)
    try
        % 文件命名处理
        [~, filename, ext] = fileparts(file_list(n).name);
        output_name = [filename '-bandpass-filtered' ext];
        disp(filename);
        
        % 读取音频文件
        [audio, fs] = audioread(fullfile(input_folder, file_list(n).name));
        
        % 设计带通滤波器
        [b, a] = butter(filter_order/2, [fc_low fc_high]/(fs/2), 'bandpass');
        
        % 多通道处理
        num_channels = size(audio, 2);
        processed_audio = zeros(size(audio));
        
        for i = 1:num_channels
            % 零相位数字滤波
            filtered_signal = filtfilt(b, a, audio(:, i));
            
            % 幅值补偿（保持原始响度）
            original_rms = rms(audio(:, i));
            filtered_rms = rms(filtered_signal);
            gain = original_rms / (filtered_rms + eps); % 防止除以零
            processed_audio(:, i) = filtered_signal * gain;
        end
        
        % 智能峰值限制
        peak_limit = 0.99; % 最大允许峰值
        current_peak = max(abs(processed_audio(:)));
        if current_peak > peak_limit
            processed_audio = processed_audio * (peak_limit / current_peak);
        end
        
        % 写入处理结果
        audiowrite(fullfile(output_folder, output_name), processed_audio, fs);
        
    catch ME
        warning('文件 %s 处理失败: %s', file_list(n).name, ME.message);
    end
end

winopen(output_folder);
disp('带通滤波处理完成！');
