%% Cut by seconds ---------------------------------------------------------
% 选择音频文件
[filename, pathname] = uigetfile({'*.wav;*.mp3;*.flac','音频文件 (*.wav,*.mp3,*.flac)'});
if isequal(filename,0)
    disp('用户取消选择');
    return;
end

% 读取完整音频
[y, Fs] = audioread(fullfile(pathname, filename));

% 设置分段时长（秒）
segment_duration = input('请输入分段时长（秒）：');
samples_per_segment = round(segment_duration * Fs);

% 计算分段参数
total_samples = size(y, 1);
num_segments = floor(total_samples / samples_per_segment);

% 拆分文件名
[~, basename, ext] = fileparts(filename);

% 创建分段并保存
for i = 1:num_segments
    start_sample = (i-1)*samples_per_segment + 1;
    end_sample = i*samples_per_segment;
    
    segment = y(start_sample:end_sample, :);
    
    new_filename = sprintf('%s-%d%s', basename, i, ext);
    audiowrite(fullfile(pathname, new_filename), segment, Fs);
end

disp(['成功分割为 ', num2str(num_segments), ' 个片段']);

%% Cut Noise
% 选择音频文件
[filename, path] = uigetfile({'*.wav;*.mp3;*.flac','音频文件 (*.wav,*.mp3,*.flac)'});
if isequal(filename,0)
    disp('用户取消选择');
    return;
end

start_seconds = 5;
end_seconds = 6;

% 读取完整音频
[y, Fs] = audioread(fullfile(path, filename));

% 计算总时长
total_seconds = size(y,1)/Fs; 

% 获取a-b秒的样本索引
start_sample = start_seconds * Fs;
end_sample = end_seconds * Fs;

% 提取单通道（默认取左声道）
channel = 1; % 修改此值选择其他声道
segmented = y(start_sample:end_sample, channel); 

% 生成新文件名
[~, name, ext] = fileparts(filename);
new_filename = [name '-noise' ext];

% 写入文件
audiowrite(fullfile(path, new_filename), segmented, Fs);
disp(['成功保存：' new_filename]);
%% Copy and paste sound section -------------------------------------------
% 选择音频文件
[filename, path] = uigetfile({'*.wav;*.mp3;*.flac','音频文件 (*.wav,*.mp3,*.flac)'});
if isequal(filename,0)
    disp('操作已取消');
    return;
end

% 读取音频数据
[y, Fs] = audioread(fullfile(path, filename));
original_length = size(y,1);  % 原始样本数

% 输入目标时长
target_seconds = 235;
target_samples = round(target_seconds * Fs);  % 目标样本数

% 计算需要复制的次数
repeat_times = ceil(target_samples / original_length);

% 生成扩展音频
extended = repmat(y, repeat_times, 1);  % 完整复制音频

% 截取到目标长度
if size(extended,1) > target_samples
    extended = extended(1:target_samples, :);
end

% 生成新文件名
[~, name, ext] = fileparts(filename);
new_filename = [name '_extended' ext];

% 写入文件
audiowrite(fullfile(path, new_filename), extended, Fs);
disp(['成功生成：' new_filename]);

% 显示时长对比
fprintf('原始时长：%.2f秒\n新文件时长：%.2f秒\n',...
        original_length/Fs, size(extended,1)/Fs);

%% Audio segmentation for an entire folder --------------------------------
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

% 分割长度
segmentation_length = 10;

% 自动创建输出目录
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

% 获取输入文件夹内所有音频文件（支持wav/mp3格式）
audio_files = dir(fullfile(input_folder, '*.wav'));  % 可扩展为{'*.wav','*.mp3'}

% 创建输出文件夹（如果不存在）
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

% 遍历处理每个音频文件
for i = 1:length(audio_files)
    % 读取音频数据
    [audio, fs] = audioread(fullfile(input_folder, audio_files(i).name));

    % 文件名
    [~, filename, ext] = fileparts(audio_files(i).name);
    disp(filename);
    
    % 计算音频总时长与可分段数
    total_samples = size(audio, 1);
    segment_length = segmentation_length * fs;  % 对应的样本数
    num_segments = floor(total_samples / segment_length);
    
    % 分割音频并写入文件
    [~, base_name, ext] = fileparts(audio_files(i).name);
    for j = 1:num_segments
        start_sample = (j-1)*segment_length + 1;
        end_sample = j*segment_length;
        clip = audio(start_sample:end_sample, :);  % 保持多声道
        
        % 生成带序号后缀的文件名
        new_filename = sprintf('%s-%d%s', base_name, j, ext);
        output_path = fullfile(output_folder, new_filename);
        
        % 写入文件（保持原采样率和比特深度）
        audiowrite(output_path, clip, fs);  
    end
end

winopen(output_folder);




