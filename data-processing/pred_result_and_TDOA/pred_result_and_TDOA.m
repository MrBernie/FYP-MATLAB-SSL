clear all;
close all;

%% Parameter Settings (User Configurable)
micDist = 0.16;          % Microphone spacing (meters)
speedOfSound = 343;     % Speed of sound (m/s)
total_time = 15;        % Total duration (seconds) 

% file name also depends on these variables!!!
ground_truth = 75;     % Manual ground truth angle (0-180) 
distance = 2;
file_index_same_pos = 1;
file_index_time_length = '(60s)';
file_index_segmentation = 1;

distance = num2str(distance);
file_index_same_pos = num2str(file_index_same_pos);
file_index_segmentation = num2str(file_index_segmentation);

% %% File Selection
% try
%     [audio_file, audio_path] = uigetfile('*.wav', 'Select raw audio file');
%     if audio_path == 0
%         error('User canceled input folder selection');
%     end
% 
%     [pred_file, pred_path] = uigetfile('*.txt', 'Select prediction result file');
%     if pred_path == 0
%         error('User canceled output folder selection');
%     end
% catch ME
%     errordlg(['Path selection error: ' ME.message], 'System Error');
%     return;
% end
% vad_dir = uigetdir('', 'Select VAD file directory');

% Select the root folder
% disp(num2str(ground_truth));
root_dir = 'D:\OneDrive\OneDrive - The Hong Kong Polytechnic University\PolyU Folder\FYP\experiment\';
% audio_path = 'D:\OneDrive\OneDrive - The Hong Kong Polytechnic University\PolyU Folder\FYP\experiment\data-mar-18-stationary-segmented-filtered\';
% audio_path = 'data-mar-18-stationary-segmented-filtered';
audio_path = 'data-mar-18-stationary-segmented-filtered';
audio_path = fullfile(root_dir,audio_path);
audio_file = [num2str(ground_truth), 'd-',distance,'m-',file_index_same_pos,file_index_time_length,'-',file_index_segmentation,'-bandpass-filtered.wav'];
disp(audio_file);

% pred_path = 'D:\OneDrive\OneDrive - The Hong Kong Polytechnic University\PolyU Folder\FYP\experiment\pred_results_mar_18_on_mar_20_stationary\pred_matrix\';
prediction_data_path = 'pred_results_mar_18_on_apr_8';
pred_path = fullfile(root_dir,prediction_data_path,'\pred_matrix');
pred_file = [num2str(ground_truth), 'd-',distance,'m-',file_index_same_pos,file_index_time_length,'-',file_index_segmentation,'-bandpass-filtered-pred-matrix.txt'];

% vad_path = 'D:\OneDrive\OneDrive - The Hong Kong Polytechnic University\PolyU Folder\FYP\experiment\pred_results_mar_18_on_mar_20_stationary\vad_out\';
vad_path = fullfile(root_dir,prediction_data_path,'\vad_out');
vad_file =[num2str(ground_truth), 'd-',distance,'m-',file_index_same_pos,file_index_time_length,'-',file_index_segmentation,'-bandpass-filtered-vad-out.txt'];

%% Data Loading
[audio_data, fs] = audioread(fullfile(audio_path, audio_file));
pred_matrix = readmatrix(fullfile(pred_path, pred_file));
[num_frames, num_angles] = size(pred_matrix);

%% VAD File Handling
% [~, pred_name] = fileparts(pred_file);
% vad_name = strrep(pred_name, '-pred-matrix', '-vad-out');
vad_data = readmatrix(fullfile(vad_path, vad_file));

%% Time Parameter Calculation
frame_duration = total_time / num_frames;
frame_samples = round(frame_duration * fs);
num_frames_audio = floor(size(audio_data,1) / frame_samples);
time_axis = (0:num_frames-1)*frame_duration + frame_duration/2;

%% Array Initialization
tdoa_angles = zeros(num_frames, 1);
dl_angles = zeros(num_frames, 1);
avg_dl_angles = zeros(num_frames,1);
avg_dl_angles_cut = zeros(num_frames,1);

% choose a cutting value between [0 1]
% choosing value -1 cut the data from raw 0
cut_value = 0.8; 

pred_matrix_norm = zeros(num_frames,num_angles);
avg_vad = zeros(num_frames, 1);

%% TDOA Calculation
channel1 = audio_data(:,1);
channel2 = audio_data(:,2);

for i = 1:num_frames
    % Frame segmentation processing
    start_idx = (i-1)*frame_samples + 1;
    end_idx = min(i*frame_samples, length(channel1));

    frame1 = channel1(start_idx:end_idx);
    frame2 = channel2(start_idx:end_idx);

    % Calculate TDOA
    [corr_values, lags] = xcorr(frame1, frame2);
    [~, max_idx] = max(corr_values);
    lag = lags(max_idx);
    tdoa = lag / fs;

    % Angle estimation
    angles = linspace(-89, 90, num_angles);
    tau = (micDist / speedOfSound) * sind(angles);
    [~, angle_idx] = min(abs(tau - tdoa));
    tdoa_angles(i) = angles(angle_idx) + 90; % Convert to 0-180 range
end

%% Deep Learning Prediction Processing
for i = 1:num_frames
    frame_pred = pred_matrix(i,:);
    
    % Directly obtain global maximum
    [~, max_idx] = max(frame_pred);
    % disp(max_idx);
    dl_angles(i) = max_idx;
    % dl_angles(i) = (max_idx - 1) * (180/(num_angles-1)); % More accurate angle conversion

    % normalized the prediction matrix to [0,1]
    min_w = min(frame_pred);
    max_w = max(frame_pred);
    frame_pred_norm = (frame_pred - min_w)/(max_w - min_w);
    zero_w = (0 - min_w)/(max_w - min_w);
    pred_matrix_norm(i,:) = frame_pred_norm;

    % Get the average predicted angle by min&max normalization 
    weights_sum_norm = frame_pred_norm / sum(frame_pred_norm);
    angle_index = 1:num_angles;
    avg_dl_angles(i) = sum(weights_sum_norm .* angle_index);

    % Get the average predicted angle by cutting off negative weight
    if cut_value == -1
        cut_value = zero_w;
    end
    filtered_pred_frame = frame_pred_norm(frame_pred_norm > cut_value) .* angle_index(frame_pred_norm > cut_value);
    filtered_pred_frame_sum = sum(filtered_pred_frame);
    filtered_pred_frame_weight_sum = sum(frame_pred_norm(frame_pred_norm > cut_value));
    if ~isempty(filtered_pred_frame)
        avg_dl_angles_cut(i) = (filtered_pred_frame_sum)/filtered_pred_frame_weight_sum;
    else
        average = 0; % cant find > 0
    end
end

%% VAD Processing
samples_per_frame = size(vad_data,1) / num_frames;
for i = 1:num_frames
    start_sample = round((i-1)*samples_per_frame + 1);
    end_sample = round(i*samples_per_frame);
    avg_vad(i) = mean(vad_data(start_sample:end_sample));
end

%% Visualization
% % figure('Position', [000 000 800 600], 'Color', 'white');
% figure;
% 
% % Angle Comparison
% % subplot(3,1,1); % with VAD
% subplot(2,1,1);
% hold on;
% scatter(time_axis, tdoa_angles, 40, 'r', 'filled', 'MarkerEdgeColor', 'k');
% scatter(time_axis, dl_angles, 40, 'b', 'filled', 'MarkerEdgeColor', 'k');
% plot([0 total_time], [ground_truth ground_truth], 'g--', 'LineWidth', 2);
% title([pred_file ' - Angle Prediction Comparison']);
% xlabel('Time (seconds)');
% ylabel('Angle (°)');
% legend({'TDOA Prediction', 'Deep Learning Prediction', 'Ground Truth'}, 'Location', 'best');
% ylim([0 180]);
% grid on;
% 
% % Original Audio Waveform
% % subplot(3,1,2); % with VAD
% subplot(2,1,2);
% t_audio = (0:length(audio_data)-1)/fs;
% plot(t_audio, audio_data(:,1), 'b');
% hold on;
% plot(t_audio, audio_data(:,2), 'r');
% title('Original Audio Signal');
% xlabel('Time (seconds)');
% ylabel('Amplitude');
% legend('Channel 1', 'Channel 2');
% xlim([0 total_time]);
% grid on;
% 
% % % VAD Results
% % subplot(3,1,3);
% % plot(time_axis, avg_vad, 'k', 'LineWidth', 1.5);
% % title('Voice Activity Detection (VAD)');
% % xlabel('Time (seconds)');
% % ylabel('VAD Average');
% % ylim([0 1]);
% % grid on;
% 
% % Synchronize zoom
% linkaxes(findall(gcf,'type','axes'), 'x');

%% Visualization 2
figure;

% % Angle Comparison with Prediction Matrix Heatmap
% subplot(2,1,1);
hold on;

% Generate angle axis (discrete angles for deep learning predictions)
angles_dl = linspace(1, 180, num_angles);

% Plot prediction matrix as a heatmap (with semi-transparent colormap)

% customized_map = colormap(spring);
% customized_map(end, :) = [1, 1, 1];
% colormap(customized_map);
% colormap(spring); % Use a high-contrast colormap
% alpha(0.35); % Set transparency for the heatmap
% caxis([0 1]); % Fix the color range (assuming predictions are probability values)

% red_green_map = zeros(256, 3);
% for i = 1:256
%     if i <= 128
%         red_green_map(i, :) = [(128-i)/128, 0, 0]; % red
%     else
%         red_green_map(i, :) = [0, (i-128)/128, 0]; % green
%     end
% end
% colormap(red_green_map);

colormap;
colorbar; % Add a colorbar to show prediction value scale
pcolor(time_axis, angles_dl, pred_matrix'); % Note: Transpose the matrix to match axis orientation

% Keep the original scatter plots
scatter(time_axis, tdoa_angles, 40, 'r', 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha',0.9);
scatter(time_axis, dl_angles, 40, 'w', 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha',0.9);
% scatter(time_axis, avg_dl_angles_cut, 40, 'g', 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha',0.9);

plot([0 total_time], [ground_truth ground_truth], 'g--', 'LineWidth', 2, 'Color', [0 1 0]);
grid on;

xlim([0 total_time]);

% Enhance plot aesthetics
title([strrep(audio_file,'_','\_') ' - Prediction Comparison']);
xlabel('Time (seconds)');
ylabel('Angle (°)');
legend({'Deep Learning Prediction Result', 'TDOA Prediction', 'Deep Learning Prediction Max','Ground Truth'}, 'Location', 'best');
ylim([0 180]);
grid on;
set(gca, 'GridLineStyle', ':', 'GridColor', 'k', 'GridLineWidth', 0.1);
set(gca, 'Layer','top'); % Ensure grid and labels are displayed above the heatmap


