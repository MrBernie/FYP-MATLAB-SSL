% function tdoa_angles = computeTDOAAngles(audio_data, fs, num_frames, num_angles, micDist, speedOfSound, frame_samples)
% % 计算每帧的TDOA角度估计
% % 输入：
% %   audio_data  - 双通道音频数据 (Nx2)
% %   fs          - 采样率
% %   num_frames  - 帧数
% %   num_angles  - 角度数量
% %   micDist     - 麦克风间距（米）
% %   speedOfSound- 声速（m/s）
% %   frame_samples - 每帧采样点数
% % 输出：
% %   tdoa_angles - 每帧对应的TDOA角度 (num_frames x 1)
% 
%     tdoa_angles = zeros(num_frames, 1);
%     channel1 = audio_data(:,1);
%     channel2 = audio_data(:,2);
%     angles = linspace(-89, 90, num_angles);
% 
%     for i = 1:num_frames
%         start_idx = (i-1)*frame_samples + 1;
%         end_idx = min(i*frame_samples, length(channel1));
%         frame1 = channel1(start_idx:end_idx);
%         frame2 = channel2(start_idx:end_idx);
% 
%         [corr_values, lags] = xcorr(frame1, frame2);
%         [~, max_idx] = max(corr_values);
%         lag = lags(max_idx);
%         tdoa = lag / fs;
% 
%         tau = (micDist / speedOfSound) * sind(angles);
%         [~, angle_idx] = min(abs(tau - tdoa));
%         tdoa_angles(i) = angles(angle_idx) + 90; % 转换到0-180度范围
%     end
% end
function tdoa_angles = computeTDOAAngles(audio_data, fs, num_frames, num_angles, micDist, speedOfSound, frame_samples)
    tdoa_angles = zeros(num_frames, 1);
    channel1 = audio_data(:,1);
    channel2 = audio_data(:,2);
    angles = linspace(-89, 90, num_angles);
    max_tdoa = micDist / speedOfSound;
    max_lag_samples = round(max_tdoa * fs);

    for i = 1:num_frames
        start_idx = (i-1)*frame_samples + 1;
        end_idx = start_idx + frame_samples - 1;
        if end_idx > length(channel1)
            break; % 或者用零填充
        end
        frame1 = channel1(start_idx:end_idx);
        frame2 = channel2(start_idx:end_idx);

        [corr_values, lags] = xcorr(frame2, frame1);
        % 限制lag范围
        valid_idx = abs(lags) <= max_lag_samples;
        [~, max_idx_rel] = max(corr_values(valid_idx));
        valid_lags = lags(valid_idx);
        lag = valid_lags(max_idx_rel);

        % 亚样本插值（抛物线拟合）
        if max_idx_rel > 1 && max_idx_rel < length(valid_lags)
            y = corr_values(valid_idx);
            x = valid_lags;
            % 取峰值及两侧点
            peak = max_idx_rel;
            denom = 2*(y(peak-1) - 2*y(peak) + y(peak+1));
            if denom ~= 0
                delta = (y(peak-1) - y(peak+1)) / denom;
                lag = lag + delta;
            end
        end

        tdoa = lag / fs;
        tau = (micDist / speedOfSound) * sind(angles);
        [~, angle_idx] = min(abs(tau - tdoa));
        tdoa_angles(i) = angles(angle_idx) + 90;
    end
end
