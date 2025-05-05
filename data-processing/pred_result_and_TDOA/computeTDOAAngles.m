function tdoa_angles = computeTDOAAngles(audio_data, fs, num_frames, num_angles, micDist, speedOfSound, frame_samples)
% 计算每帧的TDOA角度估计
% 输入：
%   audio_data  - 双通道音频数据 (Nx2)
%   fs          - 采样率
%   num_frames  - 帧数
%   num_angles  - 角度数量
%   micDist     - 麦克风间距（米）
%   speedOfSound- 声速（m/s）
%   frame_samples - 每帧采样点数
% 输出：
%   tdoa_angles - 每帧对应的TDOA角度 (num_frames x 1)

    tdoa_angles = zeros(num_frames, 1);
    channel1 = audio_data(:,1);
    channel2 = audio_data(:,2);
    angles = linspace(-89, 90, num_angles);

    for i = 1:num_frames
        start_idx = (i-1)*frame_samples + 1;
        end_idx = min(i*frame_samples, length(channel1));
        frame1 = channel1(start_idx:end_idx);
        frame2 = channel2(start_idx:end_idx);

        [corr_values, lags] = xcorr(frame1, frame2);
        [~, max_idx] = max(corr_values);
        lag = lags(max_idx);
        tdoa = lag / fs;

        tau = (micDist / speedOfSound) * sind(angles);
        [~, angle_idx] = min(abs(tau - tdoa));
        tdoa_angles(i) = angles(angle_idx) + 90; % 转换到0-180度范围
    end
end
