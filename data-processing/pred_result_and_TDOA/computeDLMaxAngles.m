function dl_angles = computeDLMaxAngles(pred_matrix, smooth)
% 计算每帧深度学习预测的最大概率对应的角度索引
% 输入：
%   pred_matrix - 预测矩阵 (num_frames x num_angles)
% 输出：
%   dl_angles   - 每帧最大概率对应的角度索引 (num_frames x 1)

    num_frames = size(pred_matrix,1);
    dl_angles = zeros(num_frames,1);
    for i = 1:num_frames
        [~, max_idx] = max(pred_matrix(i,:));
        dl_angles(i) = max_idx;
    end

    if smooth
        dl_angles = smoothdata(dl_angles,"gaussian");
    end

end
