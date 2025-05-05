function expected_angles = computeExpectedAngle(pred_matrix, smooth)
    % computeAngleError 计算每个时间帧的期望角度误差
    %
    % 输入：
    %   pred_matrix - [num_frames x num_angles] 预测概率矩阵，每行是该时间帧对所有角度的预测概率
    %   DOA_gt      - [num_frames x 1] 该时间帧的真实角度（标量或向量）
    %
    % 输出：
    %   angle_errors - [num_frames x 1] 每个时间帧的期望角度与GT角度的差值（误差）
    %
    % 说明：
    %   期望角度 = sum(角度值 * 预测概率)
    %   误差 = 期望角度 - GT角度
    
    [num_frames, num_angles] = size(pred_matrix);
    
    angles = linspace(0, 180, num_angles);
    
    % 对每行归一化，避免除零
    row_sums = sum(pred_matrix, 2);
    row_sums(row_sums == 0) = 1; % 防止除零
    
    pred_norm = pred_matrix ./ row_sums;
    
    expected_angles = pred_norm * angles';
    
    % angle_errors = expected_angles - DOA_gt(:);
    if smooth
        expected_angles = smoothdata(expected_angles,"gaussian");
    end


end


