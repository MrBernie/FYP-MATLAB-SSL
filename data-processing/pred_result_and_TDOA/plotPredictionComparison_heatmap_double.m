function fig = plotPredictionComparison_heatmap_double(time_axis, plotData, audio_file, plotOptions)
% plotPredictionComparison - 绘制预测结果对比图(上下两图)
%
% 输入参数：
%   time_axis      - 时间轴（秒）
%   plotData       - 结构体:
%                      plotData.angles          - 角度轴（0~180度）
%                      plotData.pred_matrix     - 深度学习预测概率矩阵（num_frames x num_angles）
%                      plotData.tdoa_angles     - TDOA预测角度（num_frames x 1）
%                      plotData.DOA_gt          - 真实值（num_frames x 1）
%                      plotData.dl_angles       - 深度学习预测最大角度索引（num_frames x 1）
%                      plotData.total_time      - 音频总时长（秒）
%                      plotData.expected_angles - 深度学习预测期望角度（num_frames x 1）
%   audio_file     - 音频文件名（字符串）
%   output_dir     - 输出文件夹路径
%   file_prefix    - 文件名前缀
%   plotOptions    - 结构体:
%                    plotOptions.showPredMatrix (bool)
%                    plotOptions.showTDOA (bool)
%                    plotOptions.showDLMax (bool)
%                    plotOptions.showGT (bool)
%                    plotOptions.showExpectedAngles (bool)
%                    plotOptions.saveFig (bool)
%                    
%
% 功能：
%   根据plotOptions绘制对应图形，返回fig。
    
    fig = figure('Position', [100 100 800 1200]);
    t = tiledlayout(2,1, 'TileSpacing', 'compact', 'Padding', 'compact');
    
    % 上图：预测概率热图 + 真实值
    % ax1 = subplot(2,1,1);
    ax1 = nexttile(1);
    hold on;
    
    if plotOptions.showPredMatrix
        imagesc(time_axis, plotData.angles, plotData.pred_matrix');
        axis xy;
        colormap(sky);
        colorbar;
    end
    
    if plotOptions.showGT && ~isempty(plotData.DOA_gt)
        plot(time_axis + mean(diff(time_axis))/2, plotData.DOA_gt, 'g--', 'LineWidth', 2);
    end
    
    xlim([0 plotData.total_time]);
    ylim([0 180]);
    title([strrep(audio_file, '_', '\_') ' - Prediction Heat Map & Ground Truth']);
    xlabel('Time (seconds)');
    ylabel('Angle (°)');
    grid on;
    
    legendEntriesTop = {};
    % if plotOptions.showPredMatrix
    %     legendEntriesTop{end+1} = 'Deep Learning Prediction Heatmap';
    % end
    if plotOptions.showGT && ~isempty(plotData.DOA_gt)
        legendEntriesTop{end+1} = 'Ground Truth';
    end
    if ~isempty(legendEntriesTop)
        legend(legendEntriesTop, 'Location', 'best');
    end
    
    hold off;
    
    % 下图：真实值 + 其他点（TDOA、DLMax、ExpectedAngles）
    % ax2 = subplot(2,1,2);
    ax2 = nexttile(2);
    hold on;
    
    if plotOptions.showGT && ~isempty(plotData.DOA_gt)
        plot(time_axis + mean(diff(time_axis))/2, plotData.DOA_gt, 'g--', 'LineWidth', 2);
    end
    
    if plotOptions.showTDOA && ~isempty(plotData.tdoa_angles)
        scatter(time_axis + mean(diff(time_axis))/2, plotData.tdoa_angles, ...
            5, 'w', 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', 0.9);
    end
    
    if plotOptions.showDLMax && ~isempty(plotData.dl_angles)
        % scatter(time_axis + mean(diff(time_axis))/2, plotData.dl_angles, ...
            % 25, 'r', 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', 0.9);
        plot(time_axis + mean(diff(time_axis))/2, plotData.dl_angles, 'r-', 'LineWidth', 1); 
    end
    
    if plotOptions.showExpectedAngles && ~isempty(plotData.expected_angles)
        % scatter(time_axis + mean(diff(time_axis))/2, plotData.expected_angles, ...
            % 25, 'b', 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', 0.9);
        plot(time_axis + mean(diff(time_axis))/2, plotData.expected_angles, 'b-', 'LineWidth', 1); 
    end
    
    xlim([0 plotData.total_time]);
    ylim([0 180]);
    title('Ground Truth and Other Predictions');
    xlabel('Time (seconds)');
    ylabel('Angle (°)');
    grid on;
    
    legendEntriesBottom = {};
    if plotOptions.showGT && ~isempty(plotData.DOA_gt)
        legendEntriesBottom{end+1} = 'Ground Truth';
    end
    if plotOptions.showTDOA && ~isempty(plotData.tdoa_angles)
        legendEntriesBottom{end+1} = 'TDOA Prediction';
    end
    if plotOptions.showDLMax && ~isempty(plotData.dl_angles)
        legendEntriesBottom{end+1} = 'Deep Learning Prediction Max Angles';
    end
    if plotOptions.showExpectedAngles && ~isempty(plotData.expected_angles)
        legendEntriesBottom{end+1} = 'Deep Learning Prediction Expected Angles';
    end
    if ~isempty(legendEntriesBottom)
        legend(legendEntriesBottom, 'Location', 'best');
    end
    
    % 简单链接两个子图的x轴，使得缩放和平移同步
    linkaxes([ax1, ax2], 'x');

end
