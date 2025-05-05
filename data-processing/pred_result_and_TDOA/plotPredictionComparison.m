function fig = plotPredictionComparison(time_axis, plotData, audio_file, plotOptions)
% plotPredictionComparison - 绘制预测结果对比图
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

    


    % fig = figure('Visible', 'off', 'Position', [100 100 1280 720]);
    fig = figure('Position', [100 100 800 600]);
    hold on;

    % 画深度学习预测概率热图
    if plotOptions.showPredMatrix
        imagesc(time_axis, plotData.angles, plotData.pred_matrix');
        axis xy;
        % shading interp;
        % shading flat;
        % colormap(jet);
        % colormap(hot);
        colormap(sky);
        % c = flipud(gray);
        % colormap(c);
        % colorbar;
        grid off;
        colorbar;
        % set(gca, 'GridLineStyle', ':', 'GridColor', 'k', 'GridLineWidth', 0.05);
        % set(gca, 'Layer', 'top');
    end

    % 画TDOA预测点
    if plotOptions.showTDOA
        scatter(time_axis + mean(diff(time_axis))/2, plotData.tdoa_angles, ...
            10, 'w', 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', 0.9);
    end

    % 画深度学习预测最大值点
    if plotOptions.showDLMax
        scatter(time_axis + mean(diff(time_axis))/2, plotData.dl_angles, ...
            10, 'r', 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', 0.9);
    end
    
    % 画深度学习预测期望角度
    if plotOptions.showExpectedAngles
        scatter(time_axis + mean(diff(time_axis))/2, plotData.expected_angles, ...
        10, 'b', 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', 0.9);
    end
    
    % 画真实值
    if plotOptions.showGT && ~isempty(plotData.DOA_gt)
        plot(time_axis + mean(diff(time_axis))/2, plotData.DOA_gt, 'g--', 'LineWidth', 2);
    end

    xlim([0 plotData.total_time]);
    ylim([0 180]);
    title([strrep(audio_file, '_', '\_') ' - Prediction Heat Map']);
    xlabel('Time (seconds)');
    ylabel('Angle (°)');

    legendEntries = {};
    % if plotOptions.showPredMatrix
    %     legendEntries{end+1} = 'Deep Learning Prediction Result';
    % end
    if plotOptions.showTDOA
        legendEntries{end+1} = 'TDOA Prediction';
    end
    if plotOptions.showDLMax
        legendEntries{end+1} = 'Deep Learning Prediction Max';
    end
    if plotOptions.showExpectedAngles
        legendEntries{end+1} = 'Deep Learning Prediction Expected Angles';
    end
    if plotOptions.showGT
        legendEntries{end+1} = 'Ground Truth';
    end
    if ~isempty(legendEntries)
        legend(legendEntries, 'Location', 'best');
    end

end








