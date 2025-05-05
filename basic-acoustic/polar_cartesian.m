% 清除工作区和命令窗口
clear all;
clc;

% 输入角度（单位：度）和距离
angle_deg = 180; % 角度
distance = 0.3;   % 距离

% 将角度转换为弧度
angle_rad = deg2rad(angle_deg);

% 计算 x 和 y 坐标
x = distance * cos(angle_rad);
y = distance * sin(angle_rad);

% 显示结果
disp(['对于角度 ' num2str(angle_deg) '° 和距离 ' num2str(distance) ' m:']);
disp(['x 坐标: ' num2str(x)]);
disp(['y 坐标: ' num2str(y)]);

% 可视化结果
figure;
plot(0, 0, 'ro', 'MarkerSize', 10); % 原点
hold on;
plot(x, y, 'bo', 'MarkerSize', 10); % 计算得到的点
line([0 x], [0 y], 'Color', 'k', 'LineStyle', '--'); % 连接原点和计算点的线
grid on;
xlabel('X (m)');
ylabel('Y (m)');
title('根据角度和距离计算 (x, y) 坐标');
axis equal; % 设置坐标轴比例相等
legend('Origin', 'Calculated Point');