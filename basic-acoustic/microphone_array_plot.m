
clc
clear all
close all

%% ------------------------------初始化常量-------------------------------%
c = 334;   % 声速c
fs = 16000;   % 抽样频率fs
T = 4;   % ??
t = 0:1/fs:T;  % 时间 [0, 0.1]
L = length(t); % 时间长度:101
f = 500;   % 感兴趣的频率
w = 2*pi*f;  % 角频率
k = w/c;   % 波数 k

%% ------------------------------各阵元坐标-------------------------------%
M = 2;   % 阵元个数

yi = zeros(M,1);
zi = [0; 0];
xi = [-0.08; 0.08];

figure(1)
plot(xi, zi, 'r*');
title('Microphone Array');
xlabel('x (m)');
ylabel('y (m)');
grid on;
xlim([-0.1, 0.1]);
ylim([-0.1, 0.1]);

% 添加坐标标注
for i = 1:M
    text(xi(i), zi(i)+0.01, sprintf('(%.2f, %.2f)', xi(i), zi(i)));
end

