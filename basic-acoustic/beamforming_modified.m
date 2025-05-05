clear all;
close all;

%% --------------------------------Read File-------------------------------%

% 读取音频文件
[fileName, filePath] = uigetfile('*.wav', '选择音频文件');
audioFile = fullfile(filePath, fileName);
[p, fs] = audioread(audioFile); % 读取多通道音频文件
p = double(p);
dt = 1/fs;
t = (0:length(p)-1)/fs;
t = transpose(t);

% %% ------------------------------Plot Original Sygnal-------------------------------%
figure;
title('Original Sound Signal');
tiledlayout(6,1)
for i = 1 : 6
    nexttile();
    plot(t, p(:,i));  
end

%% -------------------------------Frequency filter----------------------------------%
% 设计带通滤波器
low_cutoff = 1900; % 下限频率
high_cutoff = 2100; % 上限频率
f = 2000; % 关心的频率

% bpFilt = designfilt('bandpassfir', 'FilterOrder', 20, ...
%                     'CutoffFrequency1', low_cutoff, 'CutoffFrequency2', high_cutoff, ...
%                     'SampleRate', fs);
% p = filter(bpFilt,p);

% 使用巴特沃斯滤波器设计函数
[b, a] = butter(4, [low_cutoff high_cutoff] / (fs / 2), 'bandpass');

% 应用巴特沃斯滤波器
p = filter(b, a, p);

% pf = fft(p);


%% ------------------------------Plot Filtered Sygnal-------------------------------%
figure;
title('Filtered Sound Signal');
tiledlayout(6,1)
for i = 1 : 6
    nexttile();
    plot(t, p(:,i));  
end
%% ------------------------------Coordinate of each Mic-------------------------------%
zi = zeros(6,1); % 生成一个M*1维的零矩阵
% yi = [0.06;  0.03; -0.03; -0.06; -0.03;  0.03];
% xi = [   0; 0.052; 0.052;     0;-0.052;-0.052];

yi = [    0; 0.052; 0.052;    0;-0.052;-0.052]*1;
xi = [-0.06; -0.03; 0.03; 0.06;  0.03; -0.03]*1;


figure;
plot(xi,yi,'r*');
grid on;
% 在每个数据点旁边添加坐标标签
for i = 1:length(xi)
    % 使用 text 函数添加坐标，x(i) 和 y(i) 是点的位置
    text(xi(i), yi(i), sprintf('M%d: (%0.2f, %0.2f)', i, xi(i), yi(i)), ...
        'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
end
title('Microphone Array');
xlim([-0.1 0.1]);
ylim([-0.1 0.1]);
xlabel('x (m)');
ylabel('y (m)');
%% ------------------------------Normalize original signal-------------------------------%
p = transpose(p);

pn_max = zeros(size(p, 1), 1);
for i = 1:size(p,1)
    pn_max = max(p(i,:)); % y的第i行的最大元素的值
    pn(i,:) = p(i,:)/pn_max;
end

% for i = 1:size(p,1)
%     pn_mean = mean(p(i,:)); % y的第i行的最大元素的值
%     pn(i,:) = p(i,:)/pn_mean;
% end

% pn=p;

%% ------------------------------Beamforming-------------------------------%
T = length(p)/fs;
R = pn*pn'/T; % 接收数据的自协方差矩阵  A.'是一般转置，A'是共轭转置
w = 2*pi*f;  % 角频率
c = 334; % 声速

x2 = 0;
y2 = 0;
z2 = 0;

% 我们设置步长为0.1，扫描范围是20x20的平面，双重for循环得到M*1矢量矩阵，最后得到交叉谱矩阵（cross spectrum matrix）
% 由DSP理论，这个就是声音的功率。
step_x = 0.001;  % 步长设置为0.1
step_y = 0.001;
d_x = 1;
d_y = 1;
y = (-1*d_x:step_y:d_x);
x = (-1*d_x:step_x:d_x);  % 扫描范围
z = 0;

for k1=1:length(y)
    for k2=1:length(x)
        Ri = sqrt((x(k2)-xi).^2+(y(k1)-yi).^2+(z-zi).^2);  % 该扫描点到各阵元的聚焦距离矢量
        Ri2 = sqrt((x(k2)-x2).^2+(y(k1)-y2).^2+(z-z2).^2);
        Rn = Ri-Ri2;   % 扫描点到各阵元与参考阵元的程差矢量
        b = exp(-j*w*Rn/c); % 声压聚焦方向矢量
        Pcbf(k1,k2) = abs(b'*R*b); % CSM
    end
end

%% -------------------------------------Normalization-------------------------------------%
for k1 = 1:length(y)
    pp(k1) = max(Pcbf(k1,:)); % Pcbf 的第k1行的最大元素的值
end

Pcbf = Pcbf/max(pp);
[maxValue, linearIndex] = max(Pcbf, [], 'all');
[row, column] = ind2sub(size(Pcbf), linearIndex);
y_max = row*step_x-d_x;
x_max = column*step_y-d_y;
disp(['Maximum number ', num2str(maxValue), ', location:(', num2str(x_max), ',', num2str(y_max), ')']);
disp(['Index of the maximum ', num2str(column), ',', num2str(row)])
% 计算角度（弧度）
theta_rad = atan2(y_max, x_max);
% 转换为度
theta_deg = rad2deg(theta_rad);

disp(['Degree of Arrival ', num2str(theta_deg)])


%% -------------------------------------Plot------------------------------------%
figure;
surf(x,y,Pcbf, 'EdgeColor','none');
xlabel('x(m)'),ylabel('y(m)')
title('3D Beamforming')
colorbar

 
figure;
pcolor(x,y,Pcbf);
shading interp;
xlabel('x(m)');
ylabel('y(m)');
title('Beamforming result')
colorbar