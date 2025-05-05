% 清理环境
clear;
clc;

% % 创建数据采集会话
s = daq('ni');
% 
% % 添加模拟输入通道，假设NI设备ID为'Dev1'，通道为'ai0'，类型为'Voltage'
channel_1 = addinput(s, 'cDAQ1Mod1', 'ai0', 'IEPE');
channel_2 = addinput(s, 'cDAQ1Mod1', 'ai1', 'IEPE');
% % 配置采样率
s.Rate = 51200; % 采样率 51200 Hz，适用于音频信号
duration = 10; % 持续时间设置为5秒
fs = s.Rate;
% fs = 51200;

                                                                                                                                                                                                              
% 计划读取的时间长度
dataTime = seconds(duration);

% 开始采集数据
[data, time] = read(s, dataTime, 'OutputFormat', 'Matrix');
ch_1 = data(:,1);
ch_2 = data(:,2);

%%%%%%%%%%%%%%
data = data / max(abs(data(:))); % normalization

% Save and Load audio data
audiowrite("180d-0d-2m-5(10s).wav", data, fs)
% audiowrite("NI_recording.wav", data, fs)
% data = audioread("exp_SPL.wav");
% ch_1 = data(:,1);
% ch_2 = data(:,2);

channel1_fft = fft(ch_1);
channel2_fft = fft(ch_2);

N = size(ch_1, 1); % Length of the audio data
frequencies = (0:(N/2))*(fs/N); % Frequency vector

% validate with the COMSOL
% load rtf_power.txt
% Freq=rtf_power(:,1);
% point1=rtf_power(:,2);
% 
% load ss_spl.txt
% Freq_ss = ss_spl(:,1);
% point2 = ss_spl(:,2);

% figure
% plot(frequencies, smooth(20*log10(channel1_fft(1:N/2+1)./channel2_fft(1:N/2+1)),50),'r-','LineWidth',2);
% plot(Freq_ss, point2)



% % Plot time domain signals
% % figure
% % subplot(2, 1, 1);
% % plot((1:N)/fs, ch_1);
% % xlabel('Time (s)');
% % ylabel('Channel 1' );
% % title('Channel 1 - Time Domain');
% % subplot(2, 1, 2);
% % plot((1:N)/fs, ch_2);
% % xlabel('Time (s)');
% % ylabel('Channel 2' );
% % title('Channel 2 - Time Domain');
% % Plot frequency domain signals
% % figure
% % plot(frequencies, 20*log10(channel1_fft(1:N/2+1)./channel2_fft(1:N/2+1)));
% % hold on
% % plot(frequencies, smooth(20*log10(channel1_fft(1:N/2+1)./channel2_fft(1:N/2+1)),50),'r-','LineWidth',2);
% % plot(frequencies, 20*log10(channel1_fft(1:N/2+1)));
% % plot(frequencies, smooth(20*log10(channel1_fft(1:N/2+1)),50), "r-", "LineWidth",2);
% % xlabel('Frequency (Hz)');ylabel('Channel 1 (dB)' );
% % legend('EXP-original','EXP-smooth')
% % title('Channel 1- Frequency Domain');
% % xlim([10 2000])
% % 
% % Plot frequency domain signals
% % figure
% % plot(frequencies, 20*log10(channel2_fft(1:N/2+1))+80);
% % hold on
% % A = smooth(20*log10(channel2_fft(1:N/2+1))+80,50);
% % plot(frequencies, smooth(20*log10(channel2_fft(1:N/2+1))+80,50),'r-','LineWidth',2);
% % plot(Freq, point1, "k-","LineWidth",2);
% % xlabel('Frequency (Hz)');ylabel('Channel 2 (dB)' );
% % legend('EXP-original','EXP-smooth',"The proposed method")
% % title('Channel 2- Frequency Domain');
% % xlim([10 2000])

% 绘制时间-幅值图
figure;

% 通道1
subplot(2, 1, 1);
plot(time, ch_1, 'b', 'LineWidth',0.8);
title('Time-Amplitude Graph for Channel 1');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;
legend('Channel 1');

% 通道2
subplot(2, 1, 2);
plot(time, ch_2, 'r', 'LineWidth',0.8);
title('Time-Amplitude Graph for Channel 2');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;
legend('Channel 2');

% 调整布局
sgtitle('Time-Amplitude Graphs for Both Channels');
