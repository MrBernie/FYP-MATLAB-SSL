function [aryData,colTime] = ...
    get_audioData_store(inxOS,fs,szFrame,nMic,nBit,nSample,duration)
%% - <> SUMMARY <> -
%
% Record the audio from the readed device and store it in the output array.
% This function treats one frame as one sample.
%
% For more info, please refer to the MATLAB Doc:
% + audioDeviceReader
% 
% [IN]
% inxOS - Operating System index. 1-macOS, 2-Windows.
% fs - Sampling frequency for device.
% szFrame - Sample per frame.
% nMic - Mic number for device.
% nBit - Bit Depth for device
% nSample - Sample number.
% duration - Sampling duration for the output .wav file.
%
% [OUT]
% aryData - 3D data array with the 3rd dimension is the sample number.
% colTime - Time sequence in colunm.
%
% [UPDATE]
% May/15/2024 @LinfengLI
% May/31/2024 @LinfengLI
%
%% - <> MAIN <> -

%% ## DEVICE ##
% Refresh the list of available devices after an audio device has been
% added to or removed from the machine.
audiodevreset

% List available audio devices.
objDev = audioDeviceReader;
nameDev = getAudioDevices(objDev);

fprintf('\n # Available Devices #\n')
for iName = 1:length(nameDev)
    fprintf(' > Device %i: %s\n',iName,nameDev{iName})
end

% Select audio device.
if inxOS == 1 % macOS
    objDev.Device = 'YDM8MIC Audio';
elseif inxOS == 2 % windows
    objDev.Driver = 'DirectSound';
    objDev.Device = '麦克风阵列 (YDM6MIC Audio)';
else
    error('Assigned Operating System is not supported currently.')
end
fprintf('[*] Device %s is selected.\n',objDev.Device)

% Parameters
objDev.SampleRate = fs;
objDev.SamplesPerFrame = szFrame;
objDev.ChannelMappingSource = 'Auto';
objDev.NumChannels = nMic;
objDev.OutputDataType = 'double';

if nBit == 16
    objDev.BitDepth = '16-bit integer';
else
    error('The BitDepth set is not supported currently.')
end

% [MONITOR]
fprintf('\n # Device Properties #\n')
disp(objDev);

%% ## Acquire DATA ##
% Time sequence
colTime = transpose((0:szFrame-1)*(1/fs));

% Allocate space
aryData = zeros(szFrame,nMic,nSample);

% Store data
fprintf('\n## Speak into silicon mic now\n')

for iSample = 1:nSample
    aryData(:,:,iSample) = objDev();
end

fprintf('\n Recording Complete. ##\n')

%% ## Acquire .wav file ##
% Set filename
filename = sprintf('audio_%iMic_%.1fsec_%s.wav',nMic,duration, ...
    string(datetime('now','Format','yyMMdd_HHmm')) );

% Specify the file name and type to write
%objWav = dsp.AudioFileWriter(filename,'FileFormat','WAV');

% tic
% while toc < duration
%     mxFrame = objDev();
%     objWav(mxFrame);
% end

%% ## Release ##
% Release device
release(objDev)

% Release Writer
%release(objWav)

end % END OF FUNC


