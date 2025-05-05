function [objTim] = ...
    get_audioData_real(inxOS,fs,szFrame,nMic,nBit,delayStart,period)
%% - <> SUMMARY <> -
%
% This is a timer callback function.
% This function is only valid for the Silicon Mic Array.
%
% For more details, please refer to the MATLAB Doc:
% + timer
% + timer callback function
% 
% [IN]
% inxOS - Operating System index. 1-macOS, 2-Windows.
% fs - Sampling frequency for device.
% szFrame - Sample per frame.
% nMic - Mic number for device.
% nBit - Bit Depth for device
% delayStart - Delay between start of timer and first execution [sec].
% period - Delay between executions in timer [sec].
%
% [OUT]
% objTim - Timer object
%
% [UPDATE]
% May/15/2024 @LinfengLI
% May/16/2024 @LinfengLI
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
    objDev.Device = 'Microphone Array (3- YDM8MIC Audio)';
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

%% ## TIMER ##

colTime = (0:szFrame-1)*(1/fs);

% # Timer Info #
% Create timer object
objTim = timer;

% Timer start callback function
objTim.StartFcn = @func_start;

% Timer stop callback function
objTim.StopFcn = {@func_stop,objDev};

% Timer callback function
objTim.TimerFcn = {@func_record, objDev, nMic, colTime, szFrame};

objTim.Period = period;
objTim.ExecutionMode = 'fixedRate';
% valid only when the `ExecutionMode` is set to `fixedRate`
objTim.BusyMode = 'drop'; 

% Start Delay [sec]
objTim.StartDelay = delayStart;

% # Plot #
fig = figure;
fig.Name = 'REAL-TIME MONITOR';
fig.Units = 'centimeters';
% left bottom
fig.Position(1:2) = [2,1]*5;
% width height
fig.Position(3:4) = [3,4]*5;

end

%% - <> SUB FUNCTIONS <> -
% [ NOTE ]
% In following sub-functions, (~,~) indicates timer object and event, respectively.

% Timer callback function
function func_record(~,~,objDev,nMic,colTime,nt)
    % one frame of audio samples.
    mxFrame = objDev();
    % FFT
    mxFFT = (1/nt).*fft(mxFrame,[],1);
    nf = nt/2;
    mxFFT = mxFFT(1:nf,:);
    mxFFT(2:end,:) = 2.* mxFFT(2:end,:);
    mxFFT_abs = abs(mxFFT);
    colFreq = transpose( (0:nf-1)*(1/colTime(end)) );
    % Draw
    for iMic = 1:nMic
        % time history
        subplot(nMic,2,2*iMic-1)
        plot(colTime,mxFrame(:,iMic),'LineWidth',0.75)
        if iMic == 1
            title(sprintf('Time History Mic %i',iMic),'FontSize',8,'FontWeight','normal')
        elseif iMic == nMic
            title(sprintf('Mic %i',iMic),'FontSize',8,'FontWeight','normal')
            xlabel('Time/(sec)')
        elseif iMic == round(nMic/2)
            title(sprintf('Mic %i',iMic),'FontSize',8,'FontWeight','normal')
            ylabel('Acoustic Voltage/(V)')
        else
            title(sprintf('Mic %i',iMic),'FontSize',8,'FontWeight','normal')
        end
        grid on
        ax = gca;
        ax.XLim = [0.0,1.0];
        ax.XLimitMethod = 'padded';
        
        % frequency spectrum
        subplot(nMic,2,2*iMic)
        plot(colFreq/1e3,mxFFT_abs(:,iMic),'LineWidth',0.75)
        if iMic == 1
            title(sprintf('FFT Spectrum Mic %i',iMic),'FontSize',8,'FontWeight','normal')
        elseif iMic == nMic
            title(sprintf('Mic %i',iMic),'FontSize',8,'FontWeight','normal')
            xlabel('Frequency/(kHz)')
        elseif iMic == round(nMic/2)
            title(sprintf('Mic %i',iMic),'FontSize',8,'FontWeight','normal')
            ylabel('Acoustic Voltage/(V)')
        else
            title(sprintf('Mic %i',iMic),'FontSize',8,'FontWeight','normal')
        end
        grid on
        ax = gca; ax.XLim = [0.01,10];
        ax.XScale = 'linear';
        ax.XLimitMethod = 'padded';
    end
    % update the plot        
    drawnow
end

% Timer start callback function
function func_start(~,~)
    fprintf('\n ## Start Timer ##\n')
    fprintf(' # NOTE\n')
    fprintf([' If `stop(timer)` has no respond, click PAUSE and STOP\n' ...
             ' in the EDITOR, and type `stop(timer)` in the COMMAND \n' ...
             ' WINDOW.\n\n'] )
end

% Timer stop callback function
function func_stop(objTim,~,objDev)
    fprintf('\n ## Stop Timer ##\n')
    stop(objTim)
    delete(objTim)
    fprintf(' > Timer object has been DELETED! \n')
    release(objDev)
    fprintf(' > AudioDeviceReader object has been RELEASED! \n')
end


