function [recording, Rec] = playrecord(DeviceName,play_signal,fs,varargin)
% [recording, Rec] = playrecord(DeviceName,play_signal,fs)
% -------------------------------------------------------------------------
% File      : playrecord.m
% Creation  : 24-11-20 
% Author    : Lachlan Birnie
% Version   : 22-03-23
%
% Description
%
%   Generic function for simultaneous play and record. 
%   
%   This function uses the audioPlayerRecorder System object 
%   from the MATLAB Audio Toolbox.
%
% Inputs
%
%   DeviceName  | 'string' name of audioPlayerRecorder device.
%   play_signal  | [samples by channel] signal to play.
%   fs  | [Hz] sampling frequency.
%   
%   - optional inputs - (defaults) 
%   'FrameSize', [N] frame size, number of samples (2048)
%   'BufferSize', [N] buffer size, numbe of samples (2048)
%   'PlayerVolume', [N] playback volume, scalar 0-1 (1)
%   'PlayerChannelMapping', [1 by channel numbers] output channels ([1])
%   'PlayerTotalChannels', [N] total number of playback channels (1)
%   'RecorderChannelMapping', [1 by channel numbers] recording channels ([1])
%   'RecorderTotalChannels', [N] total number of playback channels (1)
%   'BitDepth'
%   'recordOnlyDur', [seconds] {recording only, no playback} duration (0)
%   'zeroPadStartDur', [seconds] zero-pad start of playback (0)
%   'zeroPadEndDur', [seconds] zero-pad end of playback (0)
%   'prntStatus', [bool] 0 = dont fprintf 'recording' 'done' (1)
%
% Outputs
%
%   recording   | [samples by channels] recording result.
%   Rec         | structure with recording information and data:
%     Rec.DeviceName = DeviceName input
%     Rec.FrameSize = frame size used. 
%     Rec.BufferSize = buffer size used.
%     Rec.SampleRate = sampling frequency used.
%     Rec.PlayerVolume = playback volume used. 
%     Rec.PlayerChannelMapping = channel playback was on.
%     Rec.PlayerTotalChannels = total number of playback channels (used or not)
%     Rec.RecorderChannelMapping = channel recording is from.
%     Rec.RecorderTotalChannels = total num recording channels (used or not)
%     Rec.BitDepth = ('24-bit integer')
%     Rec.recordOnlyDur = {recording only, no playback} recording duration
%     Rec.zeroPadStartDur = zero-padding time at start of playback.
%     Rec.zeroPadEndDur = zero-padding time at end of playback.
%     Rec.play_signal = signal that was played
%     Rec.recSig = raw signal that was recorded.
%     Rec.APR = audioPlayerRecorder handle.
%     Rec.signal = processed signal that was returned (removed zero-paddeding from start of recording)


    % Default inputs.
    Rec = struct();
    Rec.DeviceName = DeviceName;
    Rec.FrameSize = 2048;
    Rec.BufferSize = 2048;
    Rec.SampleRate = fs;
    Rec.PlayerVolume = 1;
    Rec.PlayerChannelMapping = 1;
    Rec.PlayerTotalChannels = 1;
    Rec.RecorderChannelMapping = 1;
    Rec.RecorderTotalChannels = 1;
    Rec.BitDepth = '24-bit integer';
    Rec.recordOnlyDur = 0;
    Rec.zeroPadStartDur = 0;
    Rec.zeroPadEndDur = 0;
    prntStatus = true;   % print words. 
    
    
    % Optional Inputs.
    if (nargin > 3)
        for i = (1:nargin-3)
            if ischar(varargin{i})
                switch varargin{i}
                    case 'FrameSize', Rec.FrameSize = varargin{i+1};
                    case 'BufferSize', Rec.BufferSize = varargin{i+1};
                    case 'PlayerVolume', Rec.PlayerVolume = varargin{i+1};
                    case 'PlayerChannelMapping', Rec.PlayerChannelMapping = varargin{i+1};
                    case 'PlayerTotalChannels', Rec.PlayerTotalChannels = varargin{i+1};
                    case 'RecorderChannelMapping', Rec.RecorderChannelMapping = varargin{i+1};
                    case 'RecorderTotalChannels', Rec.RecorderTotalChannels = varargin{i+1};
                    case 'BitDepth', Rec.BitDepth = varargin{i+1};
                    case 'recordOnlyDur', Rec.recordOnlyDur = varargin{i+1};
                    case 'zeroPadStartDur', Rec.zeroPadStartDur = varargin{i+1};
                    case 'zeroPadEndDur', Rec.zeroPadEndDur = varargin{i+1};
                    case 'prntStatus', prntStatus = varargin{i+1};
                    otherwise
                        error('Unrecognised optional input name');
                end
            end
        end
    end
    
    
    % Check channel mapping.
    if (Rec.PlayerTotalChannels < numel(Rec.PlayerChannelMapping))
        fprintf(['Rec: WARNING! PlayerTotalChannels must be >= PlayerChannelMapping \n'...
                 'Rec: Changing PlayerTotalChannels to numel(PlayerChannelMapping) \n']);
        Rec.PlayerTotalChannels = numel(Rec.PlayerChannelMapping);
    end
    
    if (Rec.RecorderTotalChannels < numel(Rec.RecorderChannelMapping))
        fprintf(['Rec: WARNING! RecorderTotalChannels must be >= RecorderChannelMapping \n'...
                 'Rec: Changing RecorderTotalChannels to numel(RecorderChannelMapping) \n']); 
        Rec.RecorderTotalChannels = numel(Rec.RecorderChannelMapping);
    end
    
    
    % Map the play signal to the player channels.
    if isempty(play_signal)
        % Only take a recording, playback zeros.
        mapped_signal = zeros(Rec.recordOnlyDur * Rec.SampleRate, Rec.PlayerTotalChannels);
        if prntStatus
            fprintf(['Rec: Record only mode, because play_signal is empty. \n' ...
                    ,'Rec: Duration of recording given by "recordOnlyDur" = %i sec \n'], Rec.recordOnlyDur);
        end
                  
    elseif (size(play_signal,2) == 1)
        % Play a single channel signal over all player channels.
        mapped_signal = repmat(play_signal, [1, Rec.PlayerTotalChannels]);
                  
    else
        % Play a multi-channel signal over equal number of channels.
        if (size(play_signal,2) ~= Rec.PlayerTotalChannels)
            error('Rec: The number of play_signal channels must match the playerChannelMapping.');
        end
        mapped_signal = play_signal;

    end
    Rec.play_signal = mapped_signal;
    
    
    % Zero pad start and end of play signal.
    Rec.play_signal = [zeros(Rec.zeroPadStartDur * Rec.SampleRate, Rec.PlayerTotalChannels)...
                      ;Rec.play_signal ...
                      ;zeros(Rec.zeroPadEndDur * Rec.SampleRate, Rec.PlayerTotalChannels)...
                      ];
    
    
    % Pad play_signal for equal number of frames.
    original_length = size(Rec.play_signal,1);
    total_frames = ceil(size(Rec.play_signal,1) / Rec.FrameSize);
    total_samples = total_frames * Rec.FrameSize;
    Rec.play_signal = [Rec.play_signal; zeros(total_samples - size(Rec.play_signal,1), Rec.PlayerTotalChannels)];
              
             
    % Adjust play volume.
    Rec.play_signal = Rec.play_signal .* Rec.PlayerVolume;
    
    
    % Make audioPlayerRecorder.
    Rec.APR = audioPlayerRecorder();
    Rec.APR.Device = Rec.DeviceName;
    Rec.APR.SampleRate = Rec.SampleRate;
    Rec.APR.BufferSize = Rec.BufferSize;
    Rec.APR.BitDepth = '24-bit integer';
    Rec.APR.PlayerChannelMapping = Rec.PlayerChannelMapping;
    Rec.APR.RecorderChannelMapping = Rec.RecorderChannelMapping;
    
    Rec.APRinfo = info(Rec.APR);
    if (Rec.PlayerTotalChannels > Rec.APRinfo.MaximumPlayerChannels)
        fprintf('Rec: WARNING! PlayerTotalChannels is > supported maximum');
    end
    if (Rec.RecorderTotalChannels > Rec.APRinfo.MaximumRecorderChannels)
        fprintf('Rec: WARNING! RecorderTotalChannels is > supported maximum');
    end
    
    
    % Pre-allocate play/record variables.
    frame_index = reshape((1 : total_samples), [Rec.FrameSize, total_frames]);
    underRun = uint32(0);
    overRun = uint32(0);
    
    
    % Allocate recording signal.
    Rec.record_signal_raw = zeros(size(Rec.play_signal,1), Rec.RecorderTotalChannels);
    
    
    % Setup audioPlayerRecorder device and send dummy frame.
    if prntStatus, fprintf('Rec: Recording Start!\n'); end
    setup(Rec.APR, zeros(Rec.FrameSize, Rec.PlayerTotalChannels));
    [Rec.record_signal_raw(frame_index(:,1),:), ~, ~] = Rec.APR(zeros(Rec.FrameSize, Rec.PlayerTotalChannels));
    
    
    % - Play & Record - 
    for i = (1 : total_frames)
        [Rec.record_signal_raw(frame_index(:,i),:), underRun, overRun] = Rec.APR(Rec.play_signal(frame_index(:,i), :));
        if (underRun || overRun)
            warning(['Rec: missed samples in playRecordSweep.m\n' ...
                     'underRun = %i\n' ...
                     'overRun  = %i\n' ...  
                     'i_frame  = %i'] ...
                     ,underRun, overRun, i);
        end
    end
    release(Rec.APR);
    if prntStatus, fprintf('Rec: Recording Finished!\n'); end
    
    
    % Post processing.
    Rec.record_signal = Rec.record_signal_raw;  % Copy raw signal.
    Rec.record_signal = Rec.record_signal(1:original_length, :);  % Remove equal number of frames padding.
    Rec.record_signal = Rec.record_signal(Rec.zeroPadStartDur*Rec.SampleRate+1:end, :);  % Remove pre-padded zeros.
    
    
    % Outputs.
    recording = Rec.record_signal;

end  % end [recording, Rec] = playrecord(DeviceName,play_signal,fs,varargin)