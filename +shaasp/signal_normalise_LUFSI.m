function [sig] = signal_normalise_LUFSI(sig, fs, target_LUFSI)
% SIGNAL_NORMALIZE_LUFSI - Normalise a signal to the target LUFSI
%
% Syntax:  [sig] = signal_normalise_LUFSI(sig, fs, target_LUFSI)
%
% Inputs:
%   sig - [samp, chn] signal to normalise
%   fs - Scalar (Hz), sampling frequency
%   target_LUFSI - Scalar (dB), normalise signals to this power
%
% Outputs:
%   sig - [samp, chn] normalised signal
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none
%
% Author: Lachlan Birnie
% Audio & Acoustic Signal Processing Group - Australian National University
% Email: Lachlan.Birnie@anu.edu.au
% Website: https://github.com/lachlanbirnie
% Creation: 03-July-2023
% Last revision: 03-July-2023

num_channels = size(sig, 2);

for ch = (1 : num_channels)
    loudness = integratedLoudness(sig(:, ch), fs);
    gaindB = target_LUFSI - loudness;
    gain = 10.^(gaindB ./ 20);
    sig(:, ch) = sig(:, ch) .* gain;
    new_loudness = integratedLoudness(sig(:, ch), fs);
    fprintf('Normalised loudness: %+2.2f -> %+2.2f LUFSI. \n', loudness, new_loudness);
    % Check for clipping +-1.
    peak = max(abs(sig(:)));
    if peak > 1
        warning('Signal is clipping after gain. \n');
    end
end

end