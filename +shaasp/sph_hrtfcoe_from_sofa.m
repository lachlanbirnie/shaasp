function [hnm_left, hnm_right] = sph_hrtfcoe_from_sofa(file_hrtf,N,wlen,nfft,hp_eq)
% Get hnm HRTF coefficients using MagLS method.
% hnm = [nm, fbin]
%
% Author:   Lachlan Birnie
% Creation: 15 September 2022

    % Load sofa file (requires API_MO package).
    run('API_MO\SOFAstart.m');
    Sofa = SOFAload(file_hrtf);

    % Constants.
    c = 343;
    fs = Sofa.Data.SamplingRate;
    f_crit = 2000;  % 2kHz

    % Source positions (r, theta, phi) == radius, elevation, azimuth.
    if strcmp(Sofa.SourcePosition_Units, 'degree, degree, metre') && any(Sofa.SourcePosition(:,2) < 0)
        r = Sofa.SourcePosition(:,3);
        t = wrapToPi(deg2rad(90 - Sofa.SourcePosition(:,2)));
        p = wrapTo2Pi(deg2rad(Sofa.SourcePosition(:,1)));
    else
        error('need to add new coordinate notation');
    end

    % Frequency bands.
    fbin = linspace(0, fs, nfft+1);
    fbin = fbin(1 : end-1);

    % Window HRIRs (rectangle)
    hrir_left = permute(Sofa.Data.IR(:, 1, :), [1,3,2]);  % [source, sample]
    hrir_right = permute(Sofa.Data.IR(:, 2, :), [1,3,2]);
    
    % Apply headphone EQ.
    if (nargin >= 5)
        fprintf('- Applying headphone EQ\n');
        hrir_left_eq = zeros(size(hrir_left,1), size(hrir_left,2)+length(hp_eq)-1);
        hrir_right_eq = zeros(size(hrir_right,1), size(hrir_right,2)+length(hp_eq)-1);
        for i = (1 : size(hrir_left,1))
            hrir_left_eq(i,:) = conv(hrir_left(i,:), hp_eq);
            hrir_right_eq(i,:) = conv(hrir_right(i,:), hp_eq);
        end
    end

    if (size(hrir_left, 2) > wlen)
        hrir_left = hrir_left(:, 1:wlen);
        hrir_right = hrir_right(:, 1:wlen);
    end

    % Get HRTFs.
    hrtf_left = fft(hrir_left, nfft, 2);  % [fbin, source]
    hrtf_right = fft(hrir_right, nfft, 2);

    % Build Ynm and invYnm matrix.
    ynm = shaasp.sph_ynm(N, t, p);  %  [source, nm]

    inv_ynm = (ynm.' * ynm)^(-1) * ynm.';
%     inv_ynm = pinv(ynm);

    % Get coefficients for each frequency bin.
    hnm_left_magls = zeros((N+1)^2, length(fbin));
    hnm_right_magls = zeros((N+1)^2, length(fbin));
    
    %hnm_left_linear = zeros((N+1)^2, length(fbin));  % linear option.
    %hnm_right_linear = zeros((N+1)^2, length(fbin));

    for ind = 1:length(fbin)

        %hnm_left_linear(:,ind) = inv_ynm * hrtf_left(:,ind);  % linear.
        %hnm_right_linear(:,ind) = inv_ynm * hrtf_right(:,ind);

        if fbin(ind) < f_crit

            hnm_left_magls(:,ind) = inv_ynm * hrtf_left(:,ind);
            hnm_right_magls(:,ind) = inv_ynm * hrtf_right(:,ind);

        else

            phi_left = angle(ynm * hnm_left_magls(:,ind-1));
            phi_right = angle(ynm * hnm_right_magls(:,ind-1));

            hnm_left_magls(:,ind) = inv_ynm * (abs(hrtf_left(:,ind)) .* exp(1i .* phi_left));
            hnm_right_magls(:,ind) = inv_ynm * (abs(hrtf_right(:,ind)) .* exp(1i .* phi_right));

        end 

    end

    hnm_left = hnm_left_magls(:, 1:nfft/2+1);
    hnm_right = hnm_right_magls(:, 1:nfft/2+1);
    
end