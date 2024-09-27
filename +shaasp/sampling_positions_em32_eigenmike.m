function [r,t,p,x,y,z,w] = sampling_positions_em32_eigenmike()
% SAMPLING_POSITIONS_EM32_EIGENMIKE - Eigenmike em32 sensor positions.
%
% Description
%
%   Properties of the EM32 EigenMike spherical microphone array.
%
%   Positions obtained from:
%   https://mhacoustics.com/sites/default/files/ReleaseNotes.pdf
%
% Inputs
%
%   N       | Truncation order of inversion matrix
%   k       | Wave numbers (frequencies) of inversions matrices
%   open    | 1 = open spherical array, 0 = [] = rigid spherical array.
%
% Outputs
%
%   r       | [Q by 1] radius of each sensor
%   t       | [Q by 1] elevation (0>pi) of each sensor
%   p       | [Q by 1] rotational position (0>2pi) of each sensor
%   x       | [Q by 1] x-position of each sensor w.r.t center
%   y       | [Q by 1] y-position ""
%   z       | [Q by 1] z-position ""
%   w       | [Q by 1] sampling weight of each sensor 
%
% Author: Lachlan Birnie
% Audio & Acoustic Signal Processing Group - Australian National University
% Email: Lachlan.Birnie@anu.edu.au
% Website: https://github.com/lachlanbirnie
% Creation: 25-March-2019
% Last revision: 02-March-2023

    % Spherical (r, theta, phi) coordinates (rads).
    r = ones(32, 1) .* 0.042;
    t = [   69  ;
            90  ;
            111 ;
            90  ;
            32  ;
            55  ;
            90  ;
            125 ;
            148 ;
            125 ;
            90  ;
            55  ;
            21  ;
            58  ;
            121 ;
            159 ;
            69  ;
            90  ;
            111 ;
            90  ;
            32  ;
            55  ;
            90  ;
            125 ;
            148 ;
            125 ;
            90  ;
            55  ;
            21  ;
            58  ;
            122 ;
            159 ;   ] .* (2*pi/360);
    p = [   0   ;
            32  ;
            0   ;
            328 ;
            0   ;
            45  ;
            69  ;
            45  ;
            0   ;
            315 ;
            291 ;
            315 ;
            91  ;
            90  ;
            90  ;
            89  ;
            180 ;
            212 ;
            180 ;
            148 ;
            180 ;
            225 ;
            249 ;
            225 ;
            180 ;
            135 ;
            111 ;
            135 ;
            269 ;
            270 ;
            270 ;
            271 ;   ] .* (2*pi/360);        
        
    % Cartesian (x,y,z) coordinates (m).
    x = (r .* sin(t) .* cos(p));
    y = (r .* sin(t) .* sin(p));
    z = (r .* cos(t));
    
    % Uniform sampling weights.
    w = ones(size(t));

end