function [r,t,p,x,y,z,w] = sampling_positions_ANU_AASP_ROSS_loudspeaker_array()
% SAMPLING_POSITIONS_ANU_AASP_BAB_LOUDSPEAKER_ARRAY - (x,y,z) (r,theta,phi)
%
% Description
%
%   Positions of the 30 channel loudspeaker array in ANU AASP lab.
%   (The old one in the BAB Lab).
%
% Outputs
%
%   r       | [Q by 1] radius of each speaker
%   t       | [Q by 1] elevation (0>pi) of each speaker
%   p       | [Q by 1] rotational position (0>2pi) of each speaker
%   x       | [Q by 1] x-position of each speaker w.r.t center
%   y       | [Q by 1] y-position ""
%   z       | [Q by 1] z-position ""
%   w       | [Q by 1] sampling weight of each speaker (unknown)
%
% Author: Lachlan Birnie
% Audio & Acoustic Signal Processing Group - Australian National University
% Email: Lachlan.Birnie@anu.edu.au
% Website: https://github.com/lachlanbirnie
% Creation: 16-Dec-2022
% Last revision: 16-Dec-2024

t = [124
    124
    124
    90
    90
    69
    21
    56
    90
    90
    69
    21
    21
    56
    56
    21
    159
    159
    111
    124
    90
    56
    56
    90
    159
    159
    111
    90
    56
    124
    124
    90] .* (2*pi/360);

p = [247.5
    202.5
    157.5
    180
    135
    112.5
    112.5
    157.5
    225
    270
    292.5
    292.5
    22.5
    247.5
    202.5
    202.5
    202.5
    292.5
    292.5
    337.5
    315
    337.5
    22.5
    0
    22.5
    112.5
    112.5
    90
    67.5
    67.5
    22.5
    45] .* (2*pi/360);

r = ones([32, 1]) .* 1.11;

% Cartesian (x,y,z) coordinates (m).
x = (r .* sin(t) .* cos(p));
y = (r .* sin(t) .* sin(p));
z = (r .* cos(t));

w = [];  % Unknown.

% Plot the positions if no output.
if ~nargout
    shaasp.sampling_positions_plot_on_sphere(x,y,z);
end

end