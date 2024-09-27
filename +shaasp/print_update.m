function [] = print_update(str)
% PRINT_UPDATE - Update a fprintf line in the matlab workspace.
%
% Description
%
%   This function fprintf's a sting, but if the sting is the same length
%   as the last string this function printed, then the function will
%   delete the previous string and replace it with the new string. 
%
%   This function is mainly used in for loop counting.
%
% Inputs 
%
%   str     | A formatted char array to be printed. Usually created by
%             a sprintf function. 
%
% Example 
%
%   for i = 1:100
%       prntUpdate(sprintf('test %4i / %4i \n',i,100));
%       pause(0.1);
%   end
%
% Author: Lachlan Birnie
% Audio & Acoustic Signal Processing Group - Australian National University
% Email: Lachlan.Birnie@anu.edu.au
% Website: https://github.com/lachlanbirnie
% Creation: 15-August-2019
% Last revision: 24-July-2024

    persistent prevLength;
    
    if (nargin == 0)
        clear('prevLength');
        return;
    end
    
    if isempty(prevLength)
        prevLength = 0;
    end
    
    if (length(str) == prevLength)
        fprintf(repmat('\b', [1, prevLength]));
    end
    
    fprintf(str);
    prevLength = length(str);
    
end