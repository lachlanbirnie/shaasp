function [Ainv] = pagepinv(A, tol)
% PAGEPINV - Pesudo inverse of 2D matrix pages in N-D matrix stack.
% NOTE: MATLAB R2024a has inbuilt pagepinv function now.
% A = (n,m, i,j,k)
% Ainv = (m,n, i,j,k)
% Where the (n,m) page is inverted for each stack (i,j,k).
%
% Syntax:  [output1,output2] = function_name(input1,input2,input3)
%
% Inputs:
%   A - N-D matrix with dimensions 1,2 being (n,m) pages.
%   tol - Tolerance.
%
% Outputs:
%   Ainv - N-D matrix with dimensions 1,2 being (m,n) inverted pages.
%   output2 - Description
%
% Example: 
%   Line 1 of example
%   Line 2 of example
%   Line 3 of example
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: OTHER_FUNCTION_NAME1,  OTHER_FUNCTION_NAME2
%
% Author: Lachlan Birnie
% Audio & Acoustic Signal Processing Group - Australian National University
% Email: Lachlan.Birnie@anu.edu.au
% Website: https://github.com/lachlanbirnie
% Creation: 11-Jul-2023
% Last revision: 11-Jul-2023

    % SVD.
    [U,S,V] = pagesvd(A);
    
    % Calulate and apply tolerance.
    if (nargin == 1) || isempty(tol)
        tol = max(size(A,[1,2])) .* eps(pagenorm(A));
    end
    S(S < tol) = 0;

    % Pseudoinverse given by Ainv = V S^-1 U^* .
    % Vh = pagectranspose(V);
    Uh = pagectranspose(U);

    oneonS = zeros(size(S));
    oneonS(S~=0) = 1./S(S~=0);
    oneonS = pagetranspose(oneonS);

    % oneonS = pageinv(S);
    
    Ainv = pagemtimes(pagemtimes(V,oneonS),Uh);
end


% function [Ainv] = pagepinv(A, tol)
% % PAGEPINV - Pesudo inverse of 2D matrix pages in N-D matrix stack.
% % A = (n,m, i,j,k)
% % Ainv = (m,n, i,j,k)
% % Where the (n,m) page is inverted for each stack (i,j,k).
% %
% % Syntax:  [output1,output2] = function_name(input1,input2,input3)
% %
% % Inputs:
% %   A - N-D matrix with dimensions 1,2 being (n,m) pages.
% %   tol - Tolerance.
% %
% % Outputs:
% %   Ainv - N-D matrix with dimensions 1,2 being (m,n) inverted pages.
% %   output2 - Description
% %
% % Example: 
% %   Line 1 of example
% %   Line 2 of example
% %   Line 3 of example
% %
% % Other m-files required: none
% % Subfunctions: none
% % MAT-files required: none
% %
% % See also: OTHER_FUNCTION_NAME1,  OTHER_FUNCTION_NAME2
% %
% % Author: Lachlan Birnie
% % Audio & Acoustic Signal Processing Group - Australian National University
% % Email: Lachlan.Birnie@anu.edu.au
% % Website: https://github.com/lachlanbirnie
% % Creation: 11-Jul-2023
% % Last revision: 11-Jul-2023
% 
%     % SVD.
%     [U,S,V] = pagesvd(A);
% 
%     % Calulate and apply tolerance.
%     if (nargin == 1) || isempty(tol)
%         tol = max(size(A,[1,2])) .* eps(pagenorm(A));
%     end
%     S(S < tol) = 0;
% 
%     % Pseudoinverse given by Ainv = V S^-1 U^* .
%     % Vh = pagectranspose(V);
%     Uh = pagectranspose(U);
% 
%     % oneonS = zeros(size(S));
%     % oneonS(S~=0) = 1./S(S~=0);
% 
%     oneonS = pageinv(S);
% 
%     Ainv = pagemtimes(pagemtimes(V,oneonS),Uh);
% end