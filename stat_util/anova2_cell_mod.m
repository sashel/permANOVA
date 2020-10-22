%% Modified version of the ANOVA2_CELL function of Arnaud Delorme, SCCN/INC/UCSD, La Jolla, 2010
% Source: http://www.mathworks.com/matlabcentral/fileexchange/27960-resampling-statistical-toolkit/content/statistics/anova2_cell.m
% 
% Header A. Delorme 
% anova2_cell() - compute F-values in cell array using ANOVA.
%
% Usage:
%    >> [FC FR FI dfc dfr dfi] = anova2_cell( data );
%
% Inputs:
%   data       = data consisting of PAIRED arrays to be compared. The last
%                dimension of the data array is used to compute ANOVA.
% Outputs:
%   FC   - F-value for columns.
%   FR   - F-value for rows.
%   FI   - F-value for interaction.
%   dfc  - degree of freedom for columns.
%   dfr  - degree of freedom for rows.
%   dfi  - degree of freedom for interaction.
%
% Note: the advantage over the ANOVA2 function of Matlab statistical
%       toolbox is that this function works on arrays (see examples). Note
%       also that you still need the statistical toolbox to assess
%       significance using the fcdf() function. The other advantage is that
%       this function will work with complex numbers.
%
% Example:
%   a = { rand(1,10) rand(1,10) rand(1,10); rand(1,10) rand(1,10) rand(1,10) }
%   [FC FR FI dfc dfr dfi] = anova2_cell(a)
%   signifC = 1-fcdf(FC, dfc(1), dfc(2))
%   signifR = 1-fcdf(FR, dfr(1), dfr(2))
%   signifI = 1-fcdf(FI, dfi(1), dfi(2))
%
%   % for comparison
%   anova2(  [ a{1,1}' a{1,2}' a{1,3}'; a{2,1}' a{2,2}' a{2,3}' ], 10)
%
%   b = { [ a{1,1}; a{1,1} ] [ a{1,2}; a{1,2} ] [ a{1,3}; a{1,3} ];
%         [ a{2,1}; a{2,1} ] [ a{2,2}; a{2,2} ] [ a{2,3}; a{2,3} ] }
%   [FC FR FI dfc dfr dfi] = anova2_cell(b)
%
%   c{1,1} = reshape(repmat(b{1,1}, [2 1]),2,2,10);
%   c{1,2} = reshape(repmat(b{1,2}, [2 1]),2,2,10);
%   c{1,3} = reshape(repmat(b{1,3}, [2 1]),2,2,10);
%   c{2,3} = reshape(repmat(b{2,3}, [2 1]),2,2,10);
%   c{2,2} = reshape(repmat(b{2,2}, [2 1]),2,2,10);
%   c{2,1} = reshape(repmat(b{2,1}, [2 1]),2,2,10)
%   [FC FR FI dfc dfr dfi] = anova2_cell(c)
%
% Author: Arnaud Delorme, SCCN/INC/UCSD, La Jolla, 2005
%
% Reference:
%   Schaum's outlines in statistics (3rd edition). 1999. Mc Graw-Hill.

% Copyright (C) Arnaud Delorme
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

function [FA, FB, FI, freeA, freeB, freeI] = anova2_cell_mod(data)

% compute all means and all std
% -----------------------------
a = size(data,1);
b = size(data,2);
c = size(data{1}, 2);



VE = zeros( size(data{1},1),1);
m  = zeros( [ size(data{1},1) size(data) ] );
for i = 1:a
    for ii = 1:b
        tmpm = mean(data{i,ii}, 2);
        m(:,i,ii) = tmpm;
        VE        = VE+sum( (data{i,ii}-repmat(tmpm, [1 size(data{i,ii},2)])).^2, 2);
    end;
end;
X  = mean(mean(m,3),2);
Xj = mean(m,3);
Xk = mean(m,2);
VR = b*c*sum( (Xj-repmat(X, [1 size(Xj,2)])).^2, 2 );
VC = a*c*sum( (Xk-repmat(X, [1 1 size(Xk,3)])).^2, 3 );

Xj = repmat(Xj, [1 1 size(m,3) ]);
Xk = repmat(Xk, [1 size(m,2)  1]);
VI = c*sum( sum( ( m - Xj - Xk + repmat(X, [1 size(m,2) size(m,3)]) ).^2, 3), 2 );
SR2 = VR/(a-1);
SC2 = VC/(b-1);
SI2 = VI/(a-1)/(b-1);
SE2 = VE/(a*b*(c-1));

FA = SR2./SE2; % rows
FB = SC2./SE2; % columns
FI = SI2./SE2; % interaction

freeA = [ a-1 a*b*(c-1) ];
freeB = [ b-1 a*b*(c-1) ];
freeI = [ (a-1)*(b-1) a*b*(c-1) ];
