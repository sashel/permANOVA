%% Adapted version of the ANOVA2RM_CELLfunction of Arnaud Delorme, SCCN/INC/UCSD, La Jolla, 2010
% Source: http://www.mathworks.com/matlabcentral/fileexchange/27960-resampling-statistical-toolkit/content/statistics/anova2rm_cell.m
% Adapted by Saskia Helbling
%% Header A. Delorme
%  anova2rm_cell() - compute F-values in cell array using repeated measurement ANOVA.

% Usage:
%    >> [FC FR FI dfc dfr dfi] = anova2rm_cell( data );
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
% Note: this function is inspired from rm_anova available at 
%       http://www.mathworks.se/matlabcentral/fileexchange/6874-two-way-rep
%       eated-measures-anova
%       It allows for fast computation of about 20 thousands ANOVA per
%       second. It is different from anova2_cell which mimics the ANOVA
%       fonction from the Matlab statistical toolbox. This function
%       computes true repeated measure ANOVA.
%
% Example:
%   a = { rand(1,10) rand(1,10) rand(1,10); rand(1,10) rand(1,10) rand(1,10) }
%   [FC FR FI dfc dfr dfi] = anova2rm_cell(a)
%   signifC = 1-fcdf(FC, dfc(1), dfc(2))
%   signifR = 1-fcdf(FR, dfr(1), dfr(2))
%   signifI = 1-fcdf(FI, dfi(1), dfi(2))
%
%   % for comparison 
%   z = zeros(10,1); o = ones(10,1); t = ones(10,1)*2;
%   rm_anova2(  [ a{1,1}';a{1,2}';a{1,3}';a{2,1}';a{2,2}';a{2,3}' ], ...
%               repmat([1:10]', [6 1]), [o;o;o;z;z;z], [z;o;t;z;o;t], {'a','b'})
%
%   c = { rand(200,400,10) rand(200,400,10); ...
%         rand(200,400,10) rand(200,400,10)};
%   [FC FR FI dfc dfr dfi] = anova2rm_cell(c) % computes 200x400 ANOVAs
%
% Author: Arnaud Delorme, SCCN/INC/UCSD, La Jolla, 2010

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

function [fA, fB, fAB, dfApair, dfBpair, dfABpair] = anova2rm_cell_mod(data)

a = size(data,1);
b = size(data,2);
n  = size(data{1},2);

    AB = zeros(size(data{1},1),a,b);
    AS = zeros(size(data{1},1),a,n);
    BS = zeros(size(data{1},1),b,n);
    sq = zeros(size(data{1},1),1);
    for ind1 = 1:a
        for ind2 = 1:b
            AB(:,ind1,ind2) = sum(data{ind1,ind2},2);
            AS(:,ind1,:)    = AS(:,ind1,:) + reshape(data{ind1,ind2},size(data{1},1),1,n);
            BS(:,ind2,:)    = BS(:,ind2,:) + reshape(data{ind1,ind2},size(data{1},1),1,n);
            sq              = sq + sum(data{ind1,ind2}.^2,2);
        end
    end

A = sum(AB,3); % sum across columns, so result is ax1 column vector
B = sum(AB,2); % sum across rows, so result is 1xb row vector
S = sum(AS,2); % sum across columns, so result is 1xs row vector
T = sum(sum(A,2),3); % could sum either A or B or S, choice is arbitrary

% degrees of freedom
dfA = a-1;
dfB = b-1;
dfAB = (a-1)*(b-1);
dfAS = (a-1)*(n-1);
dfBS = (b-1)*(n-1);
dfABS = (a-1)*(b-1)*(n-1);

% bracket terms (expected value)
expA  = sum(A.^2,2)./(b*n);
expB  = sum(B.^2,3)./(a*n);
expAB = sum(sum(AB.^2,3),2)./n;
expS  = sum(S.^2,3)./(a*b);
expAS = sum(sum(AS.^2,2),3)./b;
expBS = sum(sum(BS.^2,2),3)./a;
expY  = sq; %sum(Y.^2);
expT  = T.^2 / (a*b*n);

% sums of squares
ssA   = expA - expT;
ssB   = expB - expT;
ssAB  = expAB - expA - expB + expT;
ssAS  = expAS - expA - expS + expT;
ssBS  = expBS - expB - expS + expT;
ssABS = expY - expAB - expAS - expBS + expA + expB + expS - expT;

% mean squares
msA   = ssA / dfA;
msB   = ssB / dfB;
msAB  = ssAB / dfAB;
msAS  = ssAS / dfAS;
msBS  = ssBS / dfBS;
msABS = ssABS / dfABS;

% f statistic
fA = msA ./ msAS;
fB = msB ./ msBS;
fAB = msAB ./ msABS; 
dfApair  = [dfA dfAS];
dfBpair  = [dfB dfBS];
dfABpair = [dfAB dfABS];
