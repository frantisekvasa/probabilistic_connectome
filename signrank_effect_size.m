function [p,r,rsum_pos,rsum_neg,h] = signrank_effect_size(x,y)

% Function to calculate the rank-biserial
% correlation based on the sample sizes and the output returned by Matlab
% function "signrank" - namely, the sum of ranks of the first variable input
% into signrank.
%
% If running signrank as [p,h,stats] = signrank(x,y), then:
% Input:
% x         1st sample
% y         2nd sample
%
% Output:
% p         P value (returned by signrank)
% r         rank-biserial correlation
% rsum      sum of ranks of x (returned by signrank as stats.signedrank)
% h         Test decision - h=0/1 <-> H0 not/rejected (returned by signrank)
%
% Reference:
% Dave S. Kerby (2014) The simple difference formula: an approach to 
% teaching nonparametric correlation. Innovative Teaching 1(3).
% http://www.amsciepub.com/doi/pdf/10.2466/11.IT.3.1
%
% Frantisek Vasa, % fv247@cam.ac.uk

[p,h,stats] = signrank(x,y);    % run signrank
n1 = length(x);                 % size of first sample
rsum = sum(1:n1);               % total sum of ranks (see Kerby 2014, p.5 bottom right)

%%%

abs_rank = tiedrank(abs(x-y));
r_pos = sum(abs_rank((x-y)>0));
r_neg = sum(abs_rank((x-y)<0));

%%%

%n2 = length(y);                % size of second sample == N1, OTHERWISE SIGNRANK WON'T RUN
rsum_pos = stats.signedrank;    % sum of the ranks for the points with positive sign
rsum_neg = rsum-rsum_pos;       % sum of the ranks for the points with negative sign
r = rsum_pos/rsum - rsum_neg/rsum;            % rank-biserial correlation


end
