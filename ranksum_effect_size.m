function [p,r,U,rsum,h] = ranksum_effect_size(x,y)

% Function to calculate the Mann-Whitney U statistic and the rank-biserial
% correlation based on the sample sizes and the output returned by Matlab
% function "ranksum" - namely, the sum of ranks of the first variable input
% into ranksum.
%
% If running ranksum as [p,h,stats] = ranksum(x,y), then:
% Input:
% x         1st sample
% y         2nd sample
%
% Output:
% p         P value (returned by ranksum)
% r         rank-biserial correlation
% U         Mann-Whitney U statistic
% rsum      sum of ranks of x (returned by ranksum as stats.ranksum)
% h         Test decision - h=0/1 <-> H0 not/rejected (returned by ranksum)
%
% Reference:
% Dave S. Kerby (2014) The simple difference formula: an approach to 
% teaching nonparametric correlation. Innovative Teaching 1(3).
% http://www.amsciepub.com/doi/pdf/10.2466/11.IT.3.1
%
% Frantisek Vasa, % fv247@cam.ac.uk

[p,h,stats] = ranksum(x,y);     % run ranksum
n1 = length(x);                 % size of first sample
n2 = length(y);                 % size of second sample
rsum = stats.ranksum;           % sum of ranks of x
U = rsum-(n1*(n1+1)/2);         % U statistic
r = 1-(2*U)/(n1*n2);            % rank-biserial correlation

end
