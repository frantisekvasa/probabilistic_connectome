r% Code used to perform primary analyses described in the manuscript: 
%
% "Probabilistic thresholding of functional connectomes: application to schizophrenia"
% by Frantisek Vasa, Edward T. Bullmore and Ameera X. Patel
%
% Frantisek Vasa, October 2015 - November 2017
% fv247@cam.ac.uk
%
% Dependencies (available at https://github.com/frantisekvasa/structural_network_development):
%
% ranksum_effect_size       calculates ranksum effect size using the "simple difference formula" (Kerby, Compr. Psychol. 2014)
% signrank_effect_size      calculates signrank effect size using the "simple difference formula" (Kerby, Compr. Psychol. 2014)
% hist_box_plot             box plot combined with histograms
% hor_box_plot              horizontal box plot
%
% + variations on the boundedline_plot function (based around "boundedline.m", referenced below):
% boundedline_plot          two lines with confidence intervals - for topological differences across edge densities
% boundedline_plot_2p       two lines with confidence intervals, with two p-value markers for each density
% boundedline_plot_2p_drop  two lines with confidence intervals, with two-pvalue markers and drop-off in number of participants
% boundedline_plot_one      single line with confidence intervals - for paired group difference across edge densities
%
% Additional dependencies: 
%
% Brain Connectivity Toolbox    available from: https://sites.google.com/site/bctnet/
% boundedline                   function to plot lines with confindence interval bounds, available from Matlab Central: https://uk.mathworks.com/matlabcentral/fileexchange/27485-boundedline-m
% rgb                           color function, available from Matlab Central: https://uk.mathworks.com/matlabcentral/fileexchange/1805-rgb-m

%% set-up

base_path = pwd;                                            % path containing data
plot_path = '~/Desktop/prob_thr_fig/'; mkdir([plot_path])   % path where plots will be exported

% Wavelet scale selection
%wscale = 1;    % scale 1 = 0.125-0.260 Hz
wscale = 2;     % scale 2 = 0.060-0.125 Hz
%wscale = 3;    % scale 3 = 0.030-0.060 Hz
%wscale = 4;    % scale 4 = 0.015-0.030 Hz

sc = num2str(wscale);   % wavelet scale string
alpha = 0.01;           % significance

load([base_path 'COBRE_scale' sc '_toShare.mat']); % load data
% ctrl  healthy controls
% schz  patients with schizophrenia
%
% r     correlation matrices
% fdr   fdr-adjusted p-value matrices
% edof  (nodal) effective degrees of freedom

nroi = size(ctrl.r,1); % number of regions

nc = size(ctrl.r,3); % number of controls
np = size(schz.r,3); % number of patients   
nd = nc-np;          % difference (for certain plots)

triu_ind = find(triu(ones(nroi),1)); % upper triangular edge mask

% x-ranges for distributions:
deg_x = 0:420;         % degree
ndeg_x = 0:5:3000;     % normalised degree
str_x = 0:300;         % strength
dist_x = 0:0.5:160;    % distance
cmp_x = 0:500;         % connected components
wei_x = -0.6:0.005:1;  % full (upper triangular) weight

% colors (for plotting) 
ctrl.col = rgb('navy');
ctrl.col_light = rgb('slateblue');
ctrl.col_nonsig = rgb('purple');
ctrl.col_r = rgb('cyan');

schz.col = rgb('darkorange');
schz.col_light = rgb('orange');
schz.col_nonsig = rgb('red');
schz.col_r = rgb('gold');

% regional coordinates (as x,y,z)
load('COBRE_coords_420.mat'); dist = zeros(nroi,nroi);
for i = 1:1:nroi; for j = 1:1:nroi; dist(i,j) = sqrt(sum((coords(i,:)-coords(j,:)).^2)); end; end

% plotting parameters (some of which are re-defined in the script)
fsizet = 24;
fsize = 22;
fsizelab = 18;
fsize2 = 14;
fsizeh = 12;
lwd = 3;
lwd_ind = 1;
col_ind = 0.5;
msize = 30;
cex = 15;
alph = 0.3; % transparency for bounded plot
corr_type = 'Spearman'; % correlation type

%% threshold by p-value (independent of connectedness, or density rho)

%%% THRESHOLD MATRICES (by P, using "alpha" set above)

% controls
neg_c = [];     % count subjects with negative edges
for i = 1:1:nc
    ctrl.mat(:,:,i) = ctrl.r(:,:,i).*(ctrl.fdr(:,:,i)<alpha);
    if min(min(ctrl.mat(:,:,i))) < 0; neg_c = [neg_c i]; end
end
ctrl.mat_b = double(ctrl.mat>0);

% patients
neg_p = [];     % count subjects with negative edges
for i = 1:1:np
    schz.mat(:,:,i) = schz.r(:,:,i).*(schz.fdr(:,:,i)<alpha);
    if min(min(schz.mat(:,:,i))) < 0; neg_p = [neg_p i]; end
end
schz.mat_b = double(schz.mat>0);

%%% MEASURES

% controls
ctrl.cmp_sizes = [];
for i = 1:1:nc
    
    if mod(i,10) == 0; disp(['ctrl ' num2str(i)]); end
    
    % weight distributions + average weights
    full_mat = ctrl.r(:,:,i);
    [ctrl.full_wei(i,:),xi] = ksdensity(full_mat(triu_ind),wei_x);      % weight distribution
    ctrl.avg_wei(i) = mean(full_mat(triu_ind));     % average weight (of upper triangular part)
    ctrl.avg_wei_nod(i,:) = mean(full_mat,2);       % average nodal weight
    ctrl.std(i) = std(full_mat(triu_ind));          % standard deviation
    ctrl.skew(i) = skewness(full_mat(triu_ind));    % skewness
    ctrl.kurt(i) = kurtosis(full_mat(triu_ind));    % kurtosis
    
    % components
    [cmp,cmp_sizes] = get_components(ctrl.mat_b(:,:,i));
    ctrl.cmp_n(i) = length(unique(cmp));            % number of (dis)connected components
    ctrl.cmp_max(i) = max(cmp_sizes);               % size of the largest component
    ctrl.cmp_mean(i) = mean(cmp_sizes);             % average component size
    ctrl.single(i) = sum(cmp_sizes == 1);           % number of single nodes
    ctrl.cmp_score(i,:) = cmp_sizes(cmp)/max(cmp_sizes);    % how consistently is a node located in the largest component
    
    ctrl.nedge(i) = sum(sum(triu(ctrl.mat_b(:,:,i),1)));    % number of edges
    
    % degree
    ctrl.deg(i,:) = sum(ctrl.mat_b(:,:,i));
    [ctrl.deg_dens(i,:),xi] = ksdensity(ctrl.deg(i,:),deg_x,'function','pdf');
    ctrl.ndeg(i,:) = ctrl.deg(i,:)./(ctrl.nedge(i)/(nroi*(nroi-1)/2)); % density-normalised degree
    [ctrl.ndeg_dens(i,:),xi] = ksdensity(ctrl.ndeg(i,:),ndeg_x,'function','pdf');
    
    % distance
    temp = dist.*(ctrl.fdr(:,:,i)<alpha); temp = temp(temp>0);
    ctrl.dist_mean(i) = mean(temp);                 % average (global) distance
    [ctrl.dist(i,:),xi] = ksdensity(temp,dist_x);   % distance distribution
    dist_temp = nan(nroi);
    dist_temp(temp>0) = temp(temp>0);
    ctrl.dist_nod(i,:) = nanmean(dist_temp,2);      % nodal distance
    
end
clear temp_mat full_mat temp temp_b dist_temp

ctrl.dens = ctrl.nedge./(nroi*(nroi-1)/2); % edge density

% patients
schz.cmp_sizes = [];
for i = 1:1:np
    
    if mod(i,10) == 0; disp(['schz ' num2str(i)]); end
    
    % weight distributions + average weights
    
    % unthresholded
    full_mat = schz.r(:,:,i);
    [schz.full_wei(i,:),xi] = ksdensity(full_mat(triu_ind),wei_x);      % weight distribution
    schz.avg_wei(i) = mean(full_mat(triu_ind));     % average weight (of upper triangular part)
    schz.avg_wei_nod(i,:) = mean(full_mat,2);       % average nodal weight
    schz.std(i) = std(full_mat(triu_ind));          % standard deviation
    schz.skew(i) = skewness(full_mat(triu_ind));    % skewness
    schz.kurt(i) = kurtosis(full_mat(triu_ind));    % kurtosis
    
    % components
    [cmp,cmp_sizes] = get_components(schz.mat_b(:,:,i));
    schz.cmp_n(i) = length(unique(cmp));            % number of (dis)connected components
    schz.cmp_max(i) = max(cmp_sizes);               % size of the largest component
    schz.cmp_mean(i) = mean(cmp_sizes);             % average component size
    schz.single(i) = sum(cmp_sizes == 1);           % number of single nodes
    schz.cmp_score(i,:) = cmp_sizes(cmp)/max(cmp_sizes);    % how consistently is a node located in the largest component
    
    schz.nedge(i) = sum(sum(triu(schz.mat_b(:,:,i),1)));    % number of edges
    
    % degree
    schz.deg(i,:) = sum(schz.mat_b(:,:,i));
    [schz.deg_dens(i,:),xi] = ksdensity(schz.deg(i,:),deg_x,'function','pdf');
    schz.ndeg(i,:) = schz.deg(i,:)./(schz.nedge(i)/(nroi*(nroi-1)/2)); % density-normalised degree
    [schz.ndeg_dens(i,:),xi] = ksdensity(schz.ndeg(i,:),ndeg_x,'function','pdf');
    
    % distance
    temp = dist.*(schz.fdr(:,:,i)<alpha); temp = temp(temp>0);
    schz.dist_mean(i) = mean(temp);                 % average (global) distance
    [schz.dist(i,:),xi] = ksdensity(temp,dist_x);   % distance distribution
    dist_temp = nan(nroi);
    dist_temp(temp>0) = temp(temp>0);
    schz.dist_nod(i,:) = nanmean(dist_temp,2);      % nodal distance
    
end
clear temp_mat full_mat temp temp_b dist_temp

schz.dens = schz.nedge./(nroi*(nroi-1)/2); % edge density

% average degrees of freedom
ctrl.edof_avg = mean(ctrl.edof);
schz.edof_avg = mean(schz.edof);

%% fixed-P thresholding - figures

% full weight distributions
figure;
hold on;
for i = 1:1:nc; hc = plot(wei_x,ctrl.full_wei(i,:),'Color',ctrl.col_light,'LineWidth',lwd_ind); end
for i = 1:1:np; hp = plot(wei_x,schz.full_wei(i,:),'Color',schz.col_light,'LineWidth',lwd_ind); end
hca = plot(wei_x,mean(ctrl.full_wei,1),'Color',ctrl.col*col_ind,'LineWidth',lwd);
hpa = plot(wei_x,mean(schz.full_wei,1),'Color',schz.col*0.65,'LineWidth',lwd);
hold off
set(gca,'FontSize',fsizelab); set(gcf,'color','w'); xlim([-0.6 1]); set(gca,'XTick',[-0.6:0.4:1]);
xlabel('edge weight (unthresholded)','FontSize',fsize,'FontName','Arial');
ylabel('probability density','FontSize',fsize,'FontName','Arial');
print([plot_path 'weight_dist_full_sc' sc '.png'],'-dpng')
clear hc hp hca hpa

% boxplots / histograms

% average edge weight (full)
%ylims = [0 0.7]; yincr = 0.05;     % limits for scale 1
ylims = [0 0.7]; yincr = 0.05;      % limits for scale 2
%ylims = [0 0.8]; yincr = 0.05;     % limits for scale 3
hist_box_plot(ctrl.avg_wei,schz.avg_wei,ylims,yincr,ctrl.col,schz.col,'probability','average edge weight (full)',fsize,fsizeh);
export_fig([plot_path 'boxplot_avg_weight_full_sc' sc '.pdf'],'-pdf','-nocrop')

% edge density
%ylims = [0 0.9]; yincr = 0.1;  % limits for scale 1
ylims = [0 1]; yincr = 0.1;     % limits for scale 2
%ylims = [0 0.6]; yincr = 0.05; % limits for scale 3
hist_box_plot(ctrl.dens,schz.dens,ylims,yincr,ctrl.col,schz.col_light,'probability','edge density',fsize,fsizelab);
export_fig([plot_path 'boxplot_edge_density_sc' sc '.pdf'],'-pdf','-nocrop')

% number of connected components
%ylims = [0 24]; yincr = 2;     % limits for scale 1
ylims = [0 140]; yincr = 10;    % limits for scale 2
%ylims = [0 400]; yincr = 40;   % limits for scale 3
hist_box_plot(ctrl.cmp_n,schz.cmp_n,ylims,yincr,ctrl.col,schz.col,'probability','# connected components',fsize,fsizeh);
export_fig([plot_path 'boxplot_num_conn_comp_sc' sc '.pdf'],'-pdf','-nocrop')

% average degrees of freedom
ylims = [33.5 37.5];
yincr = 0.5;
hist_box_plot(ctrl.edof_avg,schz.edof_avg,ylims,yincr,ctrl.col,schz.col,'probability','# DOF',fsize,fsizeh);
export_fig([plot_path 'boxplot_dof_sc' sc '.pdf'],'-pdf','-nocrop')

% horizontal box plot - avg weight
hor_box_plot(ctrl.avg_wei,schz.avg_wei,[0,0.8],[0:0.2:0.8],ctrl.col,schz.col,'\mu edge weight (unthresholded)','ctrl','schz',fsize,fsizelab)
print([plot_path 'horbox_wei_full_sc' sc '.png'],'-dpng')
%ranksum(ctrl.avg_wei,schz.avg_wei)

%% compare local measures

for i = 1:1:nroi
    
    % edof
    edof_diff(i) = median(ctrl.edof(i,:))-median(schz.edof(i,:));
    edof_p(i) = ranksum(ctrl.edof(i,:),schz.edof(i,:));
    
    % avg wei full
    avg_wei_diff(i) = median(ctrl.avg_wei_nod(:,i))-median(schz.avg_wei_nod(:,i));
    avg_wei_p(i) = ranksum(ctrl.avg_wei_nod(:,i),schz.avg_wei_nod(:,i));
    
    % avg wei thresh
    avg_wei_diff(i) = median(ctrl.avg_wei_nod(:,i))-median(schz.avg_wei_nod(:,i));
    avg_wei_p(i) = ranksum(ctrl.avg_wei_nod(:,i),schz.avg_wei_nod(:,i));
    
    % degree
    deg_diff(i) = median(ctrl.deg(:,i))-median(schz.deg(:,i));
    deg_p(i) = ranksum(ctrl.deg(:,i),schz.deg(:,i));
    
    % normalised degree
    ndeg_diff(i) = median(ctrl.ndeg(:,i))-median(schz.ndeg(:,i));
    ndeg_p(i) = ranksum(ctrl.ndeg(:,i),schz.ndeg(:,i));
    
    % consistency of node in largest component
    cmp_scr_diff(i) = mean(ctrl.cmp_score(:,i))-mean(schz.cmp_score(:,i));
    cmp_scr_p(i) = ranksum(ctrl.cmp_score(:,i),schz.cmp_score(:,i));
    
    % normalised degree
    dist_diff(i) = median(ctrl.dist_nod(:,i))-median(schz.dist_nod(:,i));
    dist_p(i) = ranksum(ctrl.dist_nod(:,i),schz.dist_nod(:,i));  
    
end

loc_sig = 0.01; % significance threshold

% correct for multiple comparisons using FDR
edof_pfdr = mafdr(edof_p,'BHFDR',1);
avg_wei_pfdr = mafdr(avg_wei_p,'BHFDR',1);
avg_wei_pfdr = mafdr(avg_wei_p,'BHFDR',1);
deg_pfdr = mafdr(deg_p,'BHFDR',1);
ndeg_pfdr = mafdr(ndeg_p,'BHFDR',1);
cmp_scr_pfdr = mafdr(cmp_scr_p,'BHFDR',1);
dist_pfdr = mafdr(dist_p,'BHFDR',1);

% % numbers of regions that show differences after adjusting for FDR
% sum(edof_pfdr<loc_sig)      % degrees of freedom
% sum(avg_wei_pfdr<loc_sig)   % average nodal weight
% sum(deg_pfdr<loc_sig)       % degree
% sum(ndeg_pfdr<loc_sig)      % density-normalised degree
% sum(cmp_scr_pfdr<loc_sig)   % nodal component score
% sum(dist_pfdr<loc_sig)      % nodal distance

%% percolation threshold (positive edges only)

% add edges in order of increasing p-value
% identify p-value (and r) at which all nodes become connected

nedge_max = (nroi*(nroi-1)/2);

% controls

% first (for speed) find nearest 1000 edges (nearest to percolation threshold) 
for i = 1:1:nc
    if mod(i,10) == 0; disp(['ctrl ' num2str(i)]); end
    cmp=0;
    temp_p = ctrl.fdr(:,:,i);
    [b,ix] = sort(temp_p(triu_ind),'ascend'); % sort edges by value
    [row,col] = ind2sub([nroi,nroi],triu_ind(ix(1:nedge_max)));
    % fill new matrix
    temp_mat = zeros(nroi);
    for j = 1:nedge_max
        %temp_ind = row((1+(j-1)*100):(j*100))
        if ctrl.r(row(j),col(j),i) >0
            temp_mat(row(j),col(j)) = 1;
            temp_mat(col(j),row(j)) = 1;
        end
        if mod(j,1000) == 0;
            %disp(['edge ' num2str(j)]);
            [cmp,cmp_sizes] = get_components(temp_mat);
            if all(cmp==1);
                ctrl.percol_1000(i) = j;
                break;
            end
        end
        if j==nedge_max
            ctrl.percol_1000(i) = ceil(nedge_max/1000)*1000;
        end
    end
end

% refill matrix to (j-1)*1000 by thousands, then by individual edges
for i = 1:1:nc
    if mod(i,10) == 0; disp(['ctrl ' num2str(i)]); end
    cmp=0;
    temp_p = ctrl.fdr(:,:,i);
    [b,ix] = sort(temp_p(triu_ind),'ascend'); % sort edges by value
    [row,col] = ind2sub([nroi,nroi],triu_ind(ix(1:nedge_max)));
    % fill new matrix
    temp_mat = zeros(nroi);
    for j = 1:(ctrl.percol_1000(i)-1000)            % fill by thousands (w/out verifying connectedness)
        if ctrl.r(row(j),col(j),i)>0
            temp_mat(row(j),col(j)) = 1;
            temp_mat(col(j),row(j)) = 1;
        end
    end
    for j = (ctrl.percol_1000(i)-1000):nedge_max    % fill by individual edges (verifying connectedness)
        if ctrl.r(row(j),col(j),i)>0
            temp_mat(row(j),col(j)) = 1;
            temp_mat(col(j),row(j)) = 1;
            [cmp,cmp_sizes] = get_components(temp_mat);
            if all(cmp==1); break; end
        end
    end
    temp_p = ctrl.fdr(:,:,i).*temp_mat;
    temp_r = ctrl.r(:,:,i).*temp_mat;
    ctrl.percol_p(i) = max(temp_p(:));
    ctrl.percol_r(i) = min(temp_r(temp_r~=0));
    ctrl.percol_nedge(i) = sum(temp_mat(:))/2; % /2 as both triangulars
    ctrl.percol_mat(:,:,i) = temp_mat;
end

ctrl.percol_dens = ctrl.percol_nedge./(nroi*(nroi-1)/2);

% patients

% first (for speed) find nearest 1000 edges (nearest to percolation threshold) 
for i = 1:1:np
    if mod(i,10) == 0; disp(['schz ' num2str(i)]); end
    cmp=0;
    temp_p = schz.fdr(:,:,i); %schz.did(i)
    [b,ix] = sort(temp_p(triu_ind),'ascend'); % sort edges by value
    [row,col] = ind2sub([nroi,nroi],triu_ind(ix(1:nedge_max)));
    % fill new matrix
    temp_mat = zeros(nroi);
    for j = 1:nedge_max
        if schz.r(row(j),col(j),i)>0
            temp_mat(row(j),col(j)) = 1;
            temp_mat(col(j),row(j)) = 1;
        end
        if mod(j,1000) == 0;
            %disp(['edge ' num2str(j)]);
            [cmp,cmp_sizes] = get_components(temp_mat);
            if all(cmp==1);
                schz.percol_1000(i) = j;
                break;
            end 
        end
        if j==nedge_max
            schz.percol_1000(i) = ceil(nedge_max/1000)*1000;
        end
    end
end

% refill matrix to (j-1)*1000 by thousands, then by individual edges
for i = 1:1:np
    if mod(i,10) == 0; disp(['schz ' num2str(i)]); end
    cmp=0;
    temp_p = schz.fdr(:,:,i); %schz.did(i)
    [b,ix] = sort(temp_p(triu_ind),'ascend'); % sort edges by value
    [row,col] = ind2sub([nroi,nroi],triu_ind(ix(1:nedge_max)));
    % fill new matrix
    temp_mat = zeros(nroi);
    for j = 1:(schz.percol_1000(i)-1000)            % fill by thousands (w/out verifying connectedness)
        if schz.r(row(j),col(j),i)>0
            temp_mat(row(j),col(j)) = 1;
            temp_mat(col(j),row(j)) = 1;
        end
    end
    for j = (schz.percol_1000(i)-1000):nedge_max    % fill by individual edges (verifying connectedness)
        if schz.r(row(j),col(j),i)>0
            temp_mat(row(j),col(j)) = 1;
            temp_mat(col(j),row(j)) = 1;
            [cmp,cmp_sizes] = get_components(temp_mat);
            if all(cmp==1); break; end
        end
    end
    temp_p = schz.fdr(:,:,i).*temp_mat;
    temp_r = schz.r(:,:,i).*temp_mat;
    schz.percol_p(i) = max(temp_p(:));
    schz.percol_r(i) = min(temp_r(temp_r~=0));
    schz.percol_nedge(i) = sum(temp_mat(:))/2; % /2 as both triangulars
    schz.percol_mat(:,:,i) = temp_mat;
end

schz.percol_dens = schz.percol_nedge./(nroi*(nroi-1)/2);

% statistics
[p,r] = ranksum_effect_size(ctrl.percol_p,schz.percol_p);
prct_l = 25; prct_h = 75; % percentiles for IQR

%% threshold to fixed density by P (across range of densities)

ndens = 35; % number of edge densities (in steps of 1%)

for d = 1:1:ndens
    disp(['edge density ' num2str(d) '%'])
    
    dens_all = d/100;
    nedge_all = ceil(dens_all*(nroi*(nroi-1)/2));   % number of edges retained after thresholding
    
    % density subjects ID's - subjects whose density (from P<0.01 thresholding) is superior to current edge density
    ctrl.did{d} = find(ctrl.dens > dens_all);
    schz.did{d} = find(schz.dens > dens_all);
    
    % loop over all subjects - later distinguish significant + nonsignificant when comparing topology
    for i = 1:1:nc
        %%% this code ensures that only positive edges are considered
        % if wishing to consider negative edges, replace by commented-out line
        t_p = ctrl.fdr(:,:,i);
        t_r = ctrl.r(:,:,i);
        temp_p = nan(nroi);
        temp_p(t_r>0) = t_p(t_r>0);
        %temp_p = ctrl.fdr(:,:,i); % if wishing to consider negative edges
        %%%
        [b,ix] = sort(temp_p(triu_ind),'ascend');                   % sort edges by value
        [row,col] = ind2sub([nroi,nroi],triu_ind(ix(1:nedge_all))); % indices of edges to be added (up to nedge_all)
        % fill new matrix
        temp_mat_dens = zeros(nroi);
        for j = 1:length(row)
            temp_mat_dens(row(j),col(j)) = ctrl.r(row(j),col(j),i); % weighted
        end   
        ctrl.mat_dens(:,:,d,i) = temp_mat_dens+temp_mat_dens';
        
        % check connectedness
        [cmp,cmp_sizes] = get_components(ctrl.mat_dens(:,:,d,i));
        ctrl.dcmp_n(i,d) = length(unique(cmp));
        ctrl.dcmp_max(i,d) = max(cmp_sizes);
        
        % topological measures
        ctrl.eff(i,d) = efficiency_bin(double(logical(ctrl.mat_dens(:,:,d,i))));
        ctrl.trans(i,d) = transitivity_bu(double(logical(ctrl.mat_dens(:,:,d,i))));
        
        clear temp_p temp_mat_dens
    end
    
    for i = 1:1:np
        %%% this code ensures that only positive edges are considered
        % if wishing to consider negative edges, replace by commented-out line
        t_p = schz.fdr(:,:,i);
        t_r = schz.r(:,:,i);
        temp_p = nan(nroi);
        temp_p(t_r>0) = t_p(t_r>0);
        %temp_p = schz.fdr(:,:,i); % if wishing to consider negative edges
        %%%
        [b,ix] = sort(temp_p(triu_ind),'ascend');                   % sort edges by value
        [row,col] = ind2sub([nroi,nroi],triu_ind(ix(1:nedge_all))); % indices of edges to be added (up to nedge_all)
        % fill new matrix
        temp_mat_dens = zeros(nroi);
        for j = 1:length(row)
            temp_mat_dens(row(j),col(j)) = schz.r(row(j),col(j),i); % weighted
        end
        schz.mat_dens(:,:,d,i) = temp_mat_dens+temp_mat_dens';
        
        % check connectedness
        [cmp,cmp_sizes] = get_components(schz.mat_dens(:,:,d,i));
        schz.dcmp_n(i,d) = length(unique(cmp));
        schz.dcmp_max(i,d) = max(cmp_sizes);
        
        % topological measures
        schz.eff(i,d) = efficiency_bin(double(logical(schz.mat_dens(:,:,d,i))));
        schz.trans(i,d) = transitivity_bu(double(logical(schz.mat_dens(:,:,d,i))));
        
        clear temp_p temp_mat_dens
    end
    
end

%% sets of participants: "significant", and "P<1" 

% edge density when P->1
% ctrl
%ctrl.p1_dens = nan(nc,1);
for i = 1:nc
    %%%
    t_p = ctrl.fdr(:,:,i);
    t_r = ctrl.r(:,:,i);
    temp_p = nan(nroi);
    temp_p(t_r>0) = t_p(t_r>0);
    %%%
    temp_p(temp_p==1)=NaN;
    ctrl.p1_dens(i) = 100*sum(~isnan(temp_p(triu_ind)))/(nroi*(nroi-1)/2);
end

% schz
%schz.p1_dens = nan(np,1);
for i = 1:np
    %%%
    t_p = schz.fdr(:,:,i);
    t_r = schz.r(:,:,i);
    temp_p = nan(nroi);
    temp_p(t_r>0) = t_p(t_r>0);
    %%%
    temp_p(temp_p==1)=NaN;
    schz.p1_dens(i) = 100*sum(~isnan(temp_p(triu_ind)))/(nroi*(nroi-1)/2);
end

for d = 1:1:ndens

    % significant subjects
    ctrl.sig_id{d} = find(ctrl.dens > d/100);
    schz.sig_id{d} = find(schz.dens > d/100);
    
    % subjects with P<1
    ctrl.psub1_id{d} = find(ctrl.p1_dens > d);
    schz.psub1_id{d} = find(schz.p1_dens > d);
    
end

% numbers of subjects in each group as a function of density
ctrl.sig_count = cellfun('length',ctrl.sig_id);
ctrl.psub1_count = cellfun('length',ctrl.psub1_id);

schz.sig_count = cellfun('length',schz.sig_id);
schz.psub1_count = cellfun('length',schz.psub1_id);
    
% % find subject numbers between which ratios of sig/unsig cross 3:1 / 1:3
% %cellfun('length',ctrl.sig_id)
% ctrl.r3to1(1) = min(find(cellfun('length',ctrl.sig_id)<nc*(3/4)))-0.5;
% ctrl.r3to1(2) = min(find(cellfun('length',ctrl.sig_id)<nc*(1/4)))-0.5;
% %cellfun('length',schz.sig_id)
% schz.r3to1(1) = min(find(cellfun('length',schz.sig_id)<np*(3/4)))-0.5;
% schz.r3to1(2) = min(find(cellfun('length',schz.sig_id)<np*(1/4)))-0.5; 

%% drop off curves - include "sig" and "P<1"

figure;
hold on
% "sig"
plot(1:ndens,ctrl.sig_count(1:ndens),'Color',ctrl.col,'LineWidth',1);
plot(1:ndens,schz.sig_count(1:ndens),'Color',schz.col,'LineWidth',1);
scatter(1:ndens,ctrl.sig_count(1:ndens),26,'MarkerFaceColor',ctrl.col,'MarkerEdgeColor',ctrl.col); 
scatter(1:ndens,schz.sig_count(1:ndens),26,'MarkerFaceColor',schz.col,'MarkerEdgeColor',schz.col);
% P<1
plot(1:ndens,ctrl.psub1_count(1:ndens),'-.','Color',ctrl.col,'LineWidth',1);
plot(1:ndens,schz.psub1_count(1:ndens),'-.','Color',schz.col,'LineWidth',1);
hold off

set(gca,'FontSize',fsize2); box off
xlabel('edge density \rho (%)','FontSize',fsize,'FontName','Arial');
ylabel('# participants','FontSize',fsize,'FontName','Arial');

pos = get(gcf,'Position');
set(gcf,'Position', [pos(1) pos(2) 1.25*pos(3) 0.70*pos(4)]); set(gcf,'color','w');

xlim([0 ndens+1]); 
set(gca,'YTick',[0,40,80])

print([plot_path 'drop_off_sc' sc '_alpha' num2str(100*alpha) '.png'],'-dpng')

%% %%%%%%%%%% compare topological measures across densities - P-THRESHOLD %%%%%%%%

% m_...     mean
% ci_...    confidence interval

% eff       efficiency
% trans     transitivity
% dcmp_n    number of (dis)connected components
% dcmp_max  size of the largest component

prct_l = 25;
prct_h = 75;

for i = 1:1:ndens
    
    % all subjects
    % ctrl
    ctrl.m_eff(i) = median(ctrl.eff(:,i));
    ctrl.ci_eff(i,:) = abs(prctile(ctrl.eff(:,i),[prct_l prct_h])-ctrl.m_eff(i));
    ctrl.m_trans(i) = median(ctrl.trans(:,i));
    ctrl.ci_trans(i,:) = abs(prctile(ctrl.trans(:,i),[prct_l prct_h])-ctrl.m_trans(i));
    ctrl.m_dcmp_n(i) = median(ctrl.dcmp_n(:,i));
    ctrl.ci_dcmp_n(i,:) = abs(prctile(ctrl.dcmp_n(:,i),[prct_l prct_h])-ctrl.m_dcmp_n(i));
    ctrl.m_dcmp_max(i) = median(ctrl.dcmp_max(:,i));
    ctrl.ci_dcmp_max(i,:) = abs(prctile(ctrl.dcmp_max(:,i),[prct_l prct_h])-ctrl.m_dcmp_max(i));
    % schz
    schz.m_eff(i) = median(schz.eff(:,i));
    schz.ci_eff(i,:) = abs(prctile(schz.eff(:,i),[prct_l prct_h])-schz.m_eff(i));
    schz.m_trans(i) = median(schz.trans(:,i));
    schz.ci_trans(i,:) = abs(prctile(schz.trans(:,i),[prct_l prct_h])-schz.m_trans(i));
    schz.m_dcmp_n(i) = median(schz.dcmp_n(:,i));
    schz.ci_dcmp_n(i,:) = abs(prctile(schz.dcmp_n(:,i),[prct_l prct_h])-schz.m_dcmp_n(i));
    schz.m_dcmp_max(i) = median(schz.dcmp_max(:,i));
    schz.ci_dcmp_max(i,:) = abs(prctile(schz.dcmp_max(:,i),[prct_l prct_h])-schz.m_dcmp_max(i));
    
    % ranksum across densities
    p_eff(i) = ranksum(ctrl.eff(:,i),schz.eff(:,i));
    p_trans(i) = ranksum(ctrl.trans(:,i),schz.trans(:,i));
    p_dcmp_n(i) = ranksum(ctrl.dcmp_n(:,i),schz.dcmp_n(:,i));
    p_dcmp_max(i) = ranksum(ctrl.dcmp_max(:,i),schz.dcmp_max(:,i));
    
    % significant subjects
    % ctrl
    ctrl.m_eff_sig(i) = median(ctrl.eff(ctrl.sig_id{i},i));
    ctrl.ci_eff_sig(i,:) = abs(prctile(ctrl.eff(ctrl.sig_id{i},i),[prct_l prct_h])-ctrl.m_eff_sig(i));
    ctrl.m_trans_sig(i) = median(ctrl.trans(ctrl.sig_id{i},i));
    ctrl.ci_trans_sig(i,:) = abs(prctile(ctrl.trans(ctrl.sig_id{i},i),[prct_l prct_h])-ctrl.m_trans_sig(i));
    ctrl.m_dcmp_n_sig(i) = median(ctrl.dcmp_n(ctrl.sig_id{i},i));
    ctrl.ci_dcmp_n_sig(i,:) = abs(prctile(ctrl.dcmp_n(ctrl.sig_id{i},i),[prct_l prct_h])-ctrl.m_dcmp_n_sig(i));
    ctrl.m_dcmp_max_sig(i) = median(ctrl.dcmp_max(ctrl.sig_id{i},i));
    ctrl.ci_dcmp_max_sig(i,:) = abs(prctile(ctrl.dcmp_max(ctrl.sig_id{i},i),[prct_l prct_h])-ctrl.m_dcmp_max_sig(i));
    %%%
    ctrl.m_avg_wei_sig(i) = median(ctrl.avg_wei(ctrl.sig_id{i}));
    ctrl.ci_avg_wei_sig(i,:) = abs(prctile(ctrl.avg_wei(ctrl.sig_id{i}),[prct_l prct_h])-ctrl.m_avg_wei_sig(i));
    %%%
    % schz
    schz.m_eff_sig(i) = median(schz.eff(schz.sig_id{i},i));
    schz.ci_eff_sig(i,:) = abs(prctile(schz.eff(schz.sig_id{i},i),[prct_l prct_h])-schz.m_eff_sig(i));
    schz.m_trans_sig(i) = median(schz.trans(schz.sig_id{i},i));
    schz.ci_trans_sig(i,:) = abs(prctile(schz.trans(schz.sig_id{i},i),[prct_l prct_h])-schz.m_trans_sig(i));
    schz.m_dcmp_n_sig(i) = median(schz.dcmp_n(schz.sig_id{i},i));
    schz.ci_dcmp_n_sig(i,:) = abs(prctile(schz.dcmp_n(schz.sig_id{i},i),[prct_l prct_h])-schz.m_dcmp_n_sig(i));
    schz.m_dcmp_max_sig(i) = median(schz.dcmp_max(schz.sig_id{i},i));
    schz.ci_dcmp_max_sig(i,:) = abs(prctile(schz.dcmp_max(schz.sig_id{i},i),[prct_l prct_h])-schz.m_dcmp_max_sig(i));
    %%%
    schz.m_avg_wei_sig(i) = median(schz.avg_wei(schz.sig_id{i}));
    schz.ci_avg_wei_sig(i,:) = abs(prctile(schz.avg_wei(schz.sig_id{i}),[prct_l prct_h])-schz.m_avg_wei_sig(i));
    %%%
    
    % ranksum across densities
    p_eff_sig(i) = ranksum(ctrl.eff(ctrl.sig_id{i},i),schz.eff(schz.sig_id{i},i));
    p_trans_sig(i) = ranksum(ctrl.trans(ctrl.sig_id{i},i),schz.trans(schz.sig_id{i},i));
    p_dcmp_n_sig(i) = ranksum(ctrl.dcmp_n(ctrl.sig_id{i},i),schz.dcmp_n(schz.sig_id{i},i));
    p_dcmp_max_sig(i) = ranksum(ctrl.dcmp_max(ctrl.sig_id{i},i),schz.dcmp_max(schz.sig_id{i},i));
    %%%
    p_avg_wei_sig(i) = ranksum(ctrl.avg_wei(ctrl.sig_id{i}),schz.avg_wei(schz.sig_id{i}));
    %%%
    
    % subjects with P<1
    % ctrl
    ctrl.m_eff_psub1(i) = median(ctrl.eff(ctrl.psub1_id{i},i));
    ctrl.ci_eff_psub1(i,:) = abs(prctile(ctrl.eff(ctrl.psub1_id{i},i),[prct_l prct_h])-ctrl.m_eff_psub1(i));
    ctrl.m_trans_psub1(i) = median(ctrl.trans(ctrl.psub1_id{i},i));
    ctrl.ci_trans_psub1(i,:) = abs(prctile(ctrl.trans(ctrl.psub1_id{i},i),[prct_l prct_h])-ctrl.m_trans_psub1(i));
    ctrl.m_dcmp_n_psub1(i) = median(ctrl.dcmp_n(ctrl.psub1_id{i},i));
    ctrl.ci_dcmp_n_psub1(i,:) = abs(prctile(ctrl.dcmp_n(ctrl.psub1_id{i},i),[prct_l prct_h])-ctrl.m_dcmp_n_psub1(i));
    ctrl.m_dcmp_max_psub1(i) = median(ctrl.dcmp_max(ctrl.psub1_id{i},i));
    ctrl.ci_dcmp_max_psub1(i,:) = abs(prctile(ctrl.dcmp_max(ctrl.psub1_id{i},i),[prct_l prct_h])-ctrl.m_dcmp_max_psub1(i));
    %%%
    ctrl.m_avg_wei_psub1(i) = median(ctrl.avg_wei(ctrl.psub1_id{i}));
    ctrl.ci_avg_wei_psub1(i,:) = abs(prctile(ctrl.avg_wei(ctrl.psub1_id{i}),[prct_l prct_h])-ctrl.m_avg_wei_psub1(i));
    %%%
    % schz
    schz.m_eff_psub1(i) = median(schz.eff(schz.psub1_id{i},i));
    schz.ci_eff_psub1(i,:) = abs(prctile(schz.eff(schz.psub1_id{i},i),[prct_l prct_h])-schz.m_eff_psub1(i));
    schz.m_trans_psub1(i) = median(schz.trans(schz.psub1_id{i},i));
    schz.ci_trans_psub1(i,:) = abs(prctile(schz.trans(schz.psub1_id{i},i),[prct_l prct_h])-schz.m_trans_psub1(i));
    schz.m_dcmp_n_psub1(i) = median(schz.dcmp_n(schz.psub1_id{i},i));
    schz.ci_dcmp_n_psub1(i,:) = abs(prctile(schz.dcmp_n(schz.psub1_id{i},i),[prct_l prct_h])-schz.m_dcmp_n_psub1(i));
    schz.m_dcmp_max_psub1(i) = median(schz.dcmp_max(schz.psub1_id{i},i));
    schz.ci_dcmp_max_psub1(i,:) = abs(prctile(schz.dcmp_max(schz.psub1_id{i},i),[prct_l prct_h])-schz.m_dcmp_max_psub1(i));
    %%%
    schz.m_avg_wei_psub1(i) = median(schz.avg_wei(schz.psub1_id{i}));
    schz.ci_avg_wei_psub1(i,:) = abs(prctile(schz.avg_wei(schz.psub1_id{i}),[prct_l prct_h])-schz.m_avg_wei_psub1(i));
    %%%
    
    % ranksum across densities
    p_eff_psub1(i) = ranksum(ctrl.eff(ctrl.psub1_id{i},i),schz.eff(schz.psub1_id{i},i));
    p_trans_psub1(i) = ranksum(ctrl.trans(ctrl.psub1_id{i},i),schz.trans(schz.psub1_id{i},i));
    p_dcmp_n_psub1(i) = ranksum(ctrl.dcmp_n(ctrl.psub1_id{i},i),schz.dcmp_n(schz.psub1_id{i},i));
    p_dcmp_max_psub1(i) = ranksum(ctrl.dcmp_max(ctrl.psub1_id{i},i),schz.dcmp_max(schz.psub1_id{i},i));
    %%%
    p_avg_wei_psub1(i) = ranksum(ctrl.avg_wei(ctrl.psub1_id{i}),schz.avg_wei(schz.psub1_id{i}));
    %%%
    
    % subjects with P<1 (but P>0.01) (pint = p-intermediate)
    % ctrl
    ctrl.m_eff_pint(i) = median(ctrl.eff(setdiff(ctrl.psub1_id{i},ctrl.sig_id{i}),i));
    ctrl.ci_eff_pint(i,:) = abs(prctile(ctrl.eff(setdiff(ctrl.psub1_id{i},ctrl.sig_id{i}),i),[prct_l prct_h])-ctrl.m_eff_pint(i));
    ctrl.m_trans_pint(i) = median(ctrl.trans(setdiff(ctrl.psub1_id{i},ctrl.sig_id{i}),i));
    ctrl.ci_trans_pint(i,:) = abs(prctile(ctrl.trans(setdiff(ctrl.psub1_id{i},ctrl.sig_id{i}),i),[prct_l prct_h])-ctrl.m_trans_pint(i));
    ctrl.m_dcmp_n_pint(i) = median(ctrl.dcmp_n(setdiff(ctrl.psub1_id{i},ctrl.sig_id{i}),i));
    ctrl.ci_dcmp_n_pint(i,:) = abs(prctile(ctrl.dcmp_n(setdiff(ctrl.psub1_id{i},ctrl.sig_id{i}),i),[prct_l prct_h])-ctrl.m_dcmp_n_pint(i));
    ctrl.m_dcmp_max_pint(i) = median(ctrl.dcmp_max(setdiff(ctrl.psub1_id{i},ctrl.sig_id{i}),i));
    ctrl.ci_dcmp_max_pint(i,:) = abs(prctile(ctrl.dcmp_max(setdiff(ctrl.psub1_id{i},ctrl.sig_id{i}),i),[prct_l prct_h])-ctrl.m_dcmp_max_pint(i));
    %%%
    ctrl.m_avg_wei_pint(i) = median(ctrl.avg_wei(setdiff(ctrl.psub1_id{i},ctrl.sig_id{i})));
    ctrl.ci_avg_wei_pint(i,:) = abs(prctile(ctrl.avg_wei(setdiff(ctrl.psub1_id{i},ctrl.sig_id{i})),[prct_l prct_h])-ctrl.m_avg_wei_pint(i));
    %%%
    % schz
    schz.m_eff_pint(i) = median(schz.eff(setdiff(schz.psub1_id{i},schz.sig_id{i}),i));
    schz.ci_eff_pint(i,:) = abs(prctile(schz.eff(setdiff(schz.psub1_id{i},schz.sig_id{i}),i),[prct_l prct_h])-schz.m_eff_pint(i));
    schz.m_trans_pint(i) = median(schz.trans(setdiff(schz.psub1_id{i},schz.sig_id{i}),i));
    schz.ci_trans_pint(i,:) = abs(prctile(schz.trans(setdiff(schz.psub1_id{i},schz.sig_id{i}),i),[prct_l prct_h])-schz.m_trans_pint(i));
    schz.m_dcmp_n_pint(i) = median(schz.dcmp_n(setdiff(schz.psub1_id{i},schz.sig_id{i}),i));
    schz.ci_dcmp_n_pint(i,:) = abs(prctile(schz.dcmp_n(setdiff(schz.psub1_id{i},schz.sig_id{i}),i),[prct_l prct_h])-schz.m_dcmp_n_pint(i));
    schz.m_dcmp_max_pint(i) = median(schz.dcmp_max(setdiff(schz.psub1_id{i},schz.sig_id{i}),i));
    schz.ci_dcmp_max_pint(i,:) = abs(prctile(schz.dcmp_max(setdiff(schz.psub1_id{i},schz.sig_id{i}),i),[prct_l prct_h])-schz.m_dcmp_max_pint(i));
    %%%
    schz.m_avg_wei_pint(i) = median(schz.avg_wei(setdiff(schz.psub1_id{i},schz.sig_id{i})));
    schz.ci_avg_wei_pint(i,:) = abs(prctile(schz.avg_wei(setdiff(schz.psub1_id{i},schz.sig_id{i})),[prct_l prct_h])-schz.m_avg_wei_pint(i));
    %%%
    
    % non-significant subjects
    % ctrl
    ctrl.m_eff_nonsig(i) = median(ctrl.eff(setdiff(1:nc,ctrl.sig_id{i}),i));
    ctrl.ci_eff_nonsig(i,:) = abs(prctile(ctrl.eff(setdiff(1:nc,ctrl.sig_id{i}),i),[prct_l prct_h])-ctrl.m_eff_nonsig(i));
    ctrl.m_trans_nonsig(i) = median(ctrl.trans(setdiff(1:nc,ctrl.sig_id{i}),i));
    ctrl.ci_trans_nonsig(i,:) = abs(prctile(ctrl.trans(setdiff(1:nc,ctrl.sig_id{i}),i),[prct_l prct_h])-ctrl.m_trans_nonsig(i));
    ctrl.m_dcmp_n_nonsig(i) = median(ctrl.dcmp_n(setdiff(1:nc,ctrl.sig_id{i}),i));
    ctrl.ci_dcmp_n_nonsig(i,:) = abs(prctile(ctrl.dcmp_n(setdiff(1:nc,ctrl.sig_id{i}),i),[prct_l prct_h])-ctrl.m_dcmp_n_nonsig(i));
    ctrl.m_dcmp_max_nonsig(i) = median(ctrl.dcmp_max(setdiff(1:nc,ctrl.sig_id{i}),i));
    ctrl.ci_dcmp_max_nonsig(i,:) = abs(prctile(ctrl.dcmp_max(setdiff(1:nc,ctrl.sig_id{i}),i),[prct_l prct_h])-ctrl.m_dcmp_max_nonsig(i));
    %%%
    ctrl.m_avg_wei_nonsig(i) = median(ctrl.avg_wei(setdiff(1:nc,ctrl.sig_id{i})));
    ctrl.ci_avg_wei_nonsig(i,:) = abs(prctile(ctrl.avg_wei(setdiff(1:nc,ctrl.sig_id{i})),[prct_l prct_h])-ctrl.m_avg_wei_nonsig(i));
    %%%
    % schz
    schz.m_eff_nonsig(i) = median(schz.eff(setdiff(1:np,schz.sig_id{i}),i));
    schz.ci_eff_nonsig(i,:) = abs(prctile(schz.eff(setdiff(1:np,schz.sig_id{i}),i),[prct_l prct_h])-schz.m_eff_nonsig(i));
    schz.m_trans_nonsig(i) = median(schz.trans(setdiff(1:np,schz.sig_id{i}),i));
    schz.ci_trans_nonsig(i,:) = abs(prctile(schz.trans(setdiff(1:np,schz.sig_id{i}),i),[prct_l prct_h])-schz.m_trans_nonsig(i));
    schz.m_dcmp_n_nonsig(i) = median(schz.dcmp_n(setdiff(1:np,schz.sig_id{i}),i));
    schz.ci_dcmp_n_nonsig(i,:) = abs(prctile(schz.dcmp_n(setdiff(1:np,schz.sig_id{i}),i),[prct_l prct_h])-schz.m_dcmp_n_nonsig(i));
    schz.m_dcmp_max_nonsig(i) = median(schz.dcmp_max(setdiff(1:np,schz.sig_id{i}),i));
    schz.ci_dcmp_max_nonsig(i,:) = abs(prctile(schz.dcmp_max(setdiff(1:np,schz.sig_id{i}),i),[prct_l prct_h])-schz.m_dcmp_max_nonsig(i));
    %%%
    schz.m_avg_wei_nonsig(i) = median(schz.avg_wei(setdiff(1:np,schz.sig_id{i})));
    schz.ci_avg_wei_nonsig(i,:) = abs(prctile(schz.avg_wei(setdiff(1:np,schz.sig_id{i})),[prct_l prct_h])-schz.m_avg_wei_nonsig(i));
    %%%

    % ranksum across densities - "sig" VS "intermediate" (NEW)
    if ~isempty(setdiff(ctrl.psub1_id{i},ctrl.sig_id{i}))
        p_eff_ctrl_pint(i) = ranksum(ctrl.eff(setdiff(ctrl.psub1_id{i},ctrl.sig_id{i}),i),ctrl.eff(ctrl.sig_id{i},i));
        p_trans_ctrl_pint(i) = ranksum(ctrl.trans(setdiff(ctrl.psub1_id{i},ctrl.sig_id{i}),i),ctrl.trans(ctrl.sig_id{i},i));
        p_dcmp_n_ctrl_pint(i) = ranksum(ctrl.dcmp_n(setdiff(ctrl.psub1_id{i},ctrl.sig_id{i}),i),ctrl.dcmp_n(ctrl.sig_id{i},i));
        p_dcmp_max_ctrl_pint(i) = ranksum(ctrl.dcmp_max(setdiff(ctrl.psub1_id{i},ctrl.sig_id{i}),i),ctrl.dcmp_max(ctrl.sig_id{i},i));
        p_avg_wei_ctrl_pint(i) = ranksum(ctrl.avg_wei(setdiff(ctrl.psub1_id{i},ctrl.sig_id{i})),ctrl.avg_wei(ctrl.sig_id{i}));
    else
        p_eff_ctrl_pint(i) = NaN;
        p_trans_ctrl_pint(i) = NaN;
        p_dcmp_n_ctrl_pint(i) = NaN;
        p_dcmp_max_ctrl_pint(i) = NaN;
        p_avg_wei_ctrl_pint(i) = NaN;
    end
    
    if ~isempty(setdiff(schz.psub1_id{i},schz.sig_id{i}))
        p_eff_schz_pint(i) = ranksum(schz.eff(setdiff(schz.psub1_id{i},schz.sig_id{i}),i),schz.eff(schz.sig_id{i},i));
        p_trans_schz_pint(i) = ranksum(schz.trans(setdiff(schz.psub1_id{i},schz.sig_id{i}),i),schz.trans(schz.sig_id{i},i));
        p_dcmp_n_schz_pint(i) = ranksum(schz.dcmp_n(setdiff(schz.psub1_id{i},schz.sig_id{i}),i),schz.dcmp_n(schz.sig_id{i},i));
        p_dcmp_max_schz_pint(i) = ranksum(schz.dcmp_max(setdiff(schz.psub1_id{i},schz.sig_id{i}),i),schz.dcmp_max(schz.sig_id{i},i));
        p_avg_wei_schz_pint(i) = ranksum(schz.avg_wei(setdiff(schz.psub1_id{i},schz.sig_id{i})),schz.avg_wei(schz.sig_id{i}));
    else
        p_eff_schz_pint(i) = NaN;
        p_trans_schz_pint(i) = NaN;
        p_dcmp_n_schz_pint(i) = NaN;
        p_dcmp_max_schz_pint(i) = NaN;
        p_avg_wei_schz_pint(i) = NaN;
    end
    
end

%% permutation tests - significance of sig VS P-intermediate (0.01 <= P-int < 1)

nperm = 1000;

% efficiency
p_eff_sig_perm = ones(ndens,nperm); r_eff_sig_perm = ones(ndens,nperm);
p_eff_ctrl_perm = ones(ndens,nperm); r_eff_ctrl_perm = ones(ndens,nperm);
p_eff_schz_perm = ones(ndens,nperm); r_eff_schz_perm = ones(ndens,nperm);

% transitivity
p_trans_sig_perm = ones(ndens,nperm); r_trans_sig_perm = ones(ndens,nperm);
p_trans_ctrl_perm = ones(ndens,nperm); r_trans_ctrl_perm = ones(ndens,nperm);
p_trans_schz_perm = ones(ndens,nperm); r_trans_schz_perm = ones(ndens,nperm);

% # components
p_dcmp_n_sig_perm = ones(ndens,nperm); r_dcmp_n_sig_perm = ones(ndens,nperm);
p_dcmp_n_ctrl_perm = ones(ndens,nperm); r_dcmp_n_ctrl_perm = ones(ndens,nperm);
p_dcmp_n_schz_perm = ones(ndens,nperm); r_dcmp_n_schz_perm = ones(ndens,nperm);

% size of largest components
p_dcmp_max_sig_perm = ones(ndens,nperm); r_dcmp_max_sig_perm = ones(ndens,nperm);
p_dcmp_max_ctrl_perm = ones(ndens,nperm); r_dcmp_max_ctrl_perm = ones(ndens,nperm);
p_dcmp_max_schz_perm = ones(ndens,nperm); r_dcmp_max_schz_perm = ones(ndens,nperm);

%%%
% average weight
p_eff_sig_perm = ones(ndens,nperm); r_eff_sig_perm = ones(ndens,nperm);
p_eff_ctrl_perm = ones(ndens,nperm); r_eff_ctrl_perm = ones(ndens,nperm);
p_eff_schz_perm = ones(ndens,nperm); r_eff_schz_perm = ones(ndens,nperm);
%%%

p_eff_ctrl_pint(i) = ranksum(ctrl.eff(setdiff(ctrl.psub1_id{i},ctrl.sig_id{i}),i),ctrl.eff(ctrl.sig_id{i},i));

for d = 2:1:ndens % ignore d=1 as all controls are in the significant group (and thus the P-int group is empty)
    d
    
    % "empirical" effect sizes
    
    % efficiency
    [p_eff_sig_emp(d),r_eff_sig_emp(d)] = ranksum_effect_size(ctrl.eff(ctrl.sig_id{d},d),schz.eff(schz.sig_id{d},d)); % ctrl VS schz (sig only)
    [p_eff_ctrl_emp(d),r_eff_ctrl_emp(d)] = ranksum_effect_size(ctrl.eff(ctrl.sig_id{d},d),ctrl.eff(setdiff(ctrl.psub1_id{d},ctrl.sig_id{d}),d)); % within group (sig VS nonsig)
    [p_eff_schz_emp(d),r_eff_schz_emp(d)] = ranksum_effect_size(schz.eff(schz.sig_id{d},d),schz.eff(setdiff(schz.psub1_id{d},schz.sig_id{d}),d));
    
    % transitivity
    [p_trans_sig_emp(d),r_trans_sig_emp(d)] = ranksum_effect_size(ctrl.trans(ctrl.sig_id{d},d),schz.trans(schz.sig_id{d},d)); % ctrl VS schz (sig only)
    [p_trans_ctrl_emp(d),r_trans_ctrl_emp(d)] = ranksum_effect_size(ctrl.trans(ctrl.sig_id{d},d),ctrl.trans(setdiff(ctrl.psub1_id{d},ctrl.sig_id{d}),d)); % within group (sig VS nonsig)
    [p_trans_schz_emp(d),r_trans_schz_emp(d)] = ranksum_effect_size(schz.trans(schz.sig_id{d},d),schz.trans(setdiff(schz.psub1_id{d},schz.sig_id{d}),d));
    
    % # components
    [p_dcmp_n_sig_emp(d),r_dcmp_n_sig_emp(d)] = ranksum_effect_size(ctrl.dcmp_n(ctrl.sig_id{d},d),schz.dcmp_n(schz.sig_id{d},d)); % ctrl VS schz (sig only)
    [p_dcmp_n_ctrl_emp(d),r_dcmp_n_ctrl_emp(d)] = ranksum_effect_size(ctrl.dcmp_n(ctrl.sig_id{d},d),ctrl.dcmp_n(setdiff(ctrl.psub1_id{d},ctrl.sig_id{d}),d)); % within group (sig VS nonsig)
    [p_dcmp_n_schz_emp(d),r_dcmp_n_schz_emp(d)] = ranksum_effect_size(schz.dcmp_n(schz.sig_id{d},d),schz.dcmp_n(setdiff(schz.psub1_id{d},schz.sig_id{d}),d));
    
    % size of largest component
    [p_dcmp_max_sig_emp(d),r_dcmp_max_sig_emp(d)] = ranksum_effect_size(ctrl.dcmp_max(ctrl.sig_id{d},d),schz.dcmp_max(schz.sig_id{d},d)); % ctrl VS schz (sig only)
    [p_dcmp_max_ctrl_emp(d),r_dcmp_max_ctrl_emp(d)] = ranksum_effect_size(ctrl.dcmp_max(ctrl.sig_id{d},d),ctrl.dcmp_max(setdiff(ctrl.psub1_id{d},ctrl.sig_id{d}),d)); % within group (sig VS nonsig)
    [p_dcmp_max_schz_emp(d),r_dcmp_max_schz_emp(d)] = ranksum_effect_size(schz.dcmp_max(schz.sig_id{d},d),schz.dcmp_max(setdiff(schz.psub1_id{d},schz.sig_id{d}),d));
    
    % average weight
    [p_avg_wei_sig_emp(d),r_avg_wei_sig_emp(d)] = ranksum_effect_size(ctrl.avg_wei_full(ctrl.sig_id{d}),schz.avg_wei_full(schz.sig_id{d})); % ctrl VS schz (sig only)
    [p_avg_wei_ctrl_emp(d),r_avg_wei_ctrl_emp(d)] = ranksum_effect_size(ctrl.avg_wei_full(ctrl.sig_id{d}),ctrl.avg_wei_full(setdiff(ctrl.psub1_id{d},ctrl.sig_id{d}))); % within group (sig VS nonsig)
    [p_avg_wei_schz_emp(d),r_avg_wei_schz_emp(d)] = ranksum_effect_size(schz.avg_wei_full(schz.sig_id{d}),schz.avg_wei_full(setdiff(schz.psub1_id{d},schz.sig_id{d})));
    
    % permuted effect sizes
    for j = 1:nperm
        
        % permuted ids (length = sig)
        perm_idc = randsample(nc,length(ctrl.sig_id{d}));
        perm_idp = randsample(np,length(schz.sig_id{d}));
        
        % permuted ids P-INT "complement" (from complement of #p-sig, sample #p-int)
        perm_idc_comp = randsample(setdiff(1:nc,perm_idc),length(ctrl.psub1_id{d})-length(ctrl.sig_id{d}));
        perm_idp_comp = randsample(setdiff(1:np,perm_idp),length(ctrl.psub1_id{d})-length(ctrl.sig_id{d}));
        
        % efficiency
        [p_eff_sig_perm(d,j),r_eff_sig_perm(d,j)] = ranksum_effect_size(ctrl.eff(perm_idc,d),schz.eff(perm_idp,d)); % ctrl VS schz (sig only size)
        [p_eff_ctrl_perm(d,j),r_eff_ctrl_perm(d,j)] = ranksum_effect_size(ctrl.eff(perm_idc,d),ctrl.eff(perm_idc_comp,d)); % within group (sig VS nonsig size)
        [p_eff_schz_perm(d,j),r_eff_schz_perm(d,j)] = ranksum_effect_size(schz.eff(perm_idp,d),schz.eff(perm_idp_comp,d));
        
        % transitivity
        [p_trans_sig_perm(d,j),r_trans_sig_perm(d,j)] = ranksum_effect_size(ctrl.trans(perm_idc,d),schz.trans(perm_idp,d)); % ctrl VS schz (sig only size)
        [p_trans_ctrl_perm(d,j),r_trans_ctrl_perm(d,j)] = ranksum_effect_size(ctrl.trans(perm_idc,d),ctrl.trans(perm_idc_comp,d)); % within group (sig VS nonsig size)
        [p_trans_schz_perm(d,j),r_trans_schz_perm(d,j)] = ranksum_effect_size(schz.trans(perm_idp,d),schz.trans(perm_idp_comp,d));
        
        % # components
        [p_dcmp_n_sig_perm(d,j),r_dcmp_n_sig_perm(d,j)] = ranksum_effect_size(ctrl.dcmp_n(perm_idc,d),schz.dcmp_n(perm_idp,d)); % ctrl VS schz (sig only size)
        [p_dcmp_n_ctrl_perm(d,j),r_dcmp_n_ctrl_perm(d,j)] = ranksum_effect_size(ctrl.dcmp_n(perm_idc,d),ctrl.dcmp_n(perm_idc_comp,d)); % within group (sig VS nonsig size)
        [p_dcmp_n_schz_perm(d,j),r_dcmp_n_schz_perm(d,j)] = ranksum_effect_size(schz.dcmp_n(perm_idp,d),schz.dcmp_n(perm_idp_comp,d));
        
        % size of largest component
        [p_dcmp_max_sig_perm(d,j),r_dcmp_max_sig_perm(d,j)] = ranksum_effect_size(ctrl.dcmp_max(perm_idc,d),schz.dcmp_max(perm_idp,d)); % ctrl VS schz (sig only size)
        [p_dcmp_max_ctrl_perm(d,j),r_dcmp_max_ctrl_perm(d,j)] = ranksum_effect_size(ctrl.dcmp_max(perm_idc,d),ctrl.dcmp_max(perm_idc_comp,d)); % within group (sig VS nonsig size)
        [p_dcmp_max_schz_perm(d,j),r_dcmp_max_schz_perm(d,j)] = ranksum_effect_size(schz.dcmp_max(perm_idp,d),schz.dcmp_max(perm_idp_comp,d));

        % average weight
        [p_avg_wei_sig_perm(d,j),r_avg_wei_sig_perm(d,j)] = ranksum_effect_size(ctrl.avg_wei_full(perm_idc),schz.avg_wei_full(perm_idp)); % ctrl VS schz (sig only size)
        [p_avg_wei_ctrl_perm(d,j),r_avg_wei_ctrl_perm(d,j)] = ranksum_effect_size(ctrl.avg_wei_full(perm_idc),ctrl.avg_wei_full(perm_idc_comp)); % within group (sig VS nonsig size)
        [p_avg_wei_schz_perm(d,j),r_avg_wei_schz_perm(d,j)] = ranksum_effect_size(schz.avg_wei_full(perm_idp),schz.avg_wei_full(perm_idp_comp));
        
    end
    
    % proportions of ranksum r's greater in the empirical case
    
    % efficiency
    pperm_eff_sig(d) = sum(abs(r_eff_sig_emp(d))<abs(r_eff_sig_perm(d,:)))/nperm;
    pperm_eff_ctrl(d) = sum(abs(r_eff_ctrl_emp(d))<abs(r_eff_ctrl_perm(d,:)))/nperm;
    pperm_eff_schz(d) = sum(abs(r_eff_schz_emp(d))<abs(r_eff_schz_perm(d,:)))/nperm;
    
    % transitivity
    pperm_trans_sig(d) = sum(abs(r_trans_sig_emp(d))<abs(r_trans_sig_perm(d,:)))/nperm;
    pperm_trans_ctrl(d) = sum(abs(r_trans_ctrl_emp(d))<abs(r_trans_ctrl_perm(d,:)))/nperm;
    pperm_trans_schz(d) = sum(abs(r_trans_schz_emp(d))<abs(r_trans_schz_perm(d,:)))/nperm;
    
    % # components
    pperm_dcmp_n_sig(d) = sum(abs(r_dcmp_n_sig_emp(d))<abs(r_dcmp_n_sig_perm(d,:)))/nperm;
    pperm_dcmp_n_ctrl(d) = sum(abs(r_dcmp_n_ctrl_emp(d))<abs(r_dcmp_n_ctrl_perm(d,:)))/nperm;
    pperm_dcmp_n_schz(d) = sum(abs(r_dcmp_n_schz_emp(d))<abs(r_dcmp_n_schz_perm(d,:)))/nperm;
    
    % size of largest component
    pperm_dcmp_max_sig(d) = sum(abs(r_dcmp_max_sig_emp(d))<abs(r_dcmp_max_sig_perm(d,:)))/nperm;
    pperm_dcmp_max_ctrl(d) = sum(abs(r_dcmp_max_ctrl_emp(d))<abs(r_dcmp_max_ctrl_perm(d,:)))/nperm;
    pperm_dcmp_max_schz(d) = sum(abs(r_dcmp_max_schz_emp(d))<abs(r_dcmp_max_schz_perm(d,:)))/nperm;
    
    % average weight
    pperm_avg_wei_sig(d) = sum(abs(r_avg_wei_sig_emp(d))<abs(r_avg_wei_sig_perm(d,:)))/nperm;
    pperm_avg_wei_ctrl(d) = sum(abs(r_avg_wei_ctrl_emp(d))<abs(r_avg_wei_ctrl_perm(d,:)))/nperm;
    pperm_avg_wei_schz(d) = sum(abs(r_avg_wei_schz_emp(d))<abs(r_avg_wei_schz_perm(d,:)))/nperm;

end

%% obtain maximum p-value if keeping all subjects

maxp_dens = ndens;

% delete relevant variables
clear_fields = {'max_p'};
if all(isfield(ctrl,clear_fields)) % Will be True or False.
    ctrl = rmfield(ctrl,clear_fields);
    schz = rmfield(schz,clear_fields);
end

for d = 1:1:maxp_dens
    d
    dens_all = d/100;
    nedge_all = ceil(dens_all*(nroi*(nroi-1)/2));
    
    for i = 1:1:nc
        %%% exclude negative edges
        t_p = ctrl.fdr(:,:,i);
        t_r = ctrl.r(:,:,i);
        temp_p = nan(nroi);
        temp_p(t_r>0) = t_p(t_r>0);
        %%%
        [b,ix] = sort(temp_p(triu_ind),'ascend'); % sort edges by value
        [row,col] = ind2sub([nroi,nroi],triu_ind(ix(1:nedge_all)));
        % fill new matrix
        temp_mat_dens = zeros(nroi);
        for j = 1:length(row)
            temp_mat_dens(row(j),col(j)) = ctrl.r(row(j),col(j)); % weighted
        end
        ctrl.max_p(d,i) = max(max(temp_p(logical(temp_mat_dens))));
        clear temp_mat temp_mat_dens
    end
    
    for i = 1:1:np
        %%% exclude negative edges
        t_p = schz.fdr(:,:,i);
        t_r = schz.r(:,:,i);
        temp_p = nan(nroi);
        temp_p(t_r>0) = t_p(t_r>0);
        %%%
        [b,ix] = sort(temp_p(triu_ind),'ascend'); % sort edges by value
        [row,col] = ind2sub([nroi,nroi],triu_ind(ix(1:nedge_all)));
        % fill new matrix
        temp_mat_dens = zeros(nroi);
        for j = 1:length(row)
            temp_mat_dens(row(j),col(j)) = schz.r(row(j),col(j)); % weighted
        end
        schz.max_p(d,i) = max(max(temp_p(logical(temp_mat_dens))));
        clear temp_mat temp_mat_dens
    end
    
end

% plot, as a function of edge density
figure;
hold on
for i = 1:1:nc; hc = plot(1:maxp_dens,ctrl.max_p(:,i),'Color',ctrl.col_light,'LineWidth',lwd_ind); end
for i = 1:1:np; hp = plot(1:maxp_dens,schz.max_p(:,i),'Color',schz.col_light,'LineWidth',lwd_ind); end
hca = plot(1:maxp_dens,mean(ctrl.max_p,2),'Color',ctrl.col*col_ind,'LineWidth',lwd);
hpa = plot(1:maxp_dens,mean(schz.max_p,2),'Color',schz.col*0.8,'LineWidth',lwd);
hold off

set(gca,'FontSize',fsize2); box off
xlabel('edge density \rho (%)','FontSize',fsize,'FontName','Arial');
ylabel('max(P_{FDR}) / participant','FontSize',fsize,'FontName','Arial');

pos = get(gcf,'Position');
set(gcf,'Position', [pos(1) pos(2) 1.25*pos(3) 0.70*pos(4)]); set(gcf,'color','w');

xlim([0 ndens+1]); 
set(gca,'YTick',[0,0.5,1])

%align_Ylabels(f);

print([plot_path 'max_P_sc' sc '_alpha' num2str(100*alpha) '.png'],'-dpng')

%% histogram bars for nparticipant ratios (n(sig) VS n(non-sig))

% ctrl
figure
b = bar(1:ndens,[cellfun('length',ctrl.sig_id);cellfun('length',ctrl.psub1_id)-cellfun('length',ctrl.sig_id); nc-cellfun('length',ctrl.psub1_id)]','stack');
set(b(1),'FaceColor',ctrl.col,'FaceAlpha',0.7)
set(b(2),'FaceColor',ctrl.col_nonsig,'FaceAlpha',0.7)
set(b(3),'FaceColor',rgb('lightgrey'),'FaceAlpha',0.7)
xlim([0 ndens+1]); ylim([0 nc]); set(gca,'YTick',[0,nc]);
set(gca,'FontSize',fsize2); 
xlabel('edge density \rho (%)','FontSize',fsize,'FontName','Arial');
ylabel('# ctrl','FontSize',fsize,'FontName','Arial');
pos = get(gcf,'Position');
set(gcf,'Position', [pos(1) pos(2) 1.25*pos(3) 0.35*pos(4)]); set(gcf,'color','w');
print([plot_path 'ctrl_sig_vs_nonsig_bar_sc' sc '_alpha' num2str(100*alpha) '.png'],'-dpng')

% schz
figure
b = bar(1:ndens,[cellfun('length',schz.sig_id);cellfun('length',schz.psub1_id)-cellfun('length',schz.sig_id); np-cellfun('length',schz.psub1_id)]','stack');
set(b(1),'FaceColor',schz.col,'FaceAlpha',0.7)
set(b(2),'FaceColor',schz.col_nonsig,'FaceAlpha',0.7)
set(b(3),'FaceColor',rgb('lightgrey'),'FaceAlpha',0.7)
xlim([0 ndens+1]); ylim([0 np]); set(gca,'YTick',[0,np]);
set(gca,'FontSize',fsize2); 
xlabel('edge density \rho (%)','FontSize',fsize,'FontName','Arial');
ylabel('# schz','FontSize',fsize,'FontName','Arial');
pos = get(gcf,'Position');
set(gcf,'Position', [pos(1) pos(2) 1.25*pos(3) 0.35*pos(4)]); set(gcf,'color','w');
print([plot_path 'schz_sig_vs_nonsig_bar_sc' sc '_alpha' num2str(100*alpha) '.png'],'-dpng')

%% plots - empirical (eff, trans)

fsize = 30;
fsize2 = 20;

%%%%% ctrl VS schz

%%% all (P<1) subjects - with drop-off

% efficiency
boundedline_plot_2p_drop(1:ndens,ctrl.m_eff_psub1,ctrl.ci_eff_psub1,schz.m_eff_psub1,schz.ci_eff_psub1,ctrl.psub1_count(1:ndens),schz.psub1_count(1:ndens),...
    ctrl.col,schz.col,alph,NaN,p_eff_psub1,[],'edge density \rho (%)','efficiency','# part. (P<1)','ctrl','schz',[0 ndens+1],[0 0.8],[0 80],fsize,fsize2,...
    plot_path,['eff_ctrl_vs_schz_psub1_boundedline_sc' sc '_alpha' num2str(100*alpha)]);

% transitivity
boundedline_plot_2p_drop(1:ndens,ctrl.m_trans_psub1,ctrl.ci_trans_psub1,schz.m_trans_psub1,schz.ci_trans_psub1,ctrl.psub1_count(1:ndens),schz.psub1_count(1:ndens),...
    ctrl.col,schz.col,alph,NaN,p_trans_psub1,[],'edge density \rho (%)','transitivity','# part. (P<1)','ctrl','schz',[0 ndens+1],[0.3 0.9],[0 80],fsize,fsize2,...
    plot_path,['trans_ctrl_vs_schz_psub1_boundedline_sc' sc '_alpha' num2str(100*alpha)]);

%%% significant subjects only

% efficiency
boundedline_plot_2p_drop(1:ndens,ctrl.m_eff_sig,ctrl.ci_eff_sig,schz.m_eff_sig,schz.ci_eff_sig,ctrl.sig_count(1:ndens),schz.sig_count(1:ndens),...
    ctrl.col,schz.col,alph,NaN,p_eff_sig,pperm_eff_sig,'edge density \rho (%)','efficiency','# sig. participants','ctrl (sig)','schz (sig)',[0 ndens+1],[0 0.8],[0 80],fsize,fsize2,...
    plot_path,['eff_ctrl_vs_schz_sig_boundedline_sc' sc '_alpha' num2str(100*alpha)]);

% transitivity
boundedline_plot_2p_drop(1:ndens,ctrl.m_trans_sig,ctrl.ci_trans_sig,schz.m_trans_sig,schz.ci_trans_sig,ctrl.sig_count(1:ndens),schz.sig_count(1:ndens),...
    ctrl.col,schz.col,alph,NaN,p_trans_sig,pperm_trans_sig,'edge density \rho (%)','transitivity','# sig. participants','ctrl (sig)','schz (sig)',[0 ndens+1],[0.3 0.9],[0 80],fsize,fsize2,...
    plot_path,['trans_ctrl_vs_schz_sig_boundedline_sc' sc '_alpha' num2str(100*alpha)]);

%%%%% significant VS non-significant (P-int: 0.01 <= P-int < 1) subjects 

%%% ctrl

% efficiency
boundedline_plot_2p(1:ndens,ctrl.m_eff_sig,ctrl.ci_eff_sig,ctrl.m_eff_pint,ctrl.ci_eff_pint,...
    ctrl.col,ctrl.col_nonsig,alph,ctrl.r3to1,p_eff_ctrl_pint,pperm_eff_ctrl,'edge density \rho (%)','efficiency','significant','non-significant',[0 ndens+1],[0 0.8],fsize,fsize2,...
    plot_path,['eff_sig_vs_pint_ctrl_boundedline_sc' sc '_alpha' num2str(100*alpha)]);

% transitivity
boundedline_plot_2p(1:ndens,ctrl.m_trans_sig,ctrl.ci_trans_sig,ctrl.m_trans_pint,ctrl.ci_trans_pint,...
    ctrl.col,ctrl.col_nonsig,alph,ctrl.r3to1,p_trans_ctrl_pint,pperm_trans_ctrl,'edge density \rho (%)','transitivity','significant','non-significant',[0 ndens+1],[0.3 0.9],fsize,fsize2,...
    plot_path,['trans_sig_vs_pint_ctrl_boundedline_sc' sc '_alpha' num2str(100*alpha)]);

%%% schz

% efficiency
boundedline_plot_2p(1:ndens,schz.m_eff_sig,schz.ci_eff_sig,schz.m_eff_pint,schz.ci_eff_pint,...
    schz.col,schz.col_nonsig,alph,schz.r3to1,p_eff_schz_pint,pperm_eff_schz,'edge density \rho (%)','efficiency','significant','non-significant',[0 ndens+1],[0 0.8],fsize,fsize2,...
    plot_path,['eff_sig_vs_pint_schz_boundedline_sc' sc '_alpha' num2str(100*alpha)]);

% transitivity
boundedline_plot_2p(1:ndens,schz.m_trans_sig,schz.ci_trans_sig,schz.m_trans_pint,schz.ci_trans_pint,...
    schz.col,schz.col_nonsig,alph,schz.r3to1,p_trans_schz_pint,pperm_trans_schz,'edge density \rho (%)','transitivity','significant','non-significant',[0 ndens+1],[0.3 0.9],fsize,fsize2,...
    plot_path,['trans_sig_vs_pint_schz_boundedline_sc' sc '_alpha' num2str(100*alpha)]);

%% plots - emp (cmp, mean corr)

%%%%% ctrl VS schz

%%% all (P<1) subjects - with drop-off

% # components
boundedline_plot_2p_drop(1:ndens,ctrl.m_dcmp_n_psub1,ctrl.ci_dcmp_n_psub1,schz.m_dcmp_n_psub1,schz.ci_dcmp_n_psub1,ctrl.psub1_count(1:ndens),schz.psub1_count(1:ndens),...
    ctrl.col,schz.col,alph,NaN,p_dcmp_n_psub1,[],'edge density \rho (%)','# components','# part. (P<1)','ctrl','schz',[0 ndens+1],[0 150],[0 80],fsize,fsize2,...
    plot_path,['dcmp_n_ctrl_vs_schz_psub1_boundedline_sc' sc '_alpha' num2str(100*alpha)]);

% mean correlation
boundedline_plot_2p_drop(1:ndens,ctrl.m_avg_wei_psub1,ctrl.ci_avg_wei_psub1,schz.m_avg_wei_psub1,schz.ci_avg_wei_psub1,ctrl.psub1_count(1:ndens),schz.psub1_count(1:ndens),...
    ctrl.col,schz.col,alph,NaN,p_avg_wei_psub1,[],'edge density \rho (%)','\mu correlation','# part. (P<1)','ctrl','schz',[0 ndens+1],[0 0.8],[0 80],fsize,fsize2,...
    plot_path,['avg_wei_ctrl_vs_schz_psub1_boundedline_sc' sc '_alpha' num2str(100*alpha)]);

%%% significant subjects only

% # components
boundedline_plot_2p_drop(1:ndens,ctrl.m_dcmp_n_sig,ctrl.ci_dcmp_n_sig,schz.m_dcmp_n_sig,schz.ci_dcmp_n_sig,ctrl.sig_count(1:ndens),schz.sig_count(1:ndens),...
    ctrl.col,schz.col,alph,NaN,p_dcmp_n_sig,pperm_dcmp_n_sig,'edge density \rho (%)','# components','# sig. participants','ctrl (sig)','schz (sig)',[0 ndens+1],[0 150],[0 80],fsize,fsize2,...
    plot_path,['dcmp_n_ctrl_vs_schz_sig_boundedline_sc' sc '_alpha' num2str(100*alpha)]);

% mean correlation
boundedline_plot_2p_drop(1:ndens,ctrl.m_avg_wei_sig,ctrl.ci_avg_wei_sig,schz.m_avg_wei_sig,schz.ci_avg_wei_sig,ctrl.sig_count(1:ndens),schz.sig_count(1:ndens),...
    ctrl.col,schz.col,alph,NaN,p_avg_wei_sig,pperm_avg_wei_sig,'edge density \rho (%)','\mu correlation','# sig. participants',[],[],[0 ndens+1],[0.2 0.6],[0 80],fsize,fsize2,...
    plot_path,['avg_wei_ctrl_vs_schz_sig_boundedline_sc' sc '_alpha' num2str(100*alpha)]);

%%%%% significant VS non-significant (P-int: 0.01 <= P-int < 1) subjects 

%%% ctrl

% # components
boundedline_plot_2p(1:ndens,ctrl.m_dcmp_n_sig,ctrl.ci_dcmp_n_sig,ctrl.m_dcmp_n_pint,ctrl.ci_dcmp_n_pint,...
    ctrl.col,ctrl.col_nonsig,alph,ctrl.r3to1,p_dcmp_n_ctrl_pint,pperm_dcmp_n_ctrl,'edge density \rho (%)','# components','significant','non-significant',[0 ndens+1],[0 150],fsize,fsize2,...
    plot_path,['dcmp_n_sig_vs_pint_ctrl_boundedline_sc' sc '_alpha' num2str(100*alpha)]);

% mean correlation
boundedline_plot_2p(1:ndens,ctrl.m_avg_wei_sig,ctrl.ci_avg_wei_sig,ctrl.m_avg_wei_pint,ctrl.ci_avg_wei_pint,...
    ctrl.col,ctrl.col_nonsig,alph,ctrl.r3to1,p_avg_wei_ctrl_pint,pperm_avg_wei_ctrl,'edge density \rho (%)','\mu correlation',[],[],[0 ndens+1],[0 0.6],fsize,fsize2,...
    plot_path,['avg_wei_sig_vs_nonsig_ctrl_boundedline_sc' sc '_alpha' num2str(100*alpha)]);

%%% schz

% # components
boundedline_plot_2p(1:ndens,schz.m_dcmp_n_sig,schz.ci_dcmp_n_sig,schz.m_dcmp_n_pint,schz.ci_dcmp_n_pint,...
    schz.col,schz.col_nonsig,alph,schz.r3to1,p_dcmp_n_schz_pint,pperm_dcmp_n_schz,'edge density \rho (%)','# components','significant','non-significant',[0 ndens+1],[0 150],fsize,fsize2,...
    plot_path,['dcmp_n_sig_vs_pint_schz_boundedline_sc' sc '_alpha' num2str(100*alpha)]);

% mean correlation
boundedline_plot_2p(1:ndens,schz.m_avg_wei_sig,schz.ci_avg_wei_sig,schz.m_avg_wei_pint,schz.ci_avg_wei_pint,...
    schz.col,schz.col_nonsig,alph,schz.r3to1,p_avg_wei_schz_pint,pperm_avg_wei_schz,'edge density \rho (%)','\mu correlation',[],[],[0 ndens+1],[0 0.6],fsize,fsize2,...
    plot_path,['avg_wei_sig_vs_nonsig_schz_boundedline_sc' sc '_alpha' num2str(100*alpha)]);

%% threshold to fixed density by R (across range of densities)

%ndens = floor(100*max(max(ctrl.dens),max(schz.dens)));
ndens = 35;

for d = 1:1:ndens % density cut-off %ndens
    d
    dens_all = d/100;
    nedge_all = ceil(dens_all*(nroi*(nroi-1)/2));
    
    % loop over all subjects - later distinguish sig + nonsig when comparing topology
    for i = 1:1:nc
        temp_r = ctrl.r(:,:,i);
        [b,ix] = sort(temp_r(triu_ind),'descend'); % sort edges by value
        [row,col] = ind2sub([nroi,nroi],triu_ind(ix(1:nedge_all)));
        % fill new matrix
        temp_mat_dens = zeros(nroi);
        for j = 1:length(row)
            temp_mat_dens(row(j),col(j)) = ctrl.r(row(j),col(j),i); % weighted
        end   
        ctrl.mat_dens_r(:,:,d,i) = temp_mat_dens+temp_mat_dens';
        
        % check connectedness
        [cmp,cmp_sizes] = get_components(ctrl.mat_dens_r(:,:,d,i));
        ctrl.dcmp_n_r(i,d) = length(unique(cmp));
        ctrl.dcmp_max_r(i,d) = max(cmp_sizes);

        ctrl.eff_r(i,d) = efficiency_bin(double(logical(ctrl.mat_dens_r(:,:,d,i))));
        ctrl.trans_r(i,d) = transitivity_bu(double(logical(ctrl.mat_dens_r(:,:,d,i))));
    end
    
    for i = 1:1:np
        temp_r = schz.r(:,:,i);
        [b,ix] = sort(temp_r(triu_ind),'descend'); % sort edges by value
        [row,col] = ind2sub([nroi,nroi],triu_ind(ix(1:nedge_all)));
        % fill new matrix
        temp_mat_dens = zeros(nroi);
        for j = 1:length(row)
            temp_mat_dens(row(j),col(j)) = schz.r(row(j),col(j),i); % weighted
        end
        schz.mat_dens_r(:,:,d,i) = temp_mat_dens+temp_mat_dens';
        
        % check connectedness
        [cmp,cmp_sizes] = get_components(schz.mat_dens_r(:,:,d,i));
        schz.dcmp_n_r(i,d) = length(unique(cmp));
        schz.dcmp_max_r(i,d) = max(cmp_sizes);
        
        schz.eff_r(i,d) = efficiency_bin(double(logical(schz.mat_dens_r(:,:,d,i))));
        schz.trans_r(i,d) = transitivity_bu(double(logical(schz.mat_dens_r(:,:,d,i))));
    end
    
end

%% %%%%%%%%%% compare topological measures across densities - R-THRESHOLD %%%%%%%%

% m_...     mean
% ci_...    confidence interval

% eff       efficiency
% trans     transitivity
% dcmp_n    number of (dis)connected components
% dcmp_max  size of the largest component

% for confidence intervals (absolute value taken, as required for plotting)
prct_l = 25; % lower percentile
prct_h = 75; % upper percentile   

prct_l = 25;
prct_h = 75;

for i = 1:1:ndens
    
    % all subjects
    % ctrl
    ctrl.m_eff_r(i) = median(ctrl.eff_r(:,i));
    ctrl.ci_eff_r(i,:) = abs(prctile(ctrl.eff_r(:,i),[prct_l prct_h])-ctrl.m_eff_r(i));
    ctrl.m_trans_r(i) = median(ctrl.trans_r(:,i));
    ctrl.ci_trans_r(i,:) = abs(prctile(ctrl.trans_r(:,i),[prct_l prct_h])-ctrl.m_trans_r(i));
    ctrl.m_dcmp_n_r(i) = median(ctrl.dcmp_n_r(:,i));
    ctrl.ci_dcmp_n_r(i,:) = abs(prctile(ctrl.dcmp_n_r(:,i),[prct_l prct_h])-ctrl.m_dcmp_n_r(i));
    ctrl.m_dcmp_max_r(i) = median(ctrl.dcmp_max_r(:,i));
    ctrl.ci_dcmp_max_r(i,:) = abs(prctile(ctrl.dcmp_max_r(:,i),[prct_l prct_h])-ctrl.m_dcmp_max_r(i));
    % schz
    schz.m_eff_r(i) = median(schz.eff_r(:,i));
    schz.ci_eff_r(i,:) = abs(prctile(schz.eff_r(:,i),[prct_l prct_h])-schz.m_eff_r(i));
    schz.m_trans_r(i) = median(schz.trans_r(:,i));
    schz.ci_trans_r(i,:) = abs(prctile(schz.trans_r(:,i),[prct_l prct_h])-schz.m_trans_r(i));
    schz.m_dcmp_n_r(i) = median(schz.dcmp_n_r(:,i));
    schz.ci_dcmp_n_r(i,:) = abs(prctile(schz.dcmp_n_r(:,i),[prct_l prct_h])-schz.m_dcmp_n_r(i));
    schz.m_dcmp_max_r(i) = median(schz.dcmp_max_r(:,i));
    schz.ci_dcmp_max_r(i,:) = abs(prctile(schz.dcmp_max_r(:,i),[prct_l prct_h])-schz.m_dcmp_max_r(i));
    
    % ranksum across densities r-thr ctrl VS schz
    p_eff_r(i) = ranksum(ctrl.eff_r(:,i),schz.eff_r(:,i));
    p_trans_r(i) = ranksum(ctrl.trans_r(:,i),schz.trans_r(:,i));
    p_dcmp_n_r(i) = ranksum(ctrl.dcmp_n_r(:,i),schz.dcmp_n_r(:,i));
    p_dcmp_max_r(i) = ranksum(ctrl.dcmp_max_r(:,i),schz.dcmp_max_r(:,i));
    
    % P sub 1
    % ctrl
    ctrl.m_eff_r_psub1(i) = median(ctrl.eff_r(ctrl.psub1_id{i},i));
    ctrl.ci_eff_r_psub1(i,:) = abs(prctile(ctrl.eff_r(ctrl.psub1_id{i},i),[prct_l prct_h])-ctrl.m_eff_r_psub1(i));
    ctrl.m_trans_r_psub1(i) = median(ctrl.trans_r(ctrl.psub1_id{i},i));
    ctrl.ci_trans_r_psub1(i,:) = abs(prctile(ctrl.trans_r(ctrl.psub1_id{i},i),[prct_l prct_h])-ctrl.m_trans_r_psub1(i));
    ctrl.m_dcmp_n_r_psub1(i) = median(ctrl.dcmp_n_r(ctrl.psub1_id{i},i));
    ctrl.ci_dcmp_n_r_psub1(i,:) = abs(prctile(ctrl.dcmp_n_r(ctrl.psub1_id{i},i),[prct_l prct_h])-ctrl.m_dcmp_n_r_psub1(i));
    ctrl.m_dcmp_max_r_psub1(i) = median(ctrl.dcmp_max_r(ctrl.psub1_id{i},i));
    ctrl.ci_dcmp_max_r_psub1(i,:) = abs(prctile(ctrl.dcmp_max_r(ctrl.psub1_id{i},i),[prct_l prct_h])-ctrl.m_dcmp_max_r_psub1(i));

    % schz
    schz.m_eff_r_psub1(i) = median(schz.eff_r(schz.psub1_id{i},i));
    schz.ci_eff_r_psub1(i,:) = abs(prctile(schz.eff_r(schz.psub1_id{i},i),[prct_l prct_h])-schz.m_eff_r_psub1(i));
    schz.m_trans_r_psub1(i) = median(schz.trans_r(schz.psub1_id{i},i));
    schz.ci_trans_r_psub1(i,:) = abs(prctile(schz.trans_r(schz.psub1_id{i},i),[prct_l prct_h])-schz.m_trans_r_psub1(i));
    schz.m_dcmp_n_r_psub1(i) = median(schz.dcmp_n_r(schz.psub1_id{i},i));
    schz.ci_dcmp_n_r_psub1(i,:) = abs(prctile(schz.dcmp_n_r(schz.psub1_id{i},i),[prct_l prct_h])-schz.m_dcmp_n_r_psub1(i));
    schz.m_dcmp_max_r_psub1(i) = median(schz.dcmp_max_r(schz.psub1_id{i},i));
    schz.ci_dcmp_max_r_psub1(i,:) = abs(prctile(schz.dcmp_max_r(schz.psub1_id{i},i),[prct_l prct_h])-schz.m_dcmp_max_r_psub1(i));
    
    % ranksum across densities r-thr ctrl VS schz
    p_eff_r_psub1(i) = ranksum(ctrl.eff_r(ctrl.psub1_id{i},i),schz.eff_r(schz.psub1_id{i},i));
    p_trans_r_psub1(i) = ranksum(ctrl.trans_r(ctrl.psub1_id{i},i),schz.trans_r(schz.psub1_id{i},i));
    p_dcmp_n_r_psub1(i) = ranksum(ctrl.dcmp_n_r(ctrl.psub1_id{i},i),schz.dcmp_n_r(schz.psub1_id{i},i));
    p_dcmp_max_r_psub1(i) = ranksum(ctrl.dcmp_max_r(ctrl.psub1_id{i},i),schz.dcmp_max_r(schz.psub1_id{i},i));
    
end

ctrl.ci_eff_r = abs(ctrl.ci_eff_r); ctrl.ci_trans_r = abs(ctrl.ci_trans_r);
schz.ci_eff_r = abs(schz.ci_eff_r); schz.ci_trans_r = abs(schz.ci_trans_r);

ctrl.ci_dcmp_n_r = abs(ctrl.ci_dcmp_n_r); ctrl.ci_dcmp_max_r = abs(ctrl.ci_dcmp_max_r);
schz.ci_dcmp_max_r = abs(schz.ci_dcmp_max_r); schz.ci_dcmp_max_r = abs(schz.ci_dcmp_max_r);

ctrl.ci_eff_r_psub1 = abs(ctrl.ci_eff_r_psub1); ctrl.ci_trans_r_psub1 = abs(ctrl.ci_trans_r_psub1);
schz.ci_eff_r_psub1 = abs(schz.ci_eff_r_psub1); schz.ci_trans_r_psub1 = abs(schz.ci_trans_r_psub1);

ctrl.ci_dcmp_n_r_psub1 = abs(ctrl.ci_dcmp_n_r_psub1); ctrl.ci_dcmp_max_r_psub1 = abs(ctrl.ci_dcmp_max_r_psub1);
schz.ci_dcmp_max_r_psub1 = abs(schz.ci_dcmp_max_r_psub1); schz.ci_dcmp_max_r_psub1 = abs(schz.ci_dcmp_max_r_psub1);

for i = 1:1:ndens

    % signrank across densities r-thr VS p-thr - PSUB1
    % ctrl
    ppair_eff_rp_ctrl_psub1(i) = signrank(ctrl.eff_r(ctrl.psub1_id{i},i),ctrl.eff(ctrl.psub1_id{i},i));
    ppair_trans_rp_ctrl_psub1(i) = signrank(ctrl.trans_r(ctrl.psub1_id{i},i),ctrl.trans(ctrl.psub1_id{i},i));
    ppair_dcmp_n_rp_ctrl_psub1(i) = signrank(ctrl.dcmp_n_r(ctrl.psub1_id{i},i),ctrl.dcmp_n(ctrl.psub1_id{i},i));
    ppair_dcmp_max_rp_ctrl_psub1(i) = signrank(ctrl.dcmp_max_r(ctrl.psub1_id{i},i),ctrl.dcmp_max(ctrl.psub1_id{i},i));
    
    % schz
    ppair_eff_rp_schz_psub1(i) = signrank(schz.eff_r(schz.psub1_id{i},i),schz.eff(schz.psub1_id{i},i));
    ppair_trans_rp_schz_psub1(i) = signrank(schz.trans_r(schz.psub1_id{i},i),schz.trans(schz.psub1_id{i},i));
    ppair_dcmp_n_rp_schz_psub1(i) = signrank(schz.dcmp_n_r(schz.psub1_id{i},i),schz.dcmp_n(schz.psub1_id{i},i));
    ppair_dcmp_max_rp_schz_psub1(i) = signrank(schz.dcmp_max_r(schz.psub1_id{i},i),schz.dcmp_max(schz.psub1_id{i},i));
    
end

% within participant difference between r-thr and p-thr (subjects with P<1)
% ctrl
ctrl.d_eff_psub1 = nan(nc,ndens);
ctrl.d_trans_psub1 = nan(nc,ndens);
ctrl.d_dcmp_n_psub1 = nan(nc,ndens);
for i = 1:1:ndens
    ctrl.d_eff_psub1(ctrl.psub1_id{i},i) = ctrl.eff(ctrl.psub1_id{i},i)-ctrl.eff_r(ctrl.psub1_id{i},i);
    ctrl.d_trans_psub1(ctrl.psub1_id{i},i) = ctrl.trans(ctrl.psub1_id{i},i)-ctrl.trans_r(ctrl.psub1_id{i},i);
    ctrl.d_dcmp_n_psub1(ctrl.psub1_id{i},i) = ctrl.dcmp_n(ctrl.psub1_id{i},i)-ctrl.dcmp_n_r(ctrl.psub1_id{i},i);
    
    ctrl.m_d_eff_psub1(i) = median(ctrl.d_eff_psub1(ctrl.psub1_id{i},i));
    ctrl.m_d_trans_psub1(i) = median(ctrl.d_trans_psub1(ctrl.psub1_id{i},i));
    ctrl.m_d_dcmp_n_psub1(i) = median(ctrl.d_dcmp_n_psub1(ctrl.psub1_id{i},i));
    
    ctrl.ci_d_eff_psub1(i,:) = abs(prctile(ctrl.d_eff_psub1(ctrl.psub1_id{i},i),[prct_l prct_h])-ctrl.m_d_eff_psub1(i));
    ctrl.ci_d_trans_psub1(i,:) = abs(prctile(ctrl.d_trans_psub1(ctrl.psub1_id{i},i),[prct_l prct_h])-ctrl.m_d_trans_psub1(i));
    ctrl.ci_d_dcmp_n_psub1(i,:) = abs(prctile(ctrl.d_dcmp_n_psub1(ctrl.psub1_id{i},i),[prct_l prct_h])-ctrl.m_d_dcmp_n_psub1(i));
end

% schz
schz.d_eff_psub1 = nan(nc,ndens);
schz.d_trans_psub1 = nan(nc,ndens);
schz.d_dcmp_n_psub1 = nan(nc,ndens);
for i = 1:1:ndens
    schz.d_eff_psub1(schz.psub1_id{i},i) = schz.eff(schz.psub1_id{i},i)-schz.eff_r(schz.psub1_id{i},i);
    schz.d_trans_psub1(schz.psub1_id{i},i) = schz.trans(schz.psub1_id{i},i)-schz.trans_r(schz.psub1_id{i},i);
    schz.d_dcmp_n_psub1(schz.psub1_id{i},i) = schz.dcmp_n(schz.psub1_id{i},i)-schz.dcmp_n_r(schz.psub1_id{i},i);
    
    schz.m_d_eff_psub1(i) = median(schz.d_eff_psub1(schz.psub1_id{i},i));
    schz.m_d_trans_psub1(i) = median(schz.d_trans_psub1(schz.psub1_id{i},i));
    schz.m_d_dcmp_n_psub1(i) = median(schz.d_dcmp_n_psub1(schz.psub1_id{i},i));
    
    schz.ci_d_eff_psub1(i,:) = abs(prctile(schz.d_eff_psub1(schz.psub1_id{i},i),[prct_l prct_h])-schz.m_d_eff_psub1(i));
    schz.ci_d_trans_psub1(i,:) = abs(prctile(schz.d_trans_psub1(schz.psub1_id{i},i),[prct_l prct_h])-schz.m_d_trans_psub1(i));
    schz.ci_d_dcmp_n_psub1(i,:) = abs(prctile(schz.d_dcmp_n_psub1(schz.psub1_id{i},i),[prct_l prct_h])-schz.m_d_dcmp_n_psub1(i));
end

%% r VS P (P < 1)

fsize = 30;
fsize2 = 20;

% ctrl (r VS P)

% % efficiency
% boundedline_plot(1:ndens,ctrl.m_eff_r_psub1,ctrl.ci_eff_r_psub1,ctrl.m_eff_psub1,ctrl.ci_eff_psub1,...
%     ctrl.col_r,ctrl.col,alph,NaN,ppair_eff_rp_ctrl_psub1,'edge density \rho (%)','efficiency','ctrl r','ctrl p',[0 ndens+1],[0 0.8],fsize,fsize2,...
%     plot_path,['eff_ctrl_r_vs_p_psub1_boundedline_sc' sc '_alpha' num2str(100*alpha)]);
% 
% % transitivity
% boundedline_plot(1:ndens,ctrl.m_trans_r_psub1,ctrl.ci_trans_r_psub1,ctrl.m_trans_psub1,ctrl.ci_trans_psub1,...
%     ctrl.col_r,ctrl.col,alph,NaN,ppair_trans_rp_ctrl_psub1,'edge density \rho (%)','transitivity','ctrl r','ctrl p',[0 ndens+1],[0.3 0.9],fsize,fsize2,...
%     plot_path,['trans_ctrl_r_vs_p_psub1_boundedline_sc' sc '_alpha' num2str(100*alpha)]);

%%% difference

% efficiency
boundedline_plot_one(1:ndens,ctrl.m_d_eff_psub1,ctrl.ci_d_eff_psub1,ctrl.col,alph,ppair_eff_rp_ctrl_psub1,'edge density \rho (%)','\Delta eff. (P_{thr}-r_{thr})',[0 ndens+1],[-0.005 0.005],fsize,fsize2,...
    plot_path,['d_eff_ctrl_r_vs_p_psub1_boundedline_sc' sc '_alpha' num2str(100*alpha)]);

% transitivity
boundedline_plot_one(1:ndens,ctrl.m_d_trans_psub1,ctrl.ci_d_trans_psub1,ctrl.col,alph,ppair_trans_rp_ctrl_psub1,'edge density \rho (%)','\Delta trans. (P_{thr}-r_{thr})',[0 ndens+1],[-0.005 0.005],fsize,fsize2,...
    plot_path,['d_trans_ctrl_r_vs_p_psub1_boundedline_sc' sc '_alpha' num2str(100*alpha)]);

% schz (r VS P)

% % efficiency
% boundedline_plot(1:ndens,schz.m_eff_r_psub1,schz.ci_eff_r_psub1,schz.m_eff_psub1,schz.ci_eff_psub1,...
%     schz.col_r,schz.col,alph,NaN,ppair_eff_rp_schz_psub1,'edge density \rho (%)','efficiency','schz r','schz p',[0 ndens+1],[0 0.8],fsize,fsize2,...
%     plot_path,['eff_schz_r_vs_p_psub1_boundedline_sc' sc '_alpha' num2str(100*alpha)]);
% 
% % transitivity
% boundedline_plot(1:ndens,schz.m_trans_r_psub1,schz.ci_trans_r_psub1,schz.m_trans_psub1,schz.ci_trans_psub1,...
%     schz.col_r,schz.col,alph,NaN,ppair_trans_rp_schz_psub1,'edge density \rho (%)','transitivity','schz r','schz p',[0 ndens+1],[0.3 0.9],fsize,fsize2,...
%     plot_path,['trans_schz_r_vs_p_psub1_boundedline_sc' sc '_alpha' num2str(100*alpha)]);

%%% difference

% efficiency
boundedline_plot_one(1:ndens,schz.m_d_eff_psub1,schz.ci_d_eff_psub1,schz.col,alph,ppair_eff_rp_schz_psub1,'edge density \rho (%)','\Delta eff. (P_{thr}-r_{thr})',[0 ndens+1],[-0.005 0.005],fsize,fsize2,...
    plot_path,['d_eff_schz_r_vs_p_psub1_boundedline_sc' sc '_alpha' num2str(100*alpha)]);

% transitivity
boundedline_plot_one(1:ndens,schz.m_d_trans_psub1,schz.ci_d_trans_psub1,schz.col,alph,ppair_trans_rp_schz_psub1,'edge density \rho (%)','\Delta trans. (P_{thr}-r_{thr})',[0 ndens+1],[-0.005 0.005],fsize,fsize2,...
    plot_path,['d_trans_schz_r_vs_p_psub1_boundedline_sc' sc '_alpha' num2str(100*alpha)]);

%% compare R and P thresholded connectomes (cut-off when P=1 reached)

% number of edges that are the same
ctrl.rp_overlap = nan(ndens,nc);
schz.rp_overlap = nan(ndens,np);
for d = 1:1:ndens
    d

    % ctrl
    for i = ctrl.psub1_id{d}
        ctrl.rp_overlap(d,i) = length(intersect(intersect(find(ctrl.mat_dens_r(:,:,d,i)),triu_ind),intersect(find(ctrl.mat_dens(:,:,d,i)),triu_ind)))/length(intersect(find(ctrl.mat_dens(:,:,d,i)),triu_ind));
    end
    
    % schz
    for i = schz.psub1_id{d}
        schz.rp_overlap(d,i) = length(intersect(intersect(find(schz.mat_dens_r(:,:,d,i)),triu_ind),intersect(find(schz.mat_dens(:,:,d,i)),triu_ind)))/length(intersect(find(schz.mat_dens(:,:,d,i)),triu_ind));
    end
    
end

fsize = 30;
fsize2 = 20;

% plots

line_alpha = 0.5;

%%% separate plot for each group

%%% controls
figure;
hold on
for i = 1:1:nc; hc = plot(1:ndens,100*(1-ctrl.rp_overlap(:,i)),'Color',ctrl.col_light,'LineWidth',lwd_ind); end
hca = plot(1:ndens,100*(1-nanmean(ctrl.rp_overlap,2)),'Color',ctrl.col*col_ind,'LineWidth',lwd);
hold off

set(gca,'FontSize',fsize2); box off
xlabel('edge density \rho (%)','FontSize',fsize,'FontName','Arial');
ylabel('\Delta edges (%)','FontSize',fsize,'FontName','Arial');
xlim([0 ndens+1]); ylim([-1,15]); set(gca,'YTick',[0:5:15])

pos = get(gcf,'Position');
set(gcf,'Position', [pos(1) pos(2) 1.25*pos(3) 0.70*pos(4)]); set(gcf,'color','w');

print([plot_path 'cobre_rp_diff_ctrl.png'],'-dpng')

%%% schizophrenia
figure;
hold on
for i = 1:1:np; hp = plot(1:ndens,100*(1-schz.rp_overlap(:,i)),'Color',schz.col_light,'LineWidth',lwd_ind); end
hpa = plot(1:ndens,100*(1-nanmean(schz.rp_overlap,2)),'Color',schz.col*0.8,'LineWidth',lwd);
hold off

set(gca,'FontSize',fsize2); box off
xlabel('edge density \rho (%)','FontSize',fsize,'FontName','Arial');
ylabel('\Delta edges (%)','FontSize',fsize,'FontName','Arial');
xlim([0 ndens+1]); ylim([-1,15]); set(gca,'YTick',[0:5:15])

pos = get(gcf,'Position');
set(gcf,'Position', [pos(1) pos(2) 1.25*pos(3) 0.70*pos(4)]); set(gcf,'color','w');

print([plot_path 'cobre_rp_diff_schz.png'],'-dpng')

%% r vs P

fsize = 24;
fsize2 = 16;

figure; hold on;

ex_id = 8;
temp_r = ctrl.r(:,:,ex_id);
temp_p = ctrl.fdr(:,:,ex_id);
scatter(temp_r(triu_ind),temp_p(triu_ind),cex,ctrl.col,'filled','o');

ex_id = 2;
temp_r = schz.r(:,:,ex_id);
temp_p = schz.fdr(:,:,ex_id);
scatter(temp_r(triu_ind),temp_p(triu_ind),cex,schz.col,'filled','o');

set(gca,'FontSize',fsize2); %box off
xlabel('Pearson r','FontSize',fsize,'FontName','Arial');
ylabel('FDR-adjusted P value','FontSize',fsize,'FontName','Arial');

%pos = get(gcf,'Position');
%set(gcf,'Position', [pos(1) pos(2) 1.25*pos(3) 0.70*pos(4)]); 
set(gcf,'color','w');

xlim([-0.05 1.05])
ylim([-0.05 1.05])

% inset limit box
pi = -0.01; pf = 0.11;
ri = 0.45; rf = 0.65;

plot([ri rf],[pi pi],'k','LineWidth',2) % horiz bottom
plot([ri rf],[pf pf],'k','LineWidth',2) % horiz top
plot([ri ri],[pi pf],'k','LineWidth',2) % vert left
plot([rf rf],[pi pf],'k','LineWidth',2) % vert right

hold off;

export_fig([plot_path 'rp_rel_full.pdf'],'-pdf','-nocrop')

% inset
figure; hold on;
ex_id = 8;
temp_r = ctrl.r(:,:,ex_id);
temp_p = ctrl.fdr(:,:,ex_id);
scatter(temp_r(triu_ind),temp_p(triu_ind),cex,ctrl.col_light,'filled','o');
ex_id = 2;
temp_r = schz.r(:,:,ex_id);
temp_p = schz.fdr(:,:,ex_id);
scatter(temp_r(triu_ind),temp_p(triu_ind),cex,schz.col_light,'filled','o');
xlim([ri rf])
ylim([pi pf])

set(gca,'FontSize',fsize2); %box off
set(gca,'YTick',[0 0.05 0.1])
set(gca,'XTick',[0.5 0.6])
set(gcf,'color','w');

export_fig([plot_path 'rp_rel_inset.pdf'],'-pdf','-nocrop')

%% consistency of edges (P<1 only, for both P and r thresholds)

% within-group consistency (how many times is each edge present across the
% group, at each density and for each thresholding method?)
for d = 1:1:ndens

    % ctrl
    ctrl.r_overlap(:,:,d) = sum(squeeze(logical(ctrl.mat_dens_r(:,:,d,ctrl.psub1_id{d}))),3);
    ctrl.p_overlap(:,:,d) = sum(squeeze(logical(ctrl.mat_dens(:,:,d,ctrl.psub1_id{d}))),3);
    % schz
    schz.r_overlap(:,:,d) = sum(squeeze(logical(schz.mat_dens_r(:,:,d,schz.psub1_id{d}))),3);
    schz.p_overlap(:,:,d) = sum(squeeze(logical(schz.mat_dens(:,:,d,schz.psub1_id{d}))),3);

end

% frequencies of overlap - count the number of edges at each level of consistency
comb_tbl_r = zeros(ndens,nc+np+1);
comb_tbl_p = zeros(ndens,nc+np+1);
for d = 1:1:ndens

    %%% ctrl
    % r
    temp_ctrl = ctrl.r_overlap(:,:,d);
    tbl1_r = tabulate(temp_ctrl(triu_ind));
    ctrl.tbl_r(d,tbl1_r(:,1)+1) = tbl1_r(:,2);
    % p
    temp_ctrl = ctrl.p_overlap(:,:,d);
    tbl1_p = tabulate(temp_ctrl(triu_ind)); %
    ctrl.tbl_p(d,tbl1_p(:,1)+1) = tbl1_p(:,2);
    %%% schz
    % r
    temp_ctrl = schz.r_overlap(:,:,d);
    tbl1_r = tabulate(temp_ctrl(triu_ind));
    schz.tbl_r(d,tbl1_r(:,1)+1) = tbl1_r(:,2);
    % p
    temp_ctrl = schz.p_overlap(:,:,d);
    tbl1_p = tabulate(temp_ctrl(triu_ind));
    schz.tbl_p(d,tbl1_p(:,1)+1) = tbl1_p(:,2);

end

% difference in frequencies
ctrl.tbl_d = ctrl.tbl_p-ctrl.tbl_r; 
schz.tbl_d = schz.tbl_p-schz.tbl_r; 

% set RHS of plots (corresponding to subjects with P=1) to NaN
for d = 1:1:ndens
    
    % r
    ctrl.tbl_r(d,(length(ctrl.psub1_id{d})+2):end) = NaN;
    schz.tbl_r(d,(length(schz.psub1_id{d})+2):end) = NaN;

    % p
    ctrl.tbl_p(d,(length(ctrl.psub1_id{d})+2):end) = NaN;
    schz.tbl_p(d,(length(schz.psub1_id{d})+2):end) = NaN;

    % diff
    ctrl.tbl_d(d,(length(ctrl.psub1_id{d})+2):end) = NaN;
    schz.tbl_d(d,(length(schz.psub1_id{d})+2):end) = NaN;

end

%%%
% % frequency plots by group and method (not plotted in manuscript)
% % this is an intermediate step to the *difference of frequencies*
% % ctrl - r
% temp = log10(ctrl.tbl_r); temp(isnan(temp)) = max(log10(ctrl.tbl_r(:)))+0.1;
% figure; imagesc(temp); axis xy; c = colorbar; cmap = colormap(parula); 
% cmap(end,:) = rgb('whitesmoke'); colormap(cmap)
% hold on; for d = 1:1:ndens % add black line
%     plot([length(ctrl.psub1_id{d})+1.5 length(ctrl.psub1_id{d})+1.5],[d-0.5 d+0.5],'k','linewidth',1); % vertical segment - 1.5 extra in x-coord = 1 (count starts from 0) + 0.5 (imagesc alignment)
%     if d < 35; if length(ctrl.psub1_id{d}) > length(ctrl.psub1_id{d+1}); plot([length(ctrl.psub1_id{d})+1.5-(length(ctrl.psub1_id{d})-length(ctrl.psub1_id{d+1})) length(ctrl.psub1_id{d})+1.5],[d+0.5 d+0.5],'k','linewidth',1); end; end
% end; hold off
% set(gca,'FontSize',fsize2); set(gca,'XTick',[1:5:(nc+1)],'XTickLabel',[(1:5:(nc+1))-1])
% xlabel('consistency (# participants)','FontSize',fsize,'FontName','Arial');
% ylabel('edge density \rho (%)','FontSize',fsize,'FontName','Arial');
% ylabel(c,'frequency (# edges)','Rotation',-90,'VerticalAlignment','bottom','FontSize',fsize,'FontName','Arial'); % 
% set(gcf,'color','w');
% pos = get(gcf,'Position'); set(gcf,'Position', [pos(1) pos(2) 1.5*pos(3) pos(4)]);
% print([plot_path '/cobre_r_consistency_ctrl.png'],'-dpng')
% 
% % ctrl - p
% temp = log10(ctrl.tbl_p); temp(isnan(temp)) = max(log10(ctrl.tbl_p(:)))+0.1;
% figure; imagesc(temp); axis xy; c = colorbar; cmap = colormap(parula); 
% cmap(end,:) = rgb('whitesmoke'); colormap(cmap)
% hold on; for d = 1:1:ndens % add black line
%     plot([length(ctrl.psub1_id{d})+1.5 length(ctrl.psub1_id{d})+1.5],[d-0.5 d+0.5],'k','linewidth',1); % vertical segment - 1.5 extra in x-coord = 1 (count starts from 0) + 0.5 (imagesc alignment)
%     if d < 35; if length(ctrl.psub1_id{d}) > length(ctrl.psub1_id{d+1}); plot([length(ctrl.psub1_id{d})+1.5-(length(ctrl.psub1_id{d})-length(ctrl.psub1_id{d+1})) length(ctrl.psub1_id{d})+1.5],[d+0.5 d+0.5],'k','linewidth',1); end; end
% end; hold off
% set(gca,'FontSize',fsize2); set(gca,'XTick',[1:5:(nc+1)],'XTickLabel',[(1:5:(nc+1))-1])
% xlabel('consistency (# participants)','FontSize',fsize,'FontName','Arial');
% ylabel('edge density \rho (%)','FontSize',fsize,'FontName','Arial');
% ylabel(c,'frequency (# edges)','Rotation',-90,'VerticalAlignment','bottom','FontSize',fsize,'FontName','Arial'); % 
% set(gcf,'color','w');
% pos = get(gcf,'Position'); set(gcf,'Position', [pos(1) pos(2) 1.5*pos(3) pos(4)]);
% print([plot_path '/cobre_p_consistency_ctrl.png'],'-dpng')
% 
% % schz - r
% temp = log10(schz.tbl_r); temp(isnan(temp)) = max(log10(schz.tbl_r(:)))+0.1;
% figure; imagesc(temp); axis xy; c = colorbar; cmap = colormap(parula); 
% cmap(end,:) = rgb('whitesmoke'); colormap(cmap)
% hold on; for d = 1:1:ndens % add black line
%     plot([length(schz.psub1_id{d})+1.5 length(schz.psub1_id{d})+1.5],[d-0.5 d+0.5],'k','linewidth',1); % vertical segment - 1.5 extra in x-coord = 1 (count starts from 0) + 0.5 (imagesc alignment)
%     if d < 35; if length(schz.psub1_id{d}) > length(schz.psub1_id{d+1}); plot([length(schz.psub1_id{d})+1.5-(length(schz.psub1_id{d})-length(schz.psub1_id{d+1})) length(schz.psub1_id{d})+1.5],[d+0.5 d+0.5],'k','linewidth',1); end; end
% end; hold off
% set(gca,'FontSize',fsize2); set(gca,'XTick',[1:5:(nc+1)],'XTickLabel',[(1:5:(nc+1))-1])
% xlabel('consistency (# participants)','FontSize',fsize,'FontName','Arial');
% ylabel('edge density \rho (%)','FontSize',fsize,'FontName','Arial');
% ylabel(c,'frequency (# edges)','Rotation',-90,'VerticalAlignment','bottom','FontSize',fsize,'FontName','Arial'); % 
% set(gcf,'color','w');
% pos = get(gcf,'Position'); set(gcf,'Position', [pos(1) pos(2) 1.5*pos(3) pos(4)]);
% print([plot_path '/cobre_r_consistency_schz.png'],'-dpng')
% 
% % schz - p
% temp = log10(schz.tbl_p); temp(isnan(temp)) = max(log10(schz.tbl_p(:)))+0.1;
% figure; imagesc(temp); axis xy; c = colorbar; cmap = colormap(parula); 
% cmap(end,:) = rgb('whitesmoke'); colormap(cmap)
% hold on; for d = 1:1:ndens % add black line
%     plot([length(schz.psub1_id{d})+1.5 length(schz.psub1_id{d})+1.5],[d-0.5 d+0.5],'k','linewidth',1); % vertical segment - 1.5 extra in x-coord = 1 (count starts from 0) + 0.5 (imagesc alignment)
%     if d < 35; if length(schz.psub1_id{d}) > length(schz.psub1_id{d+1}); plot([length(schz.psub1_id{d})+1.5-(length(schz.psub1_id{d})-length(schz.psub1_id{d+1})) length(schz.psub1_id{d})+1.5],[d+0.5 d+0.5],'k','linewidth',1); end; end
% end; hold off
% set(gca,'FontSize',fsize2); set(gca,'XTick',[1:5:(nc+1)],'XTickLabel',[(1:5:(nc+1))-1])
% xlabel('consistency (# participants)','FontSize',fsize,'FontName','Arial');
% ylabel('edge density \rho (%)','FontSize',fsize,'FontName','Arial');
% ylabel(c,'frequency (# edges)','Rotation',-90,'VerticalAlignment','bottom','FontSize',fsize,'FontName','Arial'); % 
% set(gcf,'color','w');
% pos = get(gcf,'Position'); set(gcf,'Position', [pos(1) pos(2) 1.5*pos(3) pos(4)]);
% print([plot_path '/cobre_p_consistency_schz.png'],'-dpng')
%%%

% difference plots (with permutation P value markers, obtained using code "pr_thr_consistency_perm)
load('pr_thr_consistency_perm_psub1.mat')
nperm = size(ctrl_tbl_d_perm,3);

% compare empirical difference to permuted differences - ctrl
for d = 1:1:ndens 
    for i = 1:1:(length(ctrl.psub1_id{d})+1)
        if ctrl.tbl_d(d,i) >= 0
            ctrl.tbl_d_p(d,i) = sum(ctrl_tbl_d_perm(d,i,:)>ctrl.tbl_d(d,i))/nperm;
        elseif ctrl.tbl_d(d,i) < 0
            ctrl.tbl_d_p(d,i) = sum(ctrl_tbl_d_perm(d,i,:)<ctrl.tbl_d(d,i))/nperm;
        end
    end
end

% plot - ctrl
temp = ctrl.tbl_d; temp(isnan(temp)) = min(-max(abs(ctrl.tbl_d(:)))-1);
figure; imagesc(temp); axis xy; c = colorbar; cmap = colormap(redblue); 
caxis([-max(abs(ctrl.tbl_d(:))),max(abs(ctrl.tbl_d(:)))]); cmap(1,:) = rgb('whitesmoke'); colormap(cmap)
hold on; for d = 1:1:ndens % add black line
    plot([length(ctrl.psub1_id{d})+1.5 length(ctrl.psub1_id{d})+1.5],[d-0.5 d+0.5],'k','linewidth',1); % vertical segment - 1.5 extra in x-coord = 1 (count starts from 0) + 0.5 (imagesc alignment)
    if d < 35; if length(ctrl.psub1_id{d}) > length(ctrl.psub1_id{d+1}); plot([length(ctrl.psub1_id{d})+1.5-(length(ctrl.psub1_id{d})-length(ctrl.psub1_id{d+1})) length(ctrl.psub1_id{d})+1.5],[d+0.5 d+0.5],'k','linewidth',1); end; end
end
%%% P value markers
for d = 1:1:ndens
    for i = 1:1:(length(ctrl.psub1_id{d})+1)
        if ctrl.tbl_d_p(d,i) < 0.01
            scatter(i,d,5,'filled','MarkerFaceColor','black','MarkerEdgeColor','black')
        end
    end
end
%%%
hold off
set(gca,'FontSize',fsize2); set(gca,'XTick',[1:5:(nc+1)],'XTickLabel',[(1:5:(nc+1))-1])
xlabel('consistency (# participants)','FontSize',fsize,'FontName','Arial');
ylabel('edge density \rho (%)','FontSize',fsize,'FontName','Arial');
ylabel(c,'\Delta frequency (# edges)','Rotation',-90,'VerticalAlignment','bottom','FontSize',fsize,'FontName','Arial'); % 
set(gcf,'color','w');
pos = get(gcf,'Position'); set(gcf,'Position', [pos(1) pos(2) 1.5*pos(3) pos(4)]);
print([plot_path '/cobre_rp_d_consistency_ctrl_pmarkers.png'],'-dpng')

% compare empirical difference to permuted differences - schz
for d = 1:1:ndens
    for i = 1:1:(length(schz.psub1_id{d})+1)
        if schz.tbl_d(d,i) >= 0
            schz.tbl_d_p(d,i) = sum(schz_tbl_d_perm(d,i,:)>schz.tbl_d(d,i))/nperm;
        elseif schz.tbl_d(d,i) < 0
            schz.tbl_d_p(d,i) = sum(schz_tbl_d_perm(d,i,:)<schz.tbl_d(d,i))/nperm;
        end
    end
end

% plot - schz
temp = schz.tbl_d; temp(isnan(temp)) = min(-max(abs(schz.tbl_d(:)))-1);
figure; imagesc(temp); axis xy; c = colorbar; cmap = colormap(redblue); 
caxis([-max(abs(schz.tbl_d(:))),max(abs(schz.tbl_d(:)))]); cmap(1,:) = rgb('whitesmoke'); colormap(cmap)
hold on; for d = 1:1:ndens % add black line
    plot([length(schz.psub1_id{d})+1.5 length(schz.psub1_id{d})+1.5],[d-0.5 d+0.5],'k','linewidth',1); % vertical segment - 1.5 extra in x-coord = 1 (count starts from 0) + 0.5 (imagesc alignment)
    if d < 35; if length(schz.psub1_id{d}) > length(schz.psub1_id{d+1}); plot([length(schz.psub1_id{d})+1.5-(length(schz.psub1_id{d})-length(schz.psub1_id{d+1})) length(schz.psub1_id{d})+1.5],[d+0.5 d+0.5],'k','linewidth',1); end; end
end
%%% P value markers
for d = 1:1:ndens
    for i = 1:1:(length(schz.psub1_id{d})+1)
        if schz.tbl_d_p(d,i) < 0.01
            scatter(i,d,5,'filled','MarkerFaceColor','black','MarkerEdgeColor','black')
        end
    end
end
%%%
hold off
set(gca,'FontSize',fsize2); set(gca,'XTick',[1:5:(nc+1)],'XTickLabel',[(1:5:(nc+1))-1])
xlabel('consistency (# participants)','FontSize',fsize,'FontName','Arial');
ylabel('edge density \rho (%)','FontSize',fsize,'FontName','Arial');
ylabel(c,'\Delta frequency (# edges)','Rotation',-90,'VerticalAlignment','bottom','FontSize',fsize,'FontName','Arial'); % 
set(gcf,'color','w');
pos = get(gcf,'Position'); set(gcf,'Position', [pos(1) pos(2) 1.5*pos(3) pos(4)]);
print([plot_path '/cobre_rp_d_consistency_schz_pmarkers.png'],'-dpng')

%%% first column as column, with markers (overlap = 0, across densities)
% ctrl
cmap = colormap(redblue); close(gcf); x = ctrl.tbl_d(:,1); %data to be plotted
min_val = -max(abs(ctrl.tbl_d(:))); max_val = max(abs(ctrl.tbl_d(:))); y = floor(((x-min_val)/(max_val-min_val))*63)+1; 
figure; hold on; for i = 1:ndens; h = barh(i,ctrl.tbl_d(i,1)); set(h,'FaceColor',cmap(y(i),:)); end
scatter(zeros(1,sum(ctrl.tbl_d_p(:,1)<0.01)),find(ctrl.tbl_d_p(:,1)<0.01),cex,rgb('black'),'filled','o'); hold off;
%ax = gca; ax.YAxisLocation = 'Right'
ylim([0 ndens+1]); xlim([-100 350]);
set(gca,'FontSize',fsize2); %set(gca,'XTick',[1:5:(nc+1)],'XTickLabel',[(1:5:(nc+1))-1])
ylabel('edge density \rho (%)','FontSize',fsize,'FontName','Arial');
xlabel('\Delta freq. (# edges)','FontSize',fsize,'FontName','Arial');
set(gcf,'color','w');
pos = get(gcf,'Position'); set(gcf,'Position', [pos(1) pos(2) 0.5*pos(3) pos(4)]);
print([plot_path '/cobre_rp_d_consistency_ctrl_edge0_hor.png'],'-dpng')

% schz
cmap = colormap(redblue); close(gcf); x = schz.tbl_d(:,1); %data to be plotted
min_val = -max(abs(schz.tbl_d(:))); max_val = max(abs(schz.tbl_d(:))); y = floor(((x-min_val)/(max_val-min_val))*63)+1; 
figure; hold on; for i = 1:ndens; h = barh(i,schz.tbl_d(i,1)); set(h,'FaceColor',cmap(y(i),:)); end
scatter(zeros(1,sum(schz.tbl_d_p(:,1)<0.01)),find(schz.tbl_d_p(:,1)<0.01),cex,rgb('black'),'filled','o'); hold off;
%ax = gca; ax.YAxisLocation = 'Right'
ylim([0 ndens+1]); xlim([-50 550]);
set(gca,'FontSize',fsize2); %set(gca,'XTick',[1:5:(nc+1)],'XTickLabel',[(1:5:(nc+1))-1])
ylabel('edge density \rho (%)','FontSize',fsize,'FontName','Arial');
xlabel('\Delta freq. (# edges)','FontSize',fsize,'FontName','Arial');
set(gcf,'color','w');
pos = get(gcf,'Position'); set(gcf,'Position', [pos(1) pos(2) 0.5*pos(3) pos(4)]);
print([plot_path '/cobre_rp_d_consistency_schz_edge0_hor.png'],'-dpng')

%%  effect sizes for all participants

%%% controls have exception for density = 1, as all are significant at that density

for i = 1:1:ndens
    i
    % efficiency
    [p_eff_psub1(i),r_eff_psub1(i)] =           ranksum_effect_size(ctrl.eff(ctrl.psub1_id{i},i),schz.eff(schz.psub1_id{i},i));
    [p_eff_sig(i),r_eff_sig(i)] =               ranksum_effect_size(ctrl.eff(ctrl.sig_id{i},i),schz.eff(schz.sig_id{i},i));
    if i > 1
        [p_eff_ctrl_pint(i),r_eff_ctrl_pint(i)] =   ranksum_effect_size(ctrl.eff(ctrl.sig_id{i},i),ctrl.eff(setdiff(ctrl.psub1_id{i},ctrl.sig_id{i}),i));
    end
    [p_eff_schz_pint(i),r_eff_schz_pint(i)] =   ranksum_effect_size(schz.eff(schz.sig_id{i},i),schz.eff(setdiff(schz.psub1_id{i},schz.sig_id{i}),i));
    
    % efficiency normalised by random
    [p_eff_rand_psub1(i),r_eff_rand_psub1(i)] =         ranksum_effect_size(ctrl.eff_rand(ctrl.psub1_id{i},i),schz.eff_rand(schz.psub1_id{i},i));
    [p_eff_rand_sig(i),r_eff_rand_sig(i)] =             ranksum_effect_size(ctrl.eff_rand(ctrl.sig_id{i},i),schz.eff_rand(schz.sig_id{i},i)); % ctrl VS schz (sig only)
    if i > 1
        [p_eff_rand_ctrl_pint(i),r_eff_rand_ctrl_pint(i)] = ranksum_effect_size(ctrl.eff_rand(ctrl.sig_id{i},i),ctrl.eff_rand(setdiff(ctrl.psub1_id{i},ctrl.sig_id{i}),i));
    end
    [p_eff_rand_schz_pint(i),r_eff_rand_schz_pint(i)] = ranksum_effect_size(schz.eff_rand(schz.sig_id{i},i),schz.eff_rand(setdiff(schz.psub1_id{i},schz.sig_id{i}),i));
    
    % transitivity
    [p_trans_psub1(i),r_trans_psub1(i)] =           ranksum_effect_size(ctrl.trans(ctrl.psub1_id{i},i),schz.trans(schz.psub1_id{i},i));
    [p_trans_sig(i),r_trans_sig(i)] =               ranksum_effect_size(ctrl.trans(ctrl.sig_id{i},i),schz.trans(schz.sig_id{i},i));
    if i > 1
        [p_trans_ctrl_pint(i),r_trans_ctrl_pint(i)] =   ranksum_effect_size(ctrl.trans(ctrl.sig_id{i},i),ctrl.trans(setdiff(ctrl.psub1_id{i},ctrl.sig_id{i}),i));
    end
    [p_trans_schz_pint(i),r_trans_schz_pint(i)] =   ranksum_effect_size(schz.trans(schz.sig_id{i},i),schz.trans(setdiff(schz.psub1_id{i},schz.sig_id{i}),i));
    
    % transitivity normalised by random
    [p_trans_rand_psub1(i),r_trans_rand_psub1(i)] =         ranksum_effect_size(ctrl.trans_rand(ctrl.psub1_id{i},i),schz.trans_rand(schz.psub1_id{i},i));
    [p_trans_rand_sig(i),r_trans_rand_sig(i)] =             ranksum_effect_size(ctrl.trans_rand(ctrl.sig_id{i},i),schz.trans_rand(schz.sig_id{i},i)); % ctrl VS schz (sig only)
    if i > 1
        [p_trans_rand_ctrl_pint(i),r_trans_rand_ctrl_pint(i)] = ranksum_effect_size(ctrl.trans_rand(ctrl.sig_id{i},i),ctrl.trans_rand(setdiff(ctrl.psub1_id{i},ctrl.sig_id{i}),i));
    end
    [p_trans_rand_schz_pint(i),r_trans_rand_schz_pint(i)] = ranksum_effect_size(schz.trans_rand(schz.sig_id{i},i),schz.trans_rand(setdiff(schz.psub1_id{i},schz.sig_id{i}),i));
   
    % # components
    [p_dcmp_n_psub1(i),r_dcmp_n_psub1(i)] =         ranksum_effect_size(ctrl.dcmp_n(ctrl.psub1_id{i},i),schz.dcmp_n(schz.psub1_id{i},i));
    [p_dcmp_n_sig(i),r_dcmp_n_sig(i)] =             ranksum_effect_size(ctrl.dcmp_n(ctrl.sig_id{i},i),schz.dcmp_n(schz.sig_id{i},i));
    if i > 1
        [p_dcmp_n_ctrl_pint(i),r_dcmp_n_ctrl_pint(i)] = ranksum_effect_size(ctrl.dcmp_n(ctrl.sig_id{i},i),ctrl.dcmp_n(setdiff(ctrl.psub1_id{i},ctrl.sig_id{i}),i));
    end
    [p_dcmp_n_schz_pint(i),r_dcmp_n_schz_pint(i)] = ranksum_effect_size(schz.dcmp_n(schz.sig_id{i},i),schz.dcmp_n(setdiff(schz.psub1_id{i},schz.sig_id{i}),i));
    
    % mean correlation
    [p_avg_wei_psub1(i),r_avg_wei_psub1(i)] =           ranksum_effect_size(ctrl.avg_wei_full(ctrl.psub1_id{i}),schz.avg_wei_full(schz.psub1_id{i}));
    [p_avg_wei_sig(i),r_avg_wei_sig(i)] =               ranksum_effect_size(ctrl.avg_wei_full(ctrl.sig_id{i}),schz.avg_wei_full(schz.sig_id{i}));
    if i > 1
        [p_avg_wei_ctrl_pint(i),r_avg_wei_ctrl_pint(i)] =   ranksum_effect_size(ctrl.avg_wei_full(ctrl.sig_id{i}),ctrl.avg_wei_full(setdiff(ctrl.psub1_id{i},ctrl.sig_id{i})));
    end
    [p_avg_wei_schz_pint(i),r_avg_wei_schz_pint(i)] =   ranksum_effect_size(schz.avg_wei_full(schz.sig_id{i}),schz.avg_wei_full(setdiff(schz.psub1_id{i},schz.sig_id{i})));
    
    % efficiency - r VS p
    % ctrl
    [ppair_eff_rp_ctrl_psub1(i),rpair_eff_rp_ctrl_psub1(i)] = signrank_effect_size(ctrl.eff(ctrl.psub1_id{i},i),ctrl.eff_r(ctrl.psub1_id{i},i)); %
    [ppair_trans_rp_ctrl_psub1(i),rpair_trans_rp_ctrl_psub1(i)] = signrank_effect_size(ctrl.trans(ctrl.psub1_id{i},i),ctrl.trans_r(ctrl.psub1_id{i},i)); %
    
    % schz
    [ppair_eff_rp_schz_psub1(i),rpair_eff_rp_schz_psub1(i)] = signrank_effect_size(schz.eff(schz.psub1_id{i},i),schz.eff_r(schz.psub1_id{i},i)); % 
    [ppair_trans_rp_schz_psub1(i),rpair_trans_rp_schz_psub1(i)] = signrank_effect_size(schz.trans(schz.psub1_id{i},i),schz.trans_r(schz.psub1_id{i},i)); % 
    
    % efficiency - r VS p normalised by random
    % ctrl
    [ppair_eff_rand_rp_ctrl_psub1(i),rpair_eff_rand_rp_ctrl_psub1(i)] = signrank_effect_size(ctrl.eff_rand(ctrl.psub1_id{i},i),ctrl.eff_r_rand(ctrl.psub1_id{i},i)); %
    [ppair_trans_rand_rp_ctrl_psub1(i),rpair_trans_rand_rp_ctrl_psub1(i)] = signrank_effect_size(ctrl.trans_rand(ctrl.psub1_id{i},i),ctrl.trans_r_rand(ctrl.psub1_id{i},i)); %
    
    % schz
    [ppair_eff_rand_rp_schz_psub1(i),rpair_eff_rand_rp_schz_psub1(i)] = signrank_effect_size(schz.eff_rand(schz.psub1_id{i},i),schz.eff_r_rand(schz.psub1_id{i},i)); % 
    [ppair_trans_rand_rp_schz_psub1(i),rpair_trans_rand_rp_schz_psub1(i)] = signrank_effect_size(schz.trans_rand(schz.psub1_id{i},i),schz.trans_r_rand(schz.psub1_id{i},i)); % 
    
end

%% format effect sizes and P values

% number of significant figures
nrr = 2; % correlation
nrp = 2; % P value

for d = 1:1:ndens
    
    % efficiency
    eff_rp{d,1} = [num2str(roundsd(r_eff_psub1(d),nrr)) ' (' num2str(roundsd(p_eff_psub1(d),nrp)) ')'];
    eff_sig_rp{d,1} = [num2str(roundsd(r_eff_sig(d),nrr)) ' (' num2str(roundsd(p_eff_sig(d),nrp)) ')'];
    eff_ctrl_rp{d,1} = [num2str(roundsd(r_eff_ctrl_pint(d),nrr)) ' (' num2str(roundsd(p_eff_ctrl_pint(d),nrp)) ')'];
    eff_schz_rp{d,1} = [num2str(roundsd(r_eff_schz_pint(d),nrr)) ' (' num2str(roundsd(p_eff_schz_pint(d),nrp)) ')'];
    
    % efficiency normalised by random
    eff_rand_rp{d,1} = [num2str(roundsd(r_eff_rand_psub1(d),nrr)) ' (' num2str(roundsd(p_eff_rand_psub1(d),nrp)) ')'];
    eff_rand_sig_rp{d,1} = [num2str(roundsd(r_eff_rand_sig(d),nrr)) ' (' num2str(roundsd(p_eff_rand_sig(d),nrp)) ')'];
    eff_rand_ctrl_rp{d,1} = [num2str(roundsd(r_eff_rand_ctrl_pint(d),nrr)) ' (' num2str(roundsd(p_eff_rand_ctrl_pint(d),nrp)) ')'];
    eff_rand_schz_rp{d,1} = [num2str(roundsd(r_eff_rand_schz_pint(d),nrr)) ' (' num2str(roundsd(p_eff_rand_schz_pint(d),nrp)) ')'];
    
    % transitivity
    trans_rp{d,1} = [num2str(roundsd(r_trans_psub1(d),nrr)) ' (' num2str(roundsd(p_trans_psub1(d),nrp)) ')'];
    trans_sig_rp{d,1} = [num2str(roundsd(r_trans_sig(d),nrr)) ' (' num2str(roundsd(p_trans_sig(d),nrp)) ')'];
    trans_ctrl_rp{d,1} = [num2str(roundsd(r_trans_ctrl_pint(d),nrr)) ' (' num2str(roundsd(p_trans_ctrl_pint(d),nrp)) ')'];
    trans_schz_rp{d,1} = [num2str(roundsd(r_trans_schz_pint(d),nrr)) ' (' num2str(roundsd(p_trans_schz_pint(d),nrp)) ')'];
    
    % transitivity normalised by random
    trans_rand_rp{d,1} = [num2str(roundsd(r_trans_rand_psub1(d),nrr)) ' (' num2str(roundsd(p_trans_rand_psub1(d),nrp)) ')'];
    trans_rand_sig_rp{d,1} = [num2str(roundsd(r_trans_rand_sig(d),nrr)) ' (' num2str(roundsd(p_trans_rand_sig(d),nrp)) ')'];
    trans_rand_ctrl_rp{d,1} = [num2str(roundsd(r_trans_rand_ctrl_pint(d),nrr)) ' (' num2str(roundsd(p_trans_rand_ctrl_pint(d),nrp)) ')'];
    trans_rand_schz_rp{d,1} = [num2str(roundsd(r_trans_rand_schz_pint(d),nrr)) ' (' num2str(roundsd(p_trans_rand_schz_pint(d),nrp)) ')'];
    
    % # components
    dcmp_n_rp{d,1} = [num2str(roundsd(r_dcmp_n_psub1(d),nrr)) ' (' num2str(roundsd(p_dcmp_n_psub1(d),nrp)) ')'];
    dcmp_n_sig_rp{d,1} = [num2str(roundsd(r_dcmp_n_sig(d),nrr)) ' (' num2str(roundsd(p_dcmp_n_sig(d),nrp)) ')'];
    dcmp_n_ctrl_rp{d,1} = [num2str(roundsd(r_dcmp_n_ctrl_pint(d),nrr)) ' (' num2str(roundsd(p_dcmp_n_ctrl_pint(d),nrp)) ')'];
    dcmp_n_schz_rp{d,1} = [num2str(roundsd(r_dcmp_n_schz_pint(d),nrr)) ' (' num2str(roundsd(p_dcmp_n_schz_pint(d),nrp)) ')'];
    
    % mean correlation
    avg_wei_rp{d,1} = [num2str(roundsd(r_avg_wei_psub1(d),nrr)) ' (' num2str(roundsd(p_avg_wei_psub1(d),nrp)) ')'];
    avg_wei_sig_rp{d,1} = [num2str(roundsd(r_avg_wei_sig(d),nrr)) ' (' num2str(roundsd(p_avg_wei_sig(d),nrp)) ')'];
    avg_wei_ctrl_rp{d,1} = [num2str(roundsd(r_avg_wei_ctrl_pint(d),nrr)) ' (' num2str(roundsd(p_avg_wei_ctrl_pint(d),nrp)) ')'];
    avg_wei_schz_rp{d,1} = [num2str(roundsd(r_avg_wei_schz_pint(d),nrr)) ' (' num2str(roundsd(p_avg_wei_schz_pint(d),nrp)) ')'];
    
    % paired tests
    pair_eff_ctrl_rp{d,1} = [num2str(roundsd(rpair_eff_rp_ctrl_psub1(d),nrr)) ' (' num2str(roundsd(ppair_eff_rp_ctrl_psub1(d),nrp)) ')'];
    pair_eff_schz_rp{d,1} = [num2str(roundsd(rpair_eff_rp_schz_psub1(d),nrr)) ' (' num2str(roundsd(ppair_eff_rp_schz_psub1(d),nrp)) ')'];
    pair_trans_ctrl_rp{d,1} = [num2str(roundsd(rpair_trans_rp_ctrl_psub1(d),nrr)) ' (' num2str(roundsd(ppair_trans_rp_ctrl_psub1(d),nrp)) ')'];
    pair_trans_schz_rp{d,1} = [num2str(roundsd(rpair_trans_rp_schz_psub1(d),nrr)) ' (' num2str(roundsd(ppair_trans_rp_schz_psub1(d),nrp)) ')'];
    
    % paired tests normalised by random
    pair_eff_rand_ctrl_rp{d,1} = [num2str(roundsd(rpair_eff_rand_rp_ctrl_psub1(d),nrr)) ' (' num2str(roundsd(ppair_eff_rand_rp_ctrl_psub1(d),nrp)) ')'];
    pair_eff_rand_schz_rp{d,1} = [num2str(roundsd(rpair_eff_rand_rp_schz_psub1(d),nrr)) ' (' num2str(roundsd(ppair_eff_rand_rp_schz_psub1(d),nrp)) ')'];
    pair_trans_rand_ctrl_rp{d,1} = [num2str(roundsd(rpair_trans_rand_rp_ctrl_psub1(d),nrr)) ' (' num2str(roundsd(ppair_trans_rand_rp_ctrl_psub1(d),nrp)) ')'];
    pair_trans_rand_schz_rp{d,1} = [num2str(roundsd(rpair_trans_rand_rp_schz_psub1(d),nrr)) ' (' num2str(roundsd(ppair_trans_rand_rp_schz_psub1(d),nrp)) ')'];
    
end

% combine by measure
eff_all = [eff_rp eff_sig_rp eff_ctrl_rp eff_schz_rp];
eff_rand_all = [eff_rand_rp eff_rand_sig_rp eff_rand_ctrl_rp eff_rand_schz_rp];
trans_all = [trans_rp trans_sig_rp trans_ctrl_rp trans_schz_rp];
trans_rand_all = [trans_rand_rp trans_rand_sig_rp trans_rand_ctrl_rp trans_rand_schz_rp];
dcmp_n_all = [dcmp_n_rp dcmp_n_sig_rp dcmp_n_ctrl_rp dcmp_n_schz_rp];
avg_wei_all = [avg_wei_rp avg_wei_sig_rp avg_wei_ctrl_rp avg_wei_schz_rp];
pair_all = [pair_eff_ctrl_rp pair_eff_schz_rp pair_trans_ctrl_rp pair_trans_schz_rp];
pair_rand_all = [pair_eff_rand_ctrl_rp pair_eff_rand_schz_rp pair_trans_rand_ctrl_rp pair_trans_rand_schz_rp];
