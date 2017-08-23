% Code used to perform primary analyses described in the manuscript: 
%
% "Probabilistic thresholding of functional connectomes: application to schizophrenia"
% by Frantisek Vasa, Edward T. Bullmore and Ameera X. Patel
%
% Frantisek Vasa, October 2015 - August 2017
% fv247@cam.ac.uk
%
% Dependencies (available at https://github.com/frantisekvasa/structural_network_development):
%
% ranksum_effect_size       calculates ranksum effect size using the "simple difference formula" (Kerby, Compr. Psychol. 2014)
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
    temp_p = ctrl.fdr(:,:,i); %ctrl.did(i)
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
    temp_p = ctrl.fdr(:,:,i); %ctrl.did(i)
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
        temp_p = ctrl.fdr(:,:,i);
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
        temp_p = schz.fdr(:,:,i);
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

% find subject numbers between which ratios of "significant"/"unsignificant" participants cross 3:1 / 1:3
% ctrl
ctrl.r3to1(1) = min(find(cellfun('length',ctrl.did)<nc*(3/4)))-0.5;
ctrl.r3to1(2) = min(find(cellfun('length',ctrl.did)<nc*(1/4)))-0.5;
% schz
schz.r3to1(1) = min(find(cellfun('length',schz.did)<np*(3/4)))-0.5;
schz.r3to1(2) = min(find(cellfun('length',schz.did)<np*(1/4)))-0.5; 

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
    ctrl.m_eff_sig(i) = median(ctrl.eff(ctrl.did{i},i));
    ctrl.ci_eff_sig(i,:) = abs(prctile(ctrl.eff(ctrl.did{i},i),[prct_l prct_h])-ctrl.m_eff_sig(i));
    ctrl.m_trans_sig(i) = median(ctrl.trans(ctrl.did{i},i));
    ctrl.ci_trans_sig(i,:) = abs(prctile(ctrl.trans(ctrl.did{i},i),[prct_l prct_h])-ctrl.m_trans_sig(i));
    ctrl.m_dcmp_n_sig(i) = median(ctrl.dcmp_n(ctrl.did{i},i));
    ctrl.ci_dcmp_n_sig(i,:) = abs(prctile(ctrl.dcmp_n(ctrl.did{i},i),[prct_l prct_h])-ctrl.m_dcmp_n_sig(i));
    ctrl.m_dcmp_max_sig(i) = median(ctrl.dcmp_max(ctrl.did{i},i));
    ctrl.ci_dcmp_max_sig(i,:) = abs(prctile(ctrl.dcmp_max(ctrl.did{i},i),[prct_l prct_h])-ctrl.m_dcmp_max_sig(i));
    % schz
    schz.m_eff_sig(i) = median(schz.eff(schz.did{i},i));
    schz.ci_eff_sig(i,:) = abs(prctile(schz.eff(schz.did{i},i),[prct_l prct_h])-schz.m_eff_sig(i));
    schz.m_trans_sig(i) = median(schz.trans(schz.did{i},i));
    schz.ci_trans_sig(i,:) = abs(prctile(schz.trans(schz.did{i},i),[prct_l prct_h])-schz.m_trans_sig(i));
    schz.m_dcmp_n_sig(i) = median(schz.dcmp_n(schz.did{i},i));
    schz.ci_dcmp_n_sig(i,:) = abs(prctile(schz.dcmp_n(schz.did{i},i),[prct_l prct_h])-schz.m_dcmp_n_sig(i));
    schz.m_dcmp_max_sig(i) = median(schz.dcmp_max(schz.did{i},i));
    schz.ci_dcmp_max_sig(i,:) = abs(prctile(schz.dcmp_max(schz.did{i},i),[prct_l prct_h])-schz.m_dcmp_max_sig(i));
    
    % ranksum across densities
    p_eff_sig(i) = ranksum(ctrl.eff(ctrl.did{i},i),schz.eff(schz.did{i},i));
    p_trans_sig(i) = ranksum(ctrl.trans(ctrl.did{i},i),schz.trans(schz.did{i},i));
    p_dcmp_n_sig(i) = ranksum(ctrl.dcmp_n(ctrl.did{i},i),schz.dcmp_n(schz.did{i},i));
    p_dcmp_max_sig(i) = ranksum(ctrl.dcmp_max(ctrl.did{i},i),schz.dcmp_max(schz.did{i},i));
    
    % non-significant subjects
    % ctrl
    ctrl.m_eff_nonsig(i) = median(ctrl.eff(setdiff(1:nc,ctrl.did{i}),i));
    ctrl.ci_eff_nonsig(i,:) = abs(prctile(ctrl.eff(setdiff(1:nc,ctrl.did{i}),i),[prct_l prct_h])-ctrl.m_eff_nonsig(i));
    ctrl.m_trans_nonsig(i) = median(ctrl.trans(setdiff(1:nc,ctrl.did{i}),i));
    ctrl.ci_trans_nonsig(i,:) = abs(prctile(ctrl.trans(setdiff(1:nc,ctrl.did{i}),i),[prct_l prct_h])-ctrl.m_trans_nonsig(i));
    ctrl.m_dcmp_n_nonsig(i) = median(ctrl.dcmp_n(setdiff(1:nc,ctrl.did{i}),i));
    ctrl.ci_dcmp_n_nonsig(i,:) = abs(prctile(ctrl.dcmp_n(setdiff(1:nc,ctrl.did{i}),i),[prct_l prct_h])-ctrl.m_dcmp_n_nonsig(i));
    ctrl.m_dcmp_max_nonsig(i) = median(ctrl.dcmp_max(setdiff(1:nc,ctrl.did{i}),i));
    ctrl.ci_dcmp_max_nonsig(i,:) = abs(prctile(ctrl.dcmp_max(setdiff(1:nc,ctrl.did{i}),i),[prct_l prct_h])-ctrl.m_dcmp_max_nonsig(i));
    % schz
    schz.m_eff_nonsig(i) = median(schz.eff(setdiff(1:np,schz.did{i}),i));
    schz.ci_eff_nonsig(i,:) = abs(prctile(schz.eff(setdiff(1:np,schz.did{i}),i),[prct_l prct_h])-schz.m_eff_nonsig(i));
    schz.m_trans_nonsig(i) = median(schz.trans(setdiff(1:np,schz.did{i}),i));
    schz.ci_trans_nonsig(i,:) = abs(prctile(schz.trans(setdiff(1:np,schz.did{i}),i),[prct_l prct_h])-schz.m_trans_nonsig(i));
    schz.m_dcmp_n_nonsig(i) = median(schz.dcmp_n(setdiff(1:np,schz.did{i}),i));
    schz.ci_dcmp_n_nonsig(i,:) = abs(prctile(schz.dcmp_n(setdiff(1:np,schz.did{i}),i),[prct_l prct_h])-schz.m_dcmp_n_nonsig(i));
    schz.m_dcmp_max_nonsig(i) = median(schz.dcmp_max(setdiff(1:np,schz.did{i}),i));
    schz.ci_dcmp_max_nonsig(i,:) = abs(prctile(schz.dcmp_max(setdiff(1:np,schz.did{i}),i),[prct_l prct_h])-schz.m_dcmp_max_nonsig(i));
    
    % ranksum across densities
    if ~isempty(setdiff(1:nc,ctrl.did{i}))
        p_eff_ctrl(i) = ranksum(ctrl.eff(setdiff(1:nc,ctrl.did{i}),i),ctrl.eff(ctrl.did{i},i));
        p_trans_ctrl(i) = ranksum(ctrl.trans(setdiff(1:nc,ctrl.did{i}),i),ctrl.trans(ctrl.did{i},i));
        p_dcmp_n_ctrl(i) = ranksum(ctrl.dcmp_n(setdiff(1:nc,ctrl.did{i}),i),ctrl.dcmp_n(ctrl.did{i},i));
        p_dcmp_max_ctrl(i) = ranksum(ctrl.dcmp_max(setdiff(1:nc,ctrl.did{i}),i),ctrl.dcmp_max(ctrl.did{i},i));
    else
        p_eff_ctrl(i) = NaN;
        p_trans_ctrl(i) = NaN;
        p_dcmp_n_ctrl(i) = NaN;
        p_dcmp_max_ctrl(i) = NaN;
    end
    
    if ~isempty(setdiff(1:np,schz.did{i}))
        p_eff_schz(i) = ranksum(schz.eff(setdiff(1:np,schz.did{i}),i),schz.eff(schz.did{i},i));
        p_trans_schz(i) = ranksum(schz.trans(setdiff(1:np,schz.did{i}),i),schz.trans(schz.did{i},i));
        p_dcmp_n_schz(i) = ranksum(schz.dcmp_n(setdiff(1:np,schz.did{i}),i),schz.dcmp_n(schz.did{i},i));
        p_dcmp_max_schz(i) = ranksum(schz.dcmp_max(setdiff(1:np,schz.did{i}),i),schz.dcmp_max(schz.did{i},i));
    else
        p_eff_schz(i) = NaN;
        p_trans_schz(i) = NaN;
        p_dcmp_n_schz(i) = NaN;
        p_dcmp_max_schz(i) = NaN;
    end
    
end

%% significance of sig VS nonsig

nperm = 1000; % number of permutations

% efficiency
p_eff_sig_perm = zeros(ndens,nperm); r_eff_sig_perm = zeros(ndens,nperm);
p_eff_ctrl_perm = zeros(ndens,nperm); r_eff_ctrl_perm = zeros(ndens,nperm);
p_eff_schz_perm = zeros(ndens,nperm); r_eff_schz_perm = zeros(ndens,nperm);

% transitivity
p_trans_sig_perm = zeros(ndens,nperm); r_trans_sig_perm = zeros(ndens,nperm);
p_trans_ctrl_perm = zeros(ndens,nperm); r_trans_ctrl_perm = zeros(ndens,nperm);
p_trans_schz_perm = zeros(ndens,nperm); r_trans_schz_perm = zeros(ndens,nperm);

% # components
p_dcmp_n_sig_perm = zeros(ndens,nperm); r_dcmp_n_sig_perm = zeros(ndens,nperm);
p_dcmp_n_ctrl_perm = zeros(ndens,nperm); r_dcmp_n_ctrl_perm = zeros(ndens,nperm);
p_dcmp_n_schz_perm = zeros(ndens,nperm); r_dcmp_n_schz_perm = zeros(ndens,nperm);

% size of largest components
p_dcmp_max_sig_perm = zeros(ndens,nperm); r_dcmp_max_sig_perm = zeros(ndens,nperm);
p_dcmp_max_ctrl_perm = zeros(ndens,nperm); r_dcmp_max_ctrl_perm = zeros(ndens,nperm);
p_dcmp_max_schz_perm = zeros(ndens,nperm); r_dcmp_max_schz_perm = zeros(ndens,nperm);

for d = 2:1:ndens % hard-coded start at edge density 2, as at edge density 1 all controls are significant (IN COBRE DATA)
    disp(['edge density ' num2str(d) '%'])
    
    % empirical effect sizes
    
    % efficiency
    [p_eff_sig_emp(d),r_eff_sig_emp(d)] = ranksum_effect_size(ctrl.eff(ctrl.did{d},d),schz.eff(schz.did{d},d)); % ctrl VS schz (sig only)
    [p_eff_ctrl_emp(d),r_eff_ctrl_emp(d)] = ranksum_effect_size(ctrl.eff(ctrl.did{d},d),ctrl.eff(setdiff(1:nc,ctrl.did{d}),d)); % within group (sig VS nonsig)
    [p_eff_schz_emp(d),r_eff_schz_emp(d)] = ranksum_effect_size(schz.eff(schz.did{d},d),schz.eff(setdiff(1:np,schz.did{d}),d));
    
    % transitivity
    [p_trans_sig_emp(d),r_trans_sig_emp(d)] = ranksum_effect_size(ctrl.trans(ctrl.did{d},d),schz.trans(schz.did{d},d)); % ctrl VS schz (sig only)
    [p_trans_ctrl_emp(d),r_trans_ctrl_emp(d)] = ranksum_effect_size(ctrl.trans(ctrl.did{d},d),ctrl.trans(setdiff(1:nc,ctrl.did{d}),d)); % within group (sig VS nonsig)
    [p_trans_schz_emp(d),r_trans_schz_emp(d)] = ranksum_effect_size(schz.trans(schz.did{d},d),schz.trans(setdiff(1:np,schz.did{d}),d));
    
    % # components
    [p_dcmp_n_sig_emp(d),r_dcmp_n_sig_emp(d)] = ranksum_effect_size(ctrl.dcmp_n(ctrl.did{d},d),schz.dcmp_n(schz.did{d},d)); % ctrl VS schz (sig only)
    [p_dcmp_n_ctrl_emp(d),r_dcmp_n_ctrl_emp(d)] = ranksum_effect_size(ctrl.dcmp_n(ctrl.did{d},d),ctrl.dcmp_n(setdiff(1:nc,ctrl.did{d}),d)); % within group (sig VS nonsig)
    [p_dcmp_n_schz_emp(d),r_dcmp_n_schz_emp(d)] = ranksum_effect_size(schz.dcmp_n(schz.did{d},d),schz.dcmp_n(setdiff(1:np,schz.did{d}),d));
    
    % size of largest component
    [p_dcmp_max_sig_emp(d),r_dcmp_max_sig_emp(d)] = ranksum_effect_size(ctrl.dcmp_max(ctrl.did{d},d),schz.dcmp_max(schz.did{d},d)); % ctrl VS schz (sig only)
    [p_dcmp_max_ctrl_emp(d),r_dcmp_max_ctrl_emp(d)] = ranksum_effect_size(ctrl.dcmp_max(ctrl.did{d},d),ctrl.dcmp_max(setdiff(1:nc,ctrl.did{d}),d)); % within group (sig VS nonsig)
    [p_dcmp_max_schz_emp(d),r_dcmp_max_schz_emp(d)] = ranksum_effect_size(schz.dcmp_max(schz.did{d},d),schz.dcmp_max(setdiff(1:np,schz.did{d}),d));
    
    % permuted effect sizes
    for j = 1:nperm
        
        % permuted ids (length = sig)
        perm_idc = randsample(nc,length(ctrl.did{d}));
        perm_idp = randsample(np,length(schz.did{d}));
        
        % permuted ids complement (length = non-sig)
        perm_idc_comp = setdiff(1:nc,perm_idc);
        perm_idp_comp = setdiff(1:np,perm_idp);
        
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
    
end

% set ctrl "d=1" P=1 (as length(nonsig)==0)
pperm_eff_ctrl(1) = 1;
pperm_trans_ctrl(1) = 1;
pperm_dcmp_n_ctrl(1) = 1;
pperm_dcmp_max_ctrl(1) = 1;
pperm_nneg_ctrl(1) = 1;

%% invcount
% number of participants whose edges are all "significant" at a given edge density

ndens_c = floor(100*max(max(ctrl.dens)));
ndens_p = floor(100*max(max(schz.dens)));

for i = 1:1:ndens_c; ctrl.invcount(i) = sum(ctrl.dens>(i/100)); end;
for i = 1:1:ndens_p; schz.invcount(i) = sum(schz.dens>(i/100)); end;

%% obtain maximum p-value if keeping all subjects

maxp_dens = ndens;

for d = 1:1:maxp_dens
    disp(['edge density ' num2str(d) '%'])
    
    dens_all = d/100;
    nedge_all = ceil(dens_all*(nroi*(nroi-1)/2));
    
    for i = 1:1:nc
        temp_p = ctrl.fdr(:,:,i);
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
        temp_p = schz.fdr(:,:,i);
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

%% histogram bars for nparticipant ratios (n(sig) VS n(non-sig))

% ctrl
figure
b = bar(1:ndens,[cellfun('length',ctrl.did);nc-cellfun('length',ctrl.did)]','stack');
set(b(1),'FaceColor',ctrl.col,'FaceAlpha',0.7)
set(b(2),'FaceColor',ctrl.col_nonsig,'FaceAlpha',0.7)
xlim([0 ndens+1]); ylim([0 nc]); set(gca,'YTick',[0,nc]);
set(gca,'FontSize',fsize2); 
xlabel('edge density \rho (%)','FontSize',fsize,'FontName','Arial');
ylabel('# ctrl','FontSize',fsize,'FontName','Arial');
pos = get(gcf,'Position');
set(gcf,'Position', [pos(1) pos(2) 1.25*pos(3) 0.35*pos(4)]); set(gcf,'color','w');
print([plot_path 'ctrl_sig_vs_nonsig_bar_sc' sc '_alpha' num2str(100*alpha) '.png'],'-dpng')

% schz
figure
b = bar(1:ndens,[cellfun('length',schz.did);np-cellfun('length',schz.did)]','stack');
set(b(1),'FaceColor',schz.col,'FaceAlpha',0.7)
set(b(2),'FaceColor',schz.col_nonsig,'FaceAlpha',0.7)
xlim([0 ndens+1]); ylim([0 np]); set(gca,'YTick',[0,np]);
set(gca,'FontSize',fsize2); 
xlabel('edge density \rho (%)','FontSize',fsize,'FontName','Arial');
ylabel('# schz','FontSize',fsize,'FontName','Arial');
pos = get(gcf,'Position');
set(gcf,'Position', [pos(1) pos(2) 1.25*pos(3) 0.35*pos(4)]); set(gcf,'color','w');
print([plot_path 'schz_sig_vs_nonsig_bar_sc' sc '_alpha' num2str(100*alpha) '.png'],'-dpng')

% drop off curves
figure;
hold on
plot(1:ndens,ctrl.invcount(1:ndens),'Color',ctrl.col,'LineWidth',1);
plot(1:ndens,schz.invcount(1:ndens),'Color',schz.col,'LineWidth',1);
scatter(1:ndens,ctrl.invcount(1:ndens),26,'MarkerFaceColor',ctrl.col,'MarkerEdgeColor',ctrl.col); 
scatter(1:ndens,schz.invcount(1:ndens),26,'MarkerFaceColor',schz.col,'MarkerEdgeColor',schz.col);
hold off

set(gca,'FontSize',fsize2); box off
xlabel('edge density \rho (%)','FontSize',fsize,'FontName','Arial');
ylabel('# sig. participants','FontSize',fsize,'FontName','Arial');

pos = get(gcf,'Position');
set(gcf,'Position', [pos(1) pos(2) 1.25*pos(3) 0.70*pos(4)]); set(gcf,'color','w');

xlim([0 ndens+1]); 
set(gca,'YTick',[0,40,80])

print([plot_path 'drop_off_sc' sc '_alpha' num2str(100*alpha) '.png'],'-dpng')

% max(P_FDR)
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

%% plots - emp (eff, trans)

fsize = 30;
fsize2 = 20;

%%%%% ctrl VS schz

%%% all subjects

% efficiency
boundedline_plot(1:ndens,ctrl.m_eff,ctrl.ci_eff,schz.m_eff,schz.ci_eff,...
    ctrl.col,schz.col,alph,NaN,p_eff,'edge density \rho (%)','efficiency','ctrl','schz',[0 ndens+1],[0 0.8],fsize,fsize2,...
    plot_path,['eff_ctrl_vs_schz_all_boundedline_sc' sc '_alpha' num2str(100*alpha)]);

% transitivity
boundedline_plot(1:ndens,ctrl.m_trans,ctrl.ci_trans,schz.m_trans,schz.ci_trans,...
    ctrl.col,schz.col,alph,NaN,p_trans,'edge density \rho (%)','transitivity','ctrl','schz',[0 ndens+1],[0.3 0.9],fsize,fsize2,...
    plot_path,['trans_ctrl_vs_schz_all_boundedline_sc' sc '_alpha' num2str(100*alpha)]);

%%% significant subjects only

% efficiency
boundedline_plot_2p_drop(1:ndens,ctrl.m_eff_sig,ctrl.ci_eff_sig,schz.m_eff_sig,schz.ci_eff_sig,ctrl.invcount(1:ndens),schz.invcount(1:ndens),...
    ctrl.col,schz.col,alph,NaN,p_eff_sig,pperm_eff_sig,'edge density \rho (%)','efficiency','# sig. participants','ctrl (sig)','schz (sig)',[0 ndens+1],[0 0.8],[0 80],fsize,fsize2,...
    plot_path,['eff_ctrl_vs_schz_sig_boundedline_sc' sc '_alpha' num2str(100*alpha)]);

% transitivity
boundedline_plot_2p_drop(1:ndens,ctrl.m_trans_sig,ctrl.ci_trans_sig,schz.m_trans_sig,schz.ci_trans_sig,ctrl.invcount(1:ndens),schz.invcount(1:ndens),...
    ctrl.col,schz.col,alph,NaN,p_trans_sig,pperm_trans_sig,'edge density \rho (%)','transitivity','# sig. participants','ctrl (sig)','schz (sig)',[0 ndens+1],[0.3 0.9],[0 80],fsize,fsize2,...
    plot_path,['trans_ctrl_vs_schz_sig_boundedline_sc' sc '_alpha' num2str(100*alpha)]);

%%%%% significant VS non-significant subjects 

%%% ctrl

% efficiency
boundedline_plot_2p(1:ndens,ctrl.m_eff_sig,ctrl.ci_eff_sig,ctrl.m_eff_nonsig,ctrl.ci_eff_nonsig,...
    ctrl.col,ctrl.col_nonsig,alph,ctrl.r3to1,p_eff_ctrl,pperm_eff_ctrl,'edge density \rho (%)','efficiency','significant','non-significant',[0 ndens+1],[0 0.8],fsize,fsize2,...
    plot_path,['eff_sig_vs_nonsig_ctrl_boundedline_sc' sc '_alpha' num2str(100*alpha)]);

% transitivity
boundedline_plot_2p(1:ndens,ctrl.m_trans_sig,ctrl.ci_trans_sig,ctrl.m_trans_nonsig,ctrl.ci_trans_nonsig,...
    ctrl.col,ctrl.col_nonsig,alph,ctrl.r3to1,p_trans_ctrl,pperm_trans_ctrl,'edge density \rho (%)','transitivity','significant','non-significant',[0 ndens+1],[0.3 0.9],fsize,fsize2,...
    plot_path,['trans_sig_vs_nonsig_ctrl_boundedline_sc' sc '_alpha' num2str(100*alpha)]);

%%% schz

% efficiency
boundedline_plot_2p(1:ndens,schz.m_eff_sig,schz.ci_eff_sig,schz.m_eff_nonsig,schz.ci_eff_nonsig,...
    schz.col,schz.col_nonsig,alph,schz.r3to1,p_eff_schz,pperm_eff_schz,'edge density \rho (%)','efficiency','significant','non-significant',[0 ndens+1],[0 0.8],fsize,fsize2,...
    plot_path,['eff_sig_vs_nonsig_schz_boundedline_sc' sc '_alpha' num2str(100*alpha)]);

% transitivity
boundedline_plot_2p(1:ndens,schz.m_trans_sig,schz.ci_trans_sig,schz.m_trans_nonsig,schz.ci_trans_nonsig,...
    schz.col,schz.col_nonsig,alph,schz.r3to1,p_trans_schz,pperm_trans_schz,'edge density \rho (%)','transitivity','significant','non-significant',[0 ndens+1],[0.3 0.9],fsize,fsize2,...
    plot_path,['trans_sig_vs_nonsig_schz_boundedline_sc' sc '_alpha' num2str(100*alpha)]);

%% plots - emp (cmp)

%%%%% ctrl VS schz

%%% all subjects

% # components
boundedline_plot(1:ndens,ctrl.m_dcmp_n,ctrl.ci_dcmp_n,schz.m_dcmp_n,schz.ci_dcmp_n,...
    ctrl.col,schz.col,alph,NaN,p_dcmp_n,'edge density \rho (%)','# components','ctrl','schz',[0 ndens+1],[0 150],fsize,fsize2,...
    plot_path,['dcmp_n_ctrl_vs_schz_all_boundedline_sc' sc '_alpha' num2str(100*alpha)]);

% size of largest component
boundedline_plot(1:ndens,ctrl.m_dcmp_max,ctrl.ci_dcmp_max,schz.m_dcmp_max,schz.ci_dcmp_max,...
    ctrl.col,schz.col,alph,NaN,p_dcmp_max,'edge density \rho (%)','max(comp. size)','ctrl','schz',[0 ndens+1],[200 450],fsize,fsize2,...
    plot_path,['dcmp_max_ctrl_vs_schz_all_boundedline_sc' sc '_alpha' num2str(100*alpha)]);

%%% significant subjects only

% # components
boundedline_plot_2p_drop(1:ndens,ctrl.m_dcmp_n_sig,ctrl.ci_dcmp_n_sig,schz.m_dcmp_n_sig,schz.ci_dcmp_n_sig,ctrl.invcount(1:ndens),schz.invcount(1:ndens),...
    ctrl.col,schz.col,alph,NaN,p_dcmp_n_sig,pperm_dcmp_n_sig,'edge density \rho (%)','# components','# sig. participants','ctrl (sig)','schz (sig)',[0 ndens+1],[0 150],[0 80],fsize,fsize2,...
    plot_path,['dcmp_n_ctrl_vs_schz_sig_boundedline_sc' sc '_alpha' num2str(100*alpha)]);

% size of largest component
boundedline_plot_2p_drop(1:ndens,ctrl.m_dcmp_max_sig,ctrl.ci_dcmp_max_sig,schz.m_dcmp_max_sig,schz.ci_dcmp_max_sig,ctrl.invcount(1:ndens),schz.invcount(1:ndens),...
    ctrl.col,schz.col,alph,NaN,p_dcmp_max_sig,pperm_dcmp_max_sig,'edge density \rho (%)','max(comp. size)','# sig. participants','ctrl (sig)','schz (sig)',[0 ndens+1],[200 450],[0 80],fsize,fsize2,...
    plot_path,['dcmp_max_ctrl_vs_schz_sig_boundedline_sc' sc '_alpha' num2str(100*alpha)]);

%%%%% significant VS non-significant subjects 

%%% ctrl

% # components
boundedline_plot_2p(1:ndens,ctrl.m_dcmp_n_sig,ctrl.ci_dcmp_n_sig,ctrl.m_dcmp_n_nonsig,ctrl.ci_dcmp_n_nonsig,...
    ctrl.col,ctrl.col_nonsig,alph,ctrl.r3to1,p_dcmp_n_ctrl,pperm_dcmp_n_ctrl,'edge density \rho (%)','# components','significant','non-significant',[0 ndens+1],[0 150],fsize,fsize2,...
    plot_path,['dcmp_n_sig_vs_nonsig_ctrl_boundedline_sc' sc '_alpha' num2str(100*alpha)]);

% size of largest component
boundedline_plot_2p(1:ndens,ctrl.m_dcmp_max_sig,ctrl.ci_dcmp_max_sig,ctrl.m_dcmp_max_nonsig,ctrl.ci_dcmp_max_nonsig,...
    ctrl.col,ctrl.col_nonsig,alph,ctrl.r3to1,p_dcmp_max_ctrl,pperm_dcmp_max_ctrl,'edge density \rho (%)','max(comp. size)','significant','non-significant',[0 ndens+1],[200 450],fsize,fsize2,...
    plot_path,['dcmp_max_sig_vs_nonsig_ctrl_boundedline_sc' sc '_alpha' num2str(100*alpha)]);

%%% schz

% # components
boundedline_plot_2p(1:ndens,schz.m_dcmp_n_sig,schz.ci_dcmp_n_sig,schz.m_dcmp_n_nonsig,schz.ci_dcmp_n_nonsig,...
    schz.col,schz.col_nonsig,alph,schz.r3to1,p_dcmp_n_schz,pperm_dcmp_n_schz,'edge density \rho (%)','# components','significant','non-significant',[0 ndens+1],[0 150],fsize,fsize2,...
    plot_path,['dcmp_n_sig_vs_nonsig_schz_boundedline_sc' sc '_alpha' num2str(100*alpha)]);

% size of largest component
boundedline_plot_2p(1:ndens,schz.m_dcmp_max_sig,schz.ci_dcmp_max_sig,schz.m_dcmp_max_nonsig,schz.ci_dcmp_max_nonsig,...
    schz.col,schz.col_nonsig,alph,schz.r3to1,p_dcmp_max_schz,pperm_dcmp_max_schz,'edge density \rho (%)','max(comp. size)','significant','non-significant',[0 ndens+1],[200 450],fsize,fsize2,...
    plot_path,['dcmp_max_sig_vs_nonsig_schz_boundedline_sc' sc '_alpha' num2str(100*alpha)]);

%% threshold to fixed density by R (across range of densities)

ndens = 35;

for d = 1:1:ndens % density cut-off %ndens
    disp(['edge density ' num2str(d) '%'])
    
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
        
        clear temp_r temp_mat_dens
    end
    
    for i = 1:1:np %length(schz.did)
        temp_r = schz.r(:,:,i); %schz.did(i)
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
        
        clear temp_r temp_mat_dens
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
    
    % ranksum (non-parametric test of difference) across densities r-thr ctrl VS schz
    p_eff_r(i) = ranksum(ctrl.eff_r(:,i),schz.eff_r(:,i));
    p_trans_r(i) = ranksum(ctrl.trans_r(:,i),schz.trans_r(:,i));
    p_dcmp_n_r(i) = ranksum(ctrl.dcmp_n_r(:,i),schz.dcmp_n_r(:,i));
    p_dcmp_max_r(i) = ranksum(ctrl.dcmp_max_r(:,i),schz.dcmp_max_r(:,i));
    
end

for i = 1:1:ndens
    
    % ranksum across densities r-thr VS p-thr
    % ctrl
    p_eff_rp_ctrl(i) = ranksum(ctrl.eff_r(:,i),ctrl.eff(:,i));
    p_trans_rp_ctrl(i) = ranksum(ctrl.trans_r(:,i),ctrl.trans(:,i));
    p_dcmp_n_rp_ctrl(i) = ranksum(ctrl.dcmp_n_r(:,i),ctrl.dcmp_n(:,i));
    p_dcmp_max_rp_ctrl(i) = ranksum(ctrl.dcmp_max_r(:,i),ctrl.dcmp_max(:,i));
    
    % schz
    p_eff_rp_schz(i) = ranksum(schz.eff_r(:,i),schz.eff(:,i));
    p_trans_rp_schz(i) = ranksum(schz.trans_r(:,i),schz.trans(:,i));
    p_dcmp_n_rp_schz(i) = ranksum(schz.dcmp_n_r(:,i),schz.dcmp_n(:,i));
    p_dcmp_max_rp_schz(i) = ranksum(schz.dcmp_max_r(:,i),schz.dcmp_max(:,i));
    
end

for i = 1:1:ndens
    
    % signrank across densities r-thr VS p-thr
    % ctrl
    ppair_eff_rp_ctrl(i) = signrank(ctrl.eff_r(:,i),ctrl.eff(:,i));
    ppair_trans_rp_ctrl(i) = signrank(ctrl.trans_r(:,i),ctrl.trans(:,i));
    ppair_dcmp_n_rp_ctrl(i) = signrank(ctrl.dcmp_n_r(:,i),ctrl.dcmp_n(:,i));
    ppair_dcmp_max_rp_ctrl(i) = signrank(ctrl.dcmp_max_r(:,i),ctrl.dcmp_max(:,i));
    
    % schz
    ppair_eff_rp_schz(i) = signrank(schz.eff_r(:,i),schz.eff(:,i));
    ppair_trans_rp_schz(i) = signrank(schz.trans_r(:,i),schz.trans(:,i));
    ppair_dcmp_n_rp_schz(i) = signrank(schz.dcmp_n_r(:,i),schz.dcmp_n(:,i));
    ppair_dcmp_max_rp_schz(i) = signrank(schz.dcmp_max_r(:,i),schz.dcmp_max(:,i));
    
end

% within participant difference between r-thr and p-thr
% ctrl
ctrl.d_eff = ctrl.eff-ctrl.eff_r;
ctrl.d_trans = ctrl.trans-ctrl.trans_r;
ctrl.d_dcmp_n = ctrl.dcmp_n-ctrl.dcmp_n_r;

ctrl.m_d_eff = median(ctrl.d_eff);
ctrl.m_d_trans = median(ctrl.d_trans);
ctrl.m_d_dcmp_n = median(ctrl.d_dcmp_n);

for i = 1:1:ndens
    ctrl.ci_d_eff(i,:) = prctile(ctrl.d_eff(:,i),[prct_l prct_h])-ctrl.m_d_eff(i);
    ctrl.ci_d_trans(i,:) = prctile(ctrl.d_trans(:,i),[prct_l prct_h])-ctrl.m_d_trans(i);
    ctrl.ci_d_dcmp_n(i,:) = prctile(ctrl.d_dcmp_n(:,i),[prct_l prct_h])-ctrl.m_d_dcmp_n(i);
end

ctrl.ci_d_eff = abs(ctrl.ci_d_eff); ctrl.ci_d_trans = abs(ctrl.ci_d_trans); ctrl.ci_d_dcmp_n = abs(ctrl.ci_d_dcmp_n);

% schz
schz.d_eff = schz.eff-schz.eff_r;
schz.d_trans = schz.trans-schz.trans_r;
schz.d_dcmp_n = schz.dcmp_n-schz.dcmp_n_r;

schz.m_d_eff = median(schz.d_eff);
schz.m_d_trans = median(schz.d_trans);
schz.m_d_dcmp_n = median(schz.d_dcmp_n);

for i = 1:1:ndens
    schz.ci_d_eff(i,:) = prctile(schz.d_eff(:,i),[prct_l prct_h])-schz.m_d_eff(i);
    schz.ci_d_trans(i,:) = prctile(schz.d_trans(:,i),[prct_l prct_h])-schz.m_d_trans(i);
    schz.ci_d_dcmp_n(i,:) = prctile(schz.d_dcmp_n(:,i),[prct_l prct_h])-schz.m_d_dcmp_n(i);
end

schz.ci_d_eff = abs(schz.ci_d_eff); schz.ci_d_trans = abs(schz.ci_d_trans); schz.ci_d_dcmp_n = abs(schz.ci_d_dcmp_n);

%% compare P and R thresholded connectomes

% number of edges that are the same
for d = 1:1:ndens
    disp(['edge density ' num2str(d) '%'])
        
    % ctrl
    for i = 1:1:nc
        ctrl.rp_overlap(d,i) = length(intersect(intersect(find(ctrl.mat_dens_r(:,:,d,i)),triu_ind),intersect(find(ctrl.mat_dens(:,:,d,i)),triu_ind)))/length(intersect(find(ctrl.mat_dens(:,:,d,i)),triu_ind));
    end
    
    % schz
    for i = 1:1:np
        schz.rp_overlap(d,i) = length(intersect(intersect(find(schz.mat_dens_r(:,:,d,i)),triu_ind),intersect(find(schz.mat_dens(:,:,d,i)),triu_ind)))/length(intersect(find(schz.mat_dens(:,:,d,i)),triu_ind));
    end
    
end

fise = 30;
fsize2 = 20;

% plots

%%% RELATIVE

figure;
hold on
for i = 1:1:nc; hc = plot(1:ndens,100*(1-ctrl.rp_overlap(:,i)),'Color',ctrl.col_light,'LineWidth',lwd_ind); end
for i = 1:1:np; hp = plot(1:ndens,100*(1-schz.rp_overlap(:,i)),'Color',schz.col_light,'LineWidth',lwd_ind); end
hca = plot(1:ndens,100*(1-mean(ctrl.rp_overlap,2)),'Color',ctrl.col*col_ind,'LineWidth',lwd);
hpa = plot(1:ndens,100*(1-mean(schz.rp_overlap,2)),'Color',schz.col*0.8,'LineWidth',lwd);
hold off

set(gca,'FontSize',fsize2); box off
xlabel('edge density \rho (%)','FontSize',fsize,'FontName','Arial');
ylabel('r/P thr. edge diff. (%)','FontSize',fsize,'FontName','Arial');

pos = get(gcf,'Position');
set(gcf,'Position', [pos(1) pos(2) 1.25*pos(3) 0.70*pos(4)]); set(gcf,'color','w');

xlim([0 ndens+1]); 
ylim([-1,65])
set(gca,'YTick',[0:15:60])

% inset limit box
rpi = -1; rpf = 11;
rhoi = 0.95; rhof = 5.05;

hold on
plot([rhoi rhof],[rpi rpi],'k','LineWidth',2) % horiz bottom
plot([rhoi rhof],[rpf rpf],'k','LineWidth',2) % horiz top
plot([rhoi rhoi],[rpi rpf],'k','LineWidth',2) % vert left
plot([rhof rhof],[rpi rpf],'k','LineWidth',2) % vert right
hold off;

print([plot_path 'cobre_rp_diff_full.png'],'-dpng')

% inset
figure;
hold on
for i = 1:1:nc; hc = plot(1:ndens,100*(1-ctrl.rp_overlap(:,i)),'Color',ctrl.col_light,'LineWidth',lwd_ind); end
for i = 1:1:np; hp = plot(1:ndens,100*(1-schz.rp_overlap(:,i)),'Color',schz.col_light,'LineWidth',lwd_ind); end
hca = plot(1:ndens,100*(1-mean(ctrl.rp_overlap,2)),'Color',ctrl.col*col_ind,'LineWidth',lwd);
hpa = plot(1:ndens,100*(1-mean(schz.rp_overlap,2)),'Color',schz.col*0.8,'LineWidth',lwd);
hold off

set(gca,'FontSize',fsize2);

xlim([rhoi rhof])
ylim([rpi rpf])

set(gca,'FontSize',fsize2); %box off
set(gca,'YTick',[0 5 10])
set(gca,'XTick',[1 3 5])
set(gcf,'color','w');

print([plot_path 'cobre_rp_diff_inset.png'],'-dpng')

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

%% r VS P

fsize = 30;
fsize2 = 20;

% ctrl (r VS P)

% efficiency
boundedline_plot(1:ndens,ctrl.m_eff_r,ctrl.ci_eff_r,ctrl.m_eff,ctrl.ci_eff,...
    ctrl.col_r,ctrl.col,alph,NaN,ppair_eff_rp_ctrl,'edge density \rho (%)','efficiency','ctrl r','ctrl p',[0 ndens+1],[0 0.8],fsize,fsize2,...
    plot_path,['eff_ctrl_r_vs_p_boundedline_sc' sc '_alpha' num2str(100*alpha)]);

% transitivity
boundedline_plot(1:ndens,ctrl.m_trans_r,ctrl.ci_trans_r,ctrl.m_trans,ctrl.ci_trans,...
    ctrl.col_r,ctrl.col,alph,NaN,ppair_trans_rp_ctrl,'edge density \rho (%)','transitivity','ctrl r','ctrl p',[0 ndens+1],[0.3 0.9],fsize,fsize2,...
    plot_path,['trans_ctrl_r_vs_p_boundedline_sc' sc '_alpha' num2str(100*alpha)]);

% # components
boundedline_plot(1:ndens,ctrl.m_dcmp_n_r,ctrl.ci_dcmp_n_r,ctrl.m_dcmp_n,ctrl.ci_dcmp_n,...
    ctrl.col_r,ctrl.col,alph,NaN,ppair_dcmp_n_rp_ctrl,'edge density \rho (%)','# conn. comp.','ctrl r','ctrl p',[0 ndens+1],[0 150],fsize,fsize2,...
    plot_path,['dcmp_n_ctrl_r_vs_p_boundedline_sc' sc '_alpha' num2str(100*alpha)]);

% size of largest component
boundedline_plot(1:ndens,ctrl.m_dcmp_max_r,ctrl.ci_dcmp_max_r,ctrl.m_dcmp_max,ctrl.ci_dcmp_max,...
    ctrl.col_r,ctrl.col,alph,NaN,ppair_dcmp_max_rp_ctrl,'edge density \rho (%)','max(comp. size)','ctrl r','ctrl p',[0 ndens+1],[200 450],fsize,fsize2,...
    plot_path,['dcmp_max_ctrl_r_vs_p_rthr_boundedline_sc' sc '_alpha' num2str(100*alpha)]);

%%% difference

% efficiency
boundedline_plot_one(1:ndens,ctrl.m_d_eff,ctrl.ci_d_eff,rgb('darkorchid'),alph,ppair_eff_rp_ctrl,'edge density \rho (%)','\Delta efficiency',[0 ndens+1],[-0.004 0.004],fsize,fsize2,...
    plot_path,['d_eff_ctrl_r_vs_p_boundedline_sc' sc '_alpha' num2str(100*alpha)]);

% transitivity
boundedline_plot_one(1:ndens,ctrl.m_d_trans,ctrl.ci_d_trans,rgb('darkorchid'),alph,ppair_trans_rp_ctrl,'edge density \rho (%)','\Delta transitivity',[0 ndens+1],[-0.006 0.006],fsize,fsize2,...
    plot_path,['d_trans_ctrl_r_vs_p_boundedline_sc' sc '_alpha' num2str(100*alpha)]);

% # components
boundedline_plot_one(1:ndens,ctrl.m_d_dcmp_n,ctrl.ci_d_dcmp_n,rgb('darkorchid'),alph,ppair_dcmp_n_rp_ctrl,'edge density \rho (%)','\Delta # conn. comp.',[0 ndens+1],[-1.5 2.5],fsize,fsize2,...
    plot_path,['d_dcmp_n_ctrl_r_vs_p_boundedline_sc' sc '_alpha' num2str(100*alpha)]);

% schz (r VS P)

% efficiency
boundedline_plot(1:ndens,schz.m_eff_r,schz.ci_eff_r,schz.m_eff,schz.ci_eff,...
    schz.col_r,schz.col,alph,NaN,ppair_eff_rp_schz,'edge density \rho (%)','efficiency','schz r','schz p',[0 ndens+1],[0 0.8],fsize,fsize2,...
    plot_path,['eff_schz_r_vs_p_boundedline_sc' sc '_alpha' num2str(100*alpha)]);

% transitivity
boundedline_plot(1:ndens,schz.m_trans_r,schz.ci_trans_r,schz.m_trans,schz.ci_trans,...
    schz.col_r,schz.col,alph,NaN,ppair_trans_rp_schz,'edge density \rho (%)','transitivity','schz r','schz p',[0 ndens+1],[0.3 0.9],fsize,fsize2,...
    plot_path,['trans_schz_r_vs_p_boundedline_sc' sc '_alpha' num2str(100*alpha)]);

% # components
boundedline_plot(1:ndens,schz.m_dcmp_n_r,schz.ci_dcmp_n_r,schz.m_dcmp_n,schz.ci_dcmp_n,...
    schz.col_r,schz.col,alph,NaN,ppair_dcmp_n_rp_schz,'edge density \rho (%)','# conn. comp.','schz r','schz p',[0 ndens+1],[0 150],fsize,fsize2,...
    plot_path,['dcmp_n_schz_r_vs_p_boundedline_sc' sc '_alpha' num2str(100*alpha)]);

% size of largest component
boundedline_plot(1:ndens,schz.m_dcmp_max_r,schz.ci_dcmp_max_r,schz.m_dcmp_max,schz.ci_dcmp_max,...
    schz.col_r,schz.col,alph,NaN,ppair_dcmp_max_rp_schz,'edge density \rho (%)','max(comp. size)','schz r','schz p',[0 ndens+1],[200 450],fsize,fsize2,...
    plot_path,['dcmp_max_schz_r_vs_p_rthr_boundedline_sc' sc '_alpha' num2str(100*alpha)]);

%%% difference

% efficiency
boundedline_plot_one(1:ndens,schz.m_d_eff,schz.ci_d_eff,rgb('orangered'),alph,ppair_eff_rp_schz,'edge density \rho (%)','\Delta efficiency',[0 ndens+1],[-0.007 0.007],fsize,fsize2,...
    plot_path,['d_eff_schz_r_vs_p_boundedline_sc' sc '_alpha' num2str(100*alpha)]);

% transitivity
boundedline_plot_one(1:ndens,schz.m_d_trans,schz.ci_d_trans,rgb('orangered'),alph,ppair_trans_rp_schz,'edge density \rho (%)','\Delta transitivity',[0 ndens+1],[-0.01 0.01],fsize,fsize2,...
    plot_path,['d_trans_schz_r_vs_p_boundedline_sc' sc '_alpha' num2str(100*alpha)]);

% # components
boundedline_plot_one(1:ndens,schz.m_d_dcmp_n,schz.ci_d_dcmp_n,rgb('orangered'),alph,ppair_dcmp_n_rp_schz,'edge density \rho (%)','\Delta # conn. comp.',[0 ndens+1],[-1.5 3.5],fsize,fsize2,...
    plot_path,['d_dcmp_n_schz_r_vs_p_boundedline_sc' sc '_alpha' num2str(100*alpha)]);

% ctrl VS schz - r-based

% efficiency
boundedline_plot(1:ndens,ctrl.m_eff_r,ctrl.ci_eff_r,schz.m_eff_r,schz.ci_eff_r,...
    ctrl.col,schz.col,alph,NaN,p_eff_r,'edge density \rho (%)','efficiency','ctrl','schz',[0 ndens+1],[0 0.8],fsize,fsize2,...
    plot_path,['eff_ctrl_vs_schz_all_rthr_boundedline_sc' sc '_alpha' num2str(100*alpha)]);

% transitivity
boundedline_plot(1:ndens,ctrl.m_trans_r,ctrl.ci_trans_r,schz.m_trans_r,schz.ci_trans_r,...
    ctrl.col,schz.col,alph,NaN,p_trans_r,'edge density \rho (%)','transitivity','ctrl','schz',[0 ndens+1],[0.3 0.9],fsize,fsize2,...
    plot_path,['trans_ctrl_vs_schz_all_rthr_boundedline_sc' sc '_alpha' num2str(100*alpha)]);

% # components
boundedline_plot(1:ndens,ctrl.m_dcmp_n_r,ctrl.ci_dcmp_n_r,schz.m_dcmp_n_r,schz.ci_dcmp_n_r,...
    ctrl.col,schz.col,alph,NaN,p_dcmp_n_r,'edge density \rho (%)','# conn. comp.','ctrl','schz',[0 ndens+1],[0 150],fsize,fsize2,...
    plot_path,['dcmp_n_ctrl_vs_schz_all_rthr_boundedline_sc' sc '_alpha' num2str(100*alpha)]);

% size of largest component
boundedline_plot(1:ndens,ctrl.m_dcmp_max_r,ctrl.ci_dcmp_max_r,schz.m_dcmp_max_r,schz.ci_dcmp_max_r,...
    ctrl.col,schz.col,alph,NaN,p_dcmp_max_r,'edge density \rho (%)','max(comp. size)','ctrl','schz',[0 ndens+1],[200 450],fsize,fsize2,...
    plot_path,['dcmp_max_ctrl_vs_schz_all_rthr_boundedline_sc' sc '_alpha' num2str(100*alpha)]);

