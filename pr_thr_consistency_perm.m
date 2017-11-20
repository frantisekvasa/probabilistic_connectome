% permute r and p matrices within group, 10000 times

% required variables

% ctrl_psub1_id
% schz_psub1_id
% ctrl_sig_id
% schz_sig_id
% ctrl_mat_dens_r
% ctrl_mat_dens
% schz_mat_dens_r
% schz_mat_dens

%% cutoff: P=1

nperm = 10000;

ctrl_tbl_d_perm = zeros(ndens,nc+1,nperm);
schz_tbl_d_perm = zeros(ndens,np+1,nperm);
for p = 1:1:nperm
    p%if mod(p,10) == 0; disp(p); end
    for d = 1:1:ndens
        %d
        if mod(d,10) == 0; disp(d); end
        
        % number of patients and controls at this density
        nc_t = length(ctrl_psub1_id{d});
        np_t = length(schz_psub1_id{d});
        
        % randomly select "nc" (and "np") matrices from combined set of r- and P-thresholded matrices
        % controls
        set1_c = sort(ctrl_psub1_id{d}(randperm(nc_t,randi(nc_t-1))))'; % "-1" so that set 2 isn't null
        set2_c = sort(setdiff(ctrl_psub1_id{d},set1_c))';
        
        % patients
        set1_p = sort(schz_psub1_id{d}(randperm(np_t,randi(np_t-1))))'; % "-1" so that set 2 isn't null
        set2_p = sort(setdiff(schz_psub1_id{d},set1_p))';
        
        % ctrl
        ctrl_r_overlap_perm(:,:,d) = sum(squeeze(logical(ctrl_mat_dens_r(:,:,d,set1_c))),3) + sum(squeeze(logical(ctrl_mat_dens(:,:,d,set2_c))),3);%/length(ctrl_sig_id{d});
        ctrl_p_overlap_perm(:,:,d) = sum(squeeze(logical(ctrl_mat_dens(:,:,d,set1_c))),3) + sum(squeeze(logical(ctrl_mat_dens_r(:,:,d,set2_c))),3);%/length(ctrl_sig_id{d});
        % schz
        schz_r_overlap_perm(:,:,d) = sum(squeeze(logical(schz_mat_dens_r(:,:,d,set1_p))),3) + sum(squeeze(logical(schz_mat_dens(:,:,d,set2_p))),3);%/length(schz_sig_id{d});
        schz_p_overlap_perm(:,:,d) = sum(squeeze(logical(schz_mat_dens(:,:,d,set1_p))),3) + sum(squeeze(logical(schz_mat_dens_r(:,:,d,set2_p))),3);%/length(schz_sig_id{d});
        
    end
    
    % frequencies of overlap
    for d = 1:1:ndens
        %d
        
        %%% ctrl
        % r
        temp_ctrl = ctrl_r_overlap_perm(:,:,d);
        tbl1_r = tabulate(temp_ctrl(triu_ind)); % obtain frequency of frequencies (how many edges appear twice, three times, four times, ...)
        ctrl_tbl_r_perm(d,tbl1_r(:,1)+1) = tbl1_r(:,2);
        % p
        temp_ctrl = ctrl_p_overlap_perm(:,:,d);
        tbl1_p = tabulate(temp_ctrl(triu_ind)); % obtain frequency of frequencies (how many edges appear twice, three times, four times, ...)
        ctrl_tbl_p_perm(d,tbl1_p(:,1)+1) = tbl1_p(:,2);
        %%% schz
        % r
        temp_ctrl = schz_r_overlap_perm(:,:,d);
        tbl1_r = tabulate(temp_ctrl(triu_ind)); % obtain frequency of frequencies (how many edges appear twice, three times, four times, ...)
        schz_tbl_r_perm(d,tbl1_r(:,1)+1) = tbl1_r(:,2);
        % p
        temp_ctrl = schz_p_overlap_perm(:,:,d);
        tbl1_p = tabulate(temp_ctrl(triu_ind)); % obtain frequency of frequencies (how many edges appear twice, three times, four times, ...)
        schz_tbl_p_perm(d,tbl1_p(:,1)+1) = tbl1_p(:,2);

    end
    
    % difference
    ctrl_tbl_d_perm(:,:,p) = ctrl_tbl_p_perm-ctrl_tbl_r_perm;
    schz_tbl_d_perm(:,:,p) = schz_tbl_p_perm-schz_tbl_r_perm;
    
end

save('pr_thr_consistency_perm_psub1.mat','ctrl_tbl_d_perm','ctrl_tbl_d_perm')

% compare empirical difference to permuted differences

% ctrl
for d = 1:1:ndens 
    for i = 1:1:(length(ctrl_psub1_id{d})+1)
        if ctrl_tbl_d(d,i) >= 0
            ctrl_tbl_d_p(d,i) = sum(ctrl_tbl_d_perm(d,i,:)>ctrl_tbl_d(d,i))/nperm;
        elseif ctrl_tbl_d(d,i) < 0
            ctrl_tbl_d_p(d,i) = sum(ctrl_tbl_d_perm(d,i,:)<ctrl_tbl_d(d,i))/nperm;
        end
    end
end
    
% schz
for d = 1:1:ndens
    for i = 1:1:(length(schz_psub1_id{d})+1)
        if schz_tbl_d(d,i) >= 0
            schz_tbl_d_p(d,i) = sum(schz_tbl_d_perm(d,i,:)>schz_tbl_d(d,i))/nperm;
        elseif schz_tbl_d(d,i) < 0
            schz_tbl_d_p(d,i) = sum(schz_tbl_d_perm(d,i,:)<schz_tbl_d(d,i))/nperm;
        end
    end
end

save('pr_thr_consistency_perm_psub1.mat','ctrl_tbl_d_p','ctrl_tbl_d_p')

%% cutoff: Psig

nperm = 10000;

ctrl_tbl_d_perm = zeros(ndens,nc+1,nperm);
schz_tbl_d_perm = zeros(ndens,np+1,nperm);
for p = 1:1:nperm
    p%if mod(p,10) == 0; disp(p); end
    for d = 1:1:ndens
        %d
        if mod(d,10) == 0; disp(d); end
        
        % number of patients and controls at this density
        nc_t = length(ctrl_sig_id{d});
        np_t = length(schz_sig_id{d});
        
        % randomly select "nc" (and "np") matrices from combined set of r- and P-thresholded matrices
        % controls
        set1_c = sort(ctrl_sig_id{d}(randperm(nc_t,randi(nc_t-1))))'; % "-1" so that set 2 isn't null
        set2_c = sort(setdiff(ctrl_sig_id{d},set1_c))';
        
        % patients
        set1_p = sort(schz_sig_id{d}(randperm(np_t,randi(np_t-1))))'; % "-1" so that set 2 isn't null
        set2_p = sort(setdiff(schz_sig_id{d},set1_p))';
        
        % ctrl
        ctrl_r_overlap_perm(:,:,d) = sum(squeeze(logical(ctrl_mat_dens_r(:,:,d,set1_c))),3) + sum(squeeze(logical(ctrl_mat_dens(:,:,d,set2_c))),3);%/length(ctrl_sig_id{d});
        ctrl_p_overlap_perm(:,:,d) = sum(squeeze(logical(ctrl_mat_dens(:,:,d,set1_c))),3) + sum(squeeze(logical(ctrl_mat_dens_r(:,:,d,set2_c))),3);%/length(ctrl_sig_id{d});
        % schz
        schz_r_overlap_perm(:,:,d) = sum(squeeze(logical(schz_mat_dens_r(:,:,d,set1_p))),3) + sum(squeeze(logical(schz_mat_dens(:,:,d,set2_p))),3);%/length(schz_sig_id{d});
        schz_p_overlap_perm(:,:,d) = sum(squeeze(logical(schz_mat_dens(:,:,d,set1_p))),3) + sum(squeeze(logical(schz_mat_dens_r(:,:,d,set2_p))),3);%/length(schz_sig_id{d});
        
    end
    
    % frequencies of overlap
    for d = 1:1:ndens
        %d
        
        %%% ctrl
        % r
        temp_ctrl = ctrl_r_overlap_perm(:,:,d);
        tbl1_r = tabulate(temp_ctrl(triu_ind)); % obtain frequency of frequencies (how many edges appear twice, three times, four times, ...)
        ctrl_tbl_r_perm(d,tbl1_r(:,1)+1) = tbl1_r(:,2);
        % p
        temp_ctrl = ctrl_p_overlap_perm(:,:,d);
        tbl1_p = tabulate(temp_ctrl(triu_ind)); % obtain frequency of frequencies (how many edges appear twice, three times, four times, ...)
        ctrl_tbl_p_perm(d,tbl1_p(:,1)+1) = tbl1_p(:,2);
        %%% schz
        % r
        temp_ctrl = schz_r_overlap_perm(:,:,d);
        tbl1_r = tabulate(temp_ctrl(triu_ind)); % obtain frequency of frequencies (how many edges appear twice, three times, four times, ...)
        schz_tbl_r_perm(d,tbl1_r(:,1)+1) = tbl1_r(:,2);
        % p
        temp_ctrl = schz_p_overlap_perm(:,:,d);
        tbl1_p = tabulate(temp_ctrl(triu_ind)); % obtain frequency of frequencies (how many edges appear twice, three times, four times, ...)
        schz_tbl_p_perm(d,tbl1_p(:,1)+1) = tbl1_p(:,2);

    end
    
    % difference
    ctrl_tbl_d_perm(:,:,p) = ctrl_tbl_p_perm-ctrl_tbl_r_perm;
    schz_tbl_d_perm(:,:,p) = schz_tbl_p_perm-schz_tbl_r_perm;
    
end

save('pr_thr_consistency_perm_sig.mat','ctrl_tbl_d_perm','ctrl_tbl_d_perm')

% compare empirical difference to permuted differences

% ctrl
for d = 1:1:ndens 
    for i = 1:1:(length(ctrl_sig_id{d})+1)
        if ctrl_tbl_d(d,i) >= 0
            ctrl_tbl_d_p(d,i) = sum(ctrl_tbl_d_perm(d,i,:)>ctrl_tbl_d(d,i))/nperm;
        elseif ctrl_tbl_d(d,i) < 0
            ctrl_tbl_d_p(d,i) = sum(ctrl_tbl_d_perm(d,i,:)<ctrl_tbl_d(d,i))/nperm;
        end
    end
end
    
% schz
for d = 1:1:ndens
    for i = 1:1:(length(schz_sig_id{d})+1)
        if schz_tbl_d(d,i) >= 0
            schz_tbl_d_p(d,i) = sum(schz_tbl_d_perm(d,i,:)>schz_tbl_d(d,i))/nperm;
        elseif schz_tbl_d(d,i) < 0
            schz_tbl_d_p(d,i) = sum(schz_tbl_d_perm(d,i,:)<schz_tbl_d(d,i))/nperm;
        end
    end
end

save('pr_thr_consistency_perm_sig.mat','ctrl_tbl_d_p','ctrl_tbl_d_p')

