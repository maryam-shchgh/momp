%% MOMP (Motif-Only Matrix Profile)

%% Inputs:  
            % T (Input time series)
            % m (Subsequence length) 
            % [optional] verbose (Output verbosity - default = 1)
            % [optional] run_mp (run mpx for comparison - default = 1)
            % [optional] plotting (Output plots - default = 0)
            
%% Output:
            % momp_out (motifs' distance)
            % momp_loc (motif locations)

%% Supporting files/functions
            % mpx_v2_AB.m
            % ktip_v1.m
            % MASS_s2.m
            % paa.m
            % cat_timeseries.m
            
%% Main Functions

function [momp_out, momp_loc, scores] = momp_v5_AB(T, m, split, verbose, run_mp, plotting)

    
    if nargin < 3
        error('incorrect number of input arguments');
    end


    if ~exist('verbose', 'var'), verbose = 1; end
    if ~exist('run_mp', 'var'), run_mp = 0; end
    if ~exist('plotting', 'var'), plotting = 0; end
    
    % Setting optional flags as needed inside the code
    profile_steps = 0; % To show detialed timing per MOMP step
    write_video = 0; % To build and save a video of pruning process
    save_terminal_log = 0;
    
    if save_terminal_log, log_filename = ('terminal-log.txt'); diary(log_filename); end
    
    if size(T,1) == 1, T = T'; end

    T_orig = T;
    n_orig = length(T);
    if m < 512, dd = (2^floor(log2(m))) / 4; else, dd = (2^floor(log2(m))) / 64; end
    
    % Computing initial bsf value by picking random indices from 1:n-m+1 
    [val_1] = randi(split-m+1, 1,1);
    [val_2] = randi((n_orig-split)-m+1, 1,1);
    val_2 = val_2 + split - 1;

    [init_bsf] = MASS_s2(zscore(T(val_2:val_2+m-1)), zscore(T(val_1:val_1+m-1)));
    
    bsf = init_bsf; bsfloc = [val_1,val_2];
    bsf_vals = zeros(log2(dd) + 2, 2);
    bsf_vals(1,:) = [init_bsf, 0];
    
    indices = (1:n_orig)';
    
    if write_video
        v = VideoWriter(sprintf('T_%d_m_%d_initDsr_%d', n_orig, m, dd), 'MPEG-4');  % Name your video file
        v.FrameRate = 2;  % Set frame rate (frames per second)
        open(v);
    end
    
    if verbose
        fprintf('T is length %d, and m is set to %d\n', n_orig, m);
        fprintf('Starting downsampling rate: %d\n', dd);
    end
    if run_mp
        tic;
        [mp, ~, loc] = mpx_v2_AB(T, m, split);
        plot(mp)
        mp_time = toc;
        if verbose
            fprintf('====== MPx (AB) ======\n')
            fprintf('MP : %0.2f {%d, %d} | time: %0.2f\n', min(mp), loc(1,1), loc(2,1), mp_time);
        end
    end
    
    uamp_tot_time = 0;
    lbmp_tot_time = 0;
    refine_tot_time = 0;
    prune_tot_time = 0;
    
    momp_tot_tic = tic;
    
    ktip_tic = tic;
    [ktip] = ktip_v1(T, m, dd);
    ktip_tot_time = toc(ktip_tic); 
    
    if verbose, fprintf('====== MOMP ======\n'); end

    while true
        [lbmp, camp, uamp, absf, local_bsf_loc, ip, uamp_time, lbmp_time] = upsample_approximate_mp(T, m, split, dd, indices, ktip(1:end, log2(dd)));

        refine_tic = tic;
        [bsf, bsfloc, local_bsf] = refine(T_orig, m, bsf, bsfloc, local_bsf_loc, dd);
        refine_time = toc(refine_tic);

        prune_tic = tic;
        [pruned_T, pruned_indices, pruned_split] = prune(T_orig, m, indices, split, lbmp, absf, bsf);
        prune_time = toc(prune_tic);

        pruning = 1 - (length(pruned_T) / n_orig) ;
  
        if plotting , generatePlots(T, dd, uamp, ip, camp, lbmp, bsf, pruned_T); end
        
        if write_video , pruning_flow(v, T_orig, pruned_indices); end
        
        uamp_tot_time = uamp_tot_time + uamp_time;
        lbmp_tot_time = lbmp_tot_time + lbmp_time;
        refine_tot_time = refine_tot_time + refine_time;
        prune_tot_time = prune_tot_time + prune_time;
        
        
        bsf_vals(length(bsf_vals)-log2(dd), :) = [bsf, toc(momp_tot_tic)];

        if verbose, logInfo(length(T), dd, bsf, bsfloc, local_bsf,local_bsf_loc,pruning, toc(momp_tot_tic), profile_steps, uamp_time, lbmp_time,...
                                               refine_time, prune_time); end


        dd = dd/2;
        T = pruned_T;
        indices = pruned_indices;
        split = pruned_split;
        
        if dd == 1
            curr_split = find(indices == split);
            [mp, ~, loc] = mpx_v2_AB(T, m, curr_split); %for sanity check
            momp_out = min(mp); momp_loc = indices(loc(1:2,1));
            momp_time = toc(momp_tot_tic);
            fprintf('MOMP : Tpaa1in%d | BSF: %0.2f {%d, %d} | Tot Time: %0.2fs \n',...
                     dd, momp_out, momp_loc(1), momp_loc(2), momp_time);
                 
            bsf_vals(length(bsf_vals),:) = [momp_out, momp_time];
            [scores] = score_momp(bsf_vals);
            
            figure;
            subplot(2,1,1); plot(T_orig, 'b', 'LineWidth', 1); hold on;
            plot((momp_loc(1):momp_loc(1)+m), T_orig(momp_loc(1):momp_loc(1)+m), 'r', 'LineWidth', 1);
            plot((momp_loc(2):momp_loc(2)+m), T_orig(momp_loc(2):momp_loc(2)+m), 'r', 'LineWidth', 1);
            title('Input Time Series', 'FontSize',14); box off;
            subplot(2,1,2);
            plot(zscore(T_orig(momp_loc(1):momp_loc(1)+m)), 'b', 'LineWidth', 2);
            hold on;
            plot(zscore(T_orig(momp_loc(2):momp_loc(2)+m)), 'r', 'LineWidth', 1);
            title(sprintf('MOMP : T_{%d} and T_{%d} (znorm)', momp_loc(1), momp_loc(2)),'FontSize',14);
            box off;
            
            
            if verbose
                fprintf('====== Profiling Summary ======\n')
                fprintf('Tot KTIP time : %0.2fs  (%0.2f of momp time)\n',...
                                    ktip_tot_time, ktip_tot_time/momp_time);
                fprintf('Tot UAMP time : %0.2fs  (%0.2f of momp time)\n',...
                                    uamp_tot_time, uamp_tot_time/momp_time);
                fprintf('Tot LBMP time : %0.2fs  (%0.2f of momp time)\n',...
                                    lbmp_tot_time, lbmp_tot_time/momp_time);
                fprintf('Tot Refine time : %0.2fs  (%0.2f of momp time)\n',...
                                    refine_tot_time, refine_tot_time/momp_time);  
                fprintf('Tot Prune time : %0.2fs  (%0.2f of momp time)\n',...
                                    prune_tot_time, prune_tot_time/momp_time);
                                
                if run_mp, fprintf('>> Speedup: %dX \n', floor(mp_time/momp_time)); end
            end
            
            if write_video, close(v); end
            if save_terminal_log, diary off; end
            
            break;
        end
       
    end

end




%%
function [scores] = score_momp(bsf_vals)
    scores = zeros(length(bsf_vals),2);
    n = length(bsf_vals);
    bsf_init = bsf_vals(1,1);
    bsf_final = bsf_vals(n,1);
    
    for ii=1:length(bsf_vals)
        bsf_curr = bsf_vals(ii,1);
        time_curr = bsf_vals(ii,2);
        scores(ii,:) = [time_curr, 1- ((bsf_curr - bsf_final)/bsf_init)];
    end 
end


function [lbmp, camp] = compLB(n, m, amp, ktip, dd, split_ds)
    
    subsequence_count = n - m +1;
    ip = ktip(1:dd:subsequence_count); 
    n_ip = length(ip);
    amp = amp(1:n_ip);
    camp_ds = (sqrt(dd)*amp) - ip;
    camp = repelem(camp_ds, dd);
    camp = camp(1:n-m+1);
    lbmp_ds = -inf(size(ip)); 
    temp_lb = lbmp_ds;
    
    for ii=1:n_ip
        if ii < split_ds
            temp_lb(split_ds:end) = camp_ds(split_ds:end) - ip(ii);
        else
            temp_lb(1:split_ds-1) = camp_ds(1:split_ds-1) - ip(ii);
        end
        lbmp_ds = max(temp_lb, lbmp_ds);
    end
    
    lbmp = repelem(lbmp_ds, dd);
    lbmp = lbmp(1:subsequence_count);
        
end


function [lbmp, camp, uamp, absf, absf_loc, uktip, uamp_time, lbmp_time] = upsample_approximate_mp(T, m, split, dd, indices, ip)
    
    %Current assumption : T = k1*dd m = k2*dd 
    mask = indices <= length(ip);
    ktip = ip(indices(mask));
    ktip_ds = ktip(1:dd:length(T)-m +1);
    
    uamp_tic = tic;
    n = length(T);
    curr_split = find(indices == split);
    T_head = T(1:curr_split-1); T_tail = T(curr_split:end);
    pad_head = (ceil(length(T_head)/dd)*dd) - length(T_head);
    pad_tail = (ceil(length(T_tail)/dd)*dd) - length(T_tail) ;
    T_head = [T_head ; randn(pad_head,1)];
    T_tail = [T_tail ; randn(pad_tail,1)];
    [Tds_head, ~] = paa(T_head, floor(length(T_head)/dd));
    [Tds_tail, ~] = paa(T_tail, floor(length(T_tail)/dd));
        
    Tds = [Tds_head; Tds_tail];
    
    mds = floor(m/dd);
    [amp, ~, loc] = mpx_v2_AB(Tds, mds, length(Tds_head)+1);
    absf_loc = (loc(1:2,1) * dd) - dd + 1;
    absf_loc = indices(absf_loc);
    
    uamp = sqrt(dd) * repelem(amp, dd);
    uamp = uamp(1:n-m+1);
    uamp_time = toc(uamp_tic); 
    
    uktip = repelem(ktip_ds, dd);
    uktip = uktip(1:n-m+1);

    lbmp_tic = tic;
    [lbmp, camp] = compLB(n, m, amp, ktip, dd, length(Tds_head)+1);
    
    absf = min(lbmp);
    lbmp_time = toc(lbmp_tic); 

    
end



function [bsf, bsfloc, localbsf] = refine(T, m, bsf, bsfloc, amloc, dd)

    ii = amloc(1); jj = amloc(2);
    pad_size = nan;
    if (jj+m-1) > length(T), pad_size = (jj+m-1) - length(T); T = [T ; randn(pad_size,1)]; end
    
    [localbsf] = MASS_s2(zscore(T(ii:ii+m-1)), zscore(T(jj:jj+m-1)));
    if ~isnan(pad_size), T = T(1:end-pad_size); end
    
   
    st1 = max([1, ii - dd + 1]); end1 = ii + m + dd - 1;
    st2 = max([jj - dd + 1, end1]); end2 = min([length(T), jj + m + dd - 1]);


    indices = [(st1:end1), (st2:end2)];
    T_temp = [T(st1:end1); T(st2:end2)];
    
    [mp, ~, loc] = mpx_v2_AB(T_temp, m, length(T(st1:end1))+1);
    
    currbsf = min(mp);

    if currbsf <= bsf
        bsf = currbsf;
        bsfloc = indices(loc(1:2,1));

    end
   
end

function [aligned_T, split] = cat_T(T1, T2, dd)

    if isrow(T1), T1 = T1'; end
    if isrow(T2), T2 = T2'; end

    if isempty(T1), aligned_T = T2; split = nan;
    elseif isempty(T2), aligned_T = T1; split = nan;

    else

        if exist('dd','var') && dd > 1
            n = length(T1);
            pad = (ceil(n/dd)*dd) - n;
            T1 = [T1 ; randn(pad,1)];
            [T1, ~] = paa(T1, floor(n/dd));
            n = length(T2);
            pad = (ceil(n/dd)*dd) - n;
            T2 = [T2 ; randn(pad,1)];
            [T2, ~] = paa(T2, floor(n/dd));
        end

        startpoint = T2(1);
        endpoint= T1(end);

        shiftval = endpoint - startpoint;
        aligned_T = [T1 ; T2 + shiftval];
        split = length(T1)+1;
    end

end


function [pruned_T, pruned_indices, pruned_split] = prune(T, m, indices, split, uamp, absf, bsf)

    targets = find(uamp >= absf & uamp <= bsf);
    [pruned_T, pruned_indices] = idxFilter(T, targets, indices, m);
%     pruned_T = T(pruned_indices);
    pruned_split = min(pruned_indices(pruned_indices >= split));
    
end

function [pruned_T, pruned_indices] = idxFilter(T, targets, indices, m)
    margin = floor(m/4);
    gaps = find([2; diff(targets)] > 1);
    pruned = [];
    pruned_T = [];
    for ii=1:length(gaps)-1
        set = targets(gaps(ii):gaps(ii+1)-1);
        if isempty(pruned)
            ss = max([1, set(1)- margin]);
        else
            ss = max([1, set(1)- margin, pruned(end)+1]);
        end
        ee = min([set(end)+m + margin , length(indices)]);
        pruned = [pruned ; (ss:ee)'];
        [pruned_T, ~] = cat_T(pruned_T, T(indices(ss:ee)));
    end
    if isempty(pruned)
        ss = max([1 , targets(gaps(end))-margin]);
    else
        ss = max([pruned(end) + 1, targets(gaps(end))-margin]);
    end
   
    ee = min([targets(end)+ m + margin , length(indices)]);

    pruned = [pruned ; (ss:ee)'];
    [pruned_T, ~] = cat_T(pruned_T, T(indices(ss:ee)));
    pruned_indices = indices(pruned);

    
end


function logInfo(n, dd, bsf, bsfloc, local_bsf, local_bsf_loc,...
                    pruning, curr_time, profile_steps, uamp_time, lbmp_time,...
                                                    refine_time, prune_time)

    fprintf('MOMP : Tpaa1in%d | BSF: %0.2f {%d, %d} | localBSF: %0.2f {%d, %d} | pruning: %0.4f | Time: %0.2f\n',...
           dd, bsf, bsfloc(1), bsfloc(2), local_bsf, local_bsf_loc(1), local_bsf_loc(2), pruning, curr_time);
    if profile_steps
        fprintf('Tpaa1in%d - len(T) : %d \n', dd, n);
        fprintf('Tpaa1in%d  : UAMP: %0.2fs \n', dd, uamp_time);
        fprintf('Tpaa1in%d  : LBMP: %0.2fs \n', dd, lbmp_time);
        fprintf('Tpaa1in%d  : Refinement: %0.2fs \n', dd, refine_time);
        fprintf('Tpaa1in%d  : Prunning: %0.2fs \n', dd, prune_time);
    end
end


function generatePlots(T, dd, uamp, ip, camp, lbmp, bsf, pruned_T)

    figure;
    subplot(3,1,1); plot(T, 'Color', 'b', 'LineWidth', 1);
    title(sprintf('Tpaa1in%d Input T (No downsampling)', dd), 'FontSize',14)
    box off;

    subplot(3,1,2); plot(uamp, 'k', 'LineWidth', 1); 
    hold on; plot(ip, 'g', 'LineWidth', 1); 
    title('Initial UAMP', 'FontSize',14);  box off;
    hold on; plot(camp, 'b', 'LineWidth', 1); 
    hold on; plot(lbmp, 'Color', [0.3010 0.7450 0.9330], 'LineWidth', 1); 
    hold on; yline(bsf, 'Color', 'r', 'LineWidth', 1);
    legend('AMP', 'KTIP', 'CAMP', 'LBMP', 'BSF');
    title('AMP --> CAMP --> LBMP', 'FontSize',14);  box off;

    subplot(3,1,3); plot(pruned_T, 'b', 'LineWidth', 1);
    title('Pruned T', 'FontSize',14); box off;
end


function pruning_flow(v, T_orig, pruned_indices)
    video_fig = figure('Visible', 'off');
    set(video_fig, 'OuterPosition', [100, 100, 800, 450]);
    n = length(T_orig);
    mask = true(size(T_orig));
    mask(pruned_indices) = false;
    pruned_T_plot = T_orig;
    pruned_T_plot(mask) = nan;
    if isnan(pruned_T_plot(1)), pruned_T_plot(1) = 0; end
    if isnan(pruned_T_plot(n)), pruned_T_plot(n) = 0; end
    plot((1:n), pruned_T_plot, 'b', 'LineWidth', 1);
    ylim([0, max(T_orig)]);
    set(gca, 'YTick', [], 'YTickLabel', []);
    title(sprintf('n = %d | pruning = %0.2f', n, 1-(length(pruned_indices)/n)), 'FontSize',14); box off;
    frame = getframe(video_fig);
    writeVideo(v, frame);
    close(video_fig);
end
    