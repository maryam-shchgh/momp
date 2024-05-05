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
            % mpx_v2.m
            % ktip_v1.m
            % MASS_s2.m
            % paa.m
            
%% Main Functions

function [momp_out, momp_loc, scores] = momp_v4(T, m, verbose, run_mp, plotting)

    
    if nargin < 2
        error('incorrect number of input arguments');
    end


    if ~exist('verbose', 'var'), verbose = 1; end
    if ~exist('run_mp', 'var'), run_mp = 1; end
    if ~exist('plotting', 'var'), plotting = 0; end
    
    % Setting optional flags as needed inside the code
    profile_steps = 0; % To show detialed timing per MOMP step
    write_video = 0; % To build and save a video of pruning process (Not complete - Still working on it!!)
    
    
    
    if size(T,1) == 1, T = T'; end

    T_orig = T;
    n_orig = length(T);
    dd = (2^floor(log2(m))) / 64;
    
    % Computing initial bsf value by picking random indices from 1:n-m+1 
    while true
        [vals] = randi(n_orig-m+1, 2,1);
        if vals(2)-vals(1) > floor(m)
            break;
        end
    end

    [init_bsf] = MASS_s2(T(vals(1):vals(1)+m-1), T(vals(2):vals(2)+m-1)); 
    
    bsf = init_bsf; bsfloc = [vals(1),vals(2)];
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
        [mp, ~, loc, ~] = mpx_v2(T, floor(m/2), m);
        mp_time = toc;
        if verbose
            fprintf('====== MPx ======\n')
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
           
        [lbmp, camp, uamp, absf, local_bsf_loc, ip, uamp_time, lbmp_time] = upsample_approximate_mp(T, m, dd, indices, ktip(1:end, log2(dd)));
        
        refine_tic = tic;
        [bsf, bsfloc, local_bsf] = refine(T_orig, m, bsf, bsfloc, local_bsf_loc, dd);
        refine_time = toc(refine_tic);
        
        prune_tic = tic;
        [pruned_T, pruned_indices] = prune(T_orig, m, indices, lbmp, absf, bsf);
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
        
        if dd == 1
            [mp, ~, loc, ~] = mpx_v2(T, floor(m/2), m); %for sanity check
            momp_out = min(mp); momp_loc = indices(loc(1:2,1));
            momp_time = toc(momp_tot_tic);
            fprintf('MOMP : Tpaa1in%d | BSF: %0.2f {%d, %d} | Tot Time: %0.2fs \n',...
                     dd, momp_out, momp_loc(1), momp_loc(2), momp_time);
                 
            bsf_vals(length(bsf_vals),:) = [momp_out, momp_time];
            [scores] = score_momp(bsf_vals);
            
            figure;
            plot(T_orig(momp_loc(1):momp_loc(1)+m), 'b', 'LineWidth', 1);
            hold on;
            plot(T_orig(momp_loc(2):momp_loc(2)+m), 'r', 'LineWidth', 1);
            title(sprintf('MOMP : T_{%d} and T_{%d}', momp_loc(1), momp_loc(2)),'FontSize',14);
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


function [lbmp, camp] = compLB(T, m, amp, ktip, dd)
    
    subsequence_count = length(T)-m +1;
    ktip_ds = ktip(1:dd:subsequence_count);   
    camp_ds = (sqrt(dd)*amp) - ktip_ds;
    camp = repelem(camp_ds, dd);
    camp = camp(1:length(T)-m+1);
    lbmp_ds = -inf(size(ktip_ds));
    
    for ii=1:length(amp)
        ktip_ds_ii = ones(size(ktip_ds))*ktip_ds(ii);   
        ktip_ds_ii(ii) = inf;
        temp_lb = camp_ds - ktip_ds_ii;
        lbmp_ds = max(temp_lb, lbmp_ds);
    end
    
    lbmp = repelem(lbmp_ds, dd);
    lbmp = lbmp(1:subsequence_count);
        
end


function [lbmp, camp, uamp, absf, absf_loc, uktip, uamp_time, lbmp_time] = upsample_approximate_mp(T, m, dd, indices, ip)
    
    %Current assumption : T = k1*dd m = k2*dd 
    mask = indices <= length(ip);
    ktip = ip(indices(mask));
    ktip_ds = ktip(1:dd:length(T)-m +1);
    
    uamp_tic = tic;
    n = length(T);
    [Tds, ~] = paa(T, floor(n/dd));
    mds = floor(m/dd);
    [amp, ~, loc, ~] = mpx_v2(Tds, floor(mds/2), mds);
    absf_loc = (loc(1:2,1) * dd) - dd + 1;
    absf_loc = indices(absf_loc);
    
    
    uamp = sqrt(dd) * repelem(amp, dd);
    uamp = uamp(1:length(T)-m+1);
    uamp_time = toc(uamp_tic); 

    
    uktip = repelem(ktip_ds, dd);
    uktip = uktip(1:length(T)-m+1);

    lbmp_tic = tic;
    [lbmp, camp] = compLB(T, m, amp, ktip, dd);
    
    absf = min(lbmp);
    lbmp_time = toc(lbmp_tic); 

    
end



function [bsf, bsfloc, localbsf] = refine(T, m, bsf, bsfloc, amloc, dd)

    ii = amloc(1); jj = amloc(2);
    [localbsf] = MASS_s2(T(ii:ii+m-1), T(jj:jj+m-1));
    
   
    st1 = max([1, ii - dd + 1]); end1 = ii + m + dd - 1;
    st2 = max([jj - dd + 1, end1]); end2 = min([length(T), jj + m + dd - 1]);


    indices = [(st1:end1), (st2:end2)];
    T_temp = [T(st1:end1); T(st2:end2)];
    
    [mp, ~, loc, ~] = mpx_v2(T_temp, floor(m/2), m);
    
    currbsf = min(mp);

    if currbsf <= bsf
        bsf = currbsf;
        bsfloc = indices(loc(1:2,1));

    end
   
end


function [pruned_T, pruned_indices] = prune(T, m, indices, uamp, absf, bsf)

    targets = find(uamp >= absf & uamp <= bsf);
    [pruned_indices] = idxFilter(targets, indices, m);
    pruned_T = T(pruned_indices);
    
end

function [pruned_indices] = idxFilter(targets, indices, m)
    margin = floor(m/4);
    gaps = find([2; diff(targets)] > 1);
    pruned = [];
    for ii=1:length(gaps)-1
        set = targets(gaps(ii):gaps(ii+1)-1);
        if isempty(pruned)
            ss = max([1, set(1)- margin]);
        else
            ss = max([1, set(1)- margin, pruned(end)+1]);
        end
        ee = min([set(end)+m + margin , length(indices)]);
        pruned = [pruned ; (ss:ee)'];
    end
    if isempty(pruned)
        ss = max([1 , targets(gaps(end))-margin]);
    else
        ss = max([pruned(end) + 1, targets(gaps(end))-margin]);
    end
   
    ee = min([targets(end)+ m + margin , length(indices)]);

    pruned = [pruned ; (ss:ee)'];
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
    