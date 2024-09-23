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
            % cat_timeseries.m
            
%% Main Functions

function [momp_out, momp_loc] = mist_v02(T_set, m, verbose, run_mpx)

    if nargin < 2
        error('incorrect number of input arguments');
    end


    if ~exist('verbose', 'var'), verbose = 1; end
    if ~exist('run_mpx', 'var'), run_mpx = 0; end
    plotting = 0;
    
    T_set_count = size(T_set, 2);

    dd = 64;
    global_bsf = inf; bsf = inf;
    global_bsf_T_index = nan;
    bsfloc = nan;
    Torders = 1:T_set_count;
    tot_times = zeros(1,4);

    run_mp_time = tic;
    if run_mpx
        for tt=1:T_set_count
            [mp, ~, loc, ~] = mpx_v2(T_set(1:end, tt), floor(m/2), m);
            fprintf('MPX - T%d : Motif = %0.2f {%d,%d}\n', tt, min(mp), loc(1), loc(2));
        end
    end
    run_mp_time = toc(run_mp_time);
    if ~run_mpx, run_mp_time = nan; end
    fprintf('Tot MPx time: %0.2f\n', run_mp_time);
      
    if verbose
        fprintf('Input set includes %d timeseries.\n', T_set_count);
        fprintf('Target motif length is set to %d.\n', m);
        fprintf('Initial downsampling rate is set to %d.\n', dd);
        for tt=Torders
            fprintf('>> T%d size: [%d X 1]\n', tt, length(T_set(:,tt)));
        end
    end

    if verbose, fprintf('====== MIST ======\n'); end
    momp_tot_tic = tic;

    ktip_tic = tic;
    [T_struct_arrary] = ktip_handler(T_set, Torders, m, dd, bsf, bsfloc);
    ktip_tot_time = toc(ktip_tic); 
    
    while dd >= 1
        for tt=Torders
            if all(isnan(T_set(1:end, tt)))
                continue;
            end
            if dd > 1
                curr_ktip = T_struct_arrary(tt).ktip; 
                curr_ktip = curr_ktip(1:end, log2(dd));
            else
                curr_ktip = zeros(T_struct_arrary(tt).ktip_rows,1); 
            end

            [T_struct_arrary(tt), tot_times] = momp_step(T_set(1:end, tt), T_struct_arrary(tt), m, dd, tt, global_bsf, curr_ktip, plotting, tot_times, momp_tot_tic);
            if isempty(T_struct_arrary(tt).T)
                T_set(1:end, tt) = nan(size(T_set(1:end, tt)));
                fprintf('# T%d : fully pruned!\n', tt)
            end
        
        end

        best_so_fars = [T_struct_arrary(Torders).bsf];

        % Sort the struct array based on the 'age' field
        [~, sortIdx] = sort(best_so_fars);
        
        Torders = Torders(sortIdx);
        global_bsf = T_struct_arrary(Torders(1)).bsf;
        global_bsf_T_index = Torders(1);
        dd = dd/2;
        fprintf('\n')
    end
    
    momp_out = global_bsf;
    momp_loc = T_struct_arrary(global_bsf_T_index).bsfloc;
    momp_time = toc(momp_tot_tic);
    final_output_and_plot(T_set, m, global_bsf_T_index, global_bsf_T_index, momp_out, momp_loc, ...
                                        ktip_tot_time, momp_time, tot_times, run_mp_time);
    
end




%%


function [T_struct, tot_times] = momp_step(T_orig, T_struct, m, dd, tt, gbsf, curr_ktip, plotting, tot_times,momp_tot_tic)

    [lbmp, ~, ~, ~, local_bsf_loc, ~, uamp_time, lbmp_time] = upsample_approximate_mp(T_struct.T, m, dd, T_struct.indices, curr_ktip, plotting);
    refine_tic = tic;
    [bsf, bsfloc] = refine(T_orig, m, T_struct.bsf, T_struct.bsfloc, local_bsf_loc, dd);
    refine_time = toc(refine_tic);
    
    
    prune_tic = tic;
    if dd > 1 
        [pruned_T, pruned_indices] = prune(T_orig, m, T_struct.indices, lbmp, min([bsf, gbsf]));
    else
        pruned_T = T_struct.T; pruned_indices = T_struct.indices;
    end
    prune_time = toc(prune_tic);

    pruning = 1 - (length(pruned_T) / length(T_orig)) ;

    uamp_tot_time = tot_times(1) + uamp_time;
    lbmp_tot_time = tot_times(2)+ lbmp_time;
    refine_tot_time = tot_times(3) + refine_time;
    prune_tot_time = tot_times(4) + prune_time;

    tot_times(1:4) = [uamp_tot_time, lbmp_tot_time, refine_tot_time, prune_tot_time];
    logInfo(length(T_struct.T), dd, tt, bsf, bsfloc,pruning, toc(momp_tot_tic), 0, uamp_time, lbmp_time,...
                                           refine_time, prune_time); 

    
    T_struct.T = pruned_T;
    T_struct.indices = pruned_indices;
    T_struct.bsf = bsf;
    T_struct.bsfloc = bsfloc;

end

function [T_struct_arrary] = ktip_handler(T_set, Torders, m, dd, bsf, bsf_loc)
    T_set_count = size(T_set, 2);
    T_struct_arrary(T_set_count) = struct('T', [],  'indices', [], 'dd', [], 'bsf', [], 'ktip', [], 'ktip_rows',[]);
    for tt=Torders
        [ktip] = ktip_v1(T_set(:,tt), m, dd);
        T_struct_arrary(tt).ktip = ktip;
        T_struct_arrary(tt).ktip_rows = size(ktip,1);
        T_struct_arrary(tt).T = T_set(:,tt);
        T_struct_arrary(tt).indices = (1:length(T_set(:,tt)))';
        T_struct_arrary(tt).dd = dd;
        T_struct_arrary(tt).bsf = bsf;
        T_struct_arrary(tt).bsfloc = bsf_loc;
    end
end


function final_output_and_plot(T_set, m, tt, gbsf_Ti, momp_out, momp_loc, ktip_tot_time, momp_time, tot_times, mpx_time)
    fprintf('====== Profiling Summary ======\n')
    fprintf('Tot KTIP time : %0.2fs  (%0.2f of momp time)\n',...
                        ktip_tot_time, ktip_tot_time/momp_time);
    fprintf('Tot UAMP time : %0.2fs  (%0.2f of momp time)\n',...
                        tot_times(1), tot_times(1)/momp_time);
    fprintf('Tot LBMP time : %0.2fs  (%0.2f of momp time)\n',...
                        tot_times(2), tot_times(2)/momp_time);
    fprintf('Tot Refine time : %0.2fs  (%0.2f of momp time)\n',...
                        tot_times(3), tot_times(3)/momp_time);  
    fprintf('Tot Prune time : %0.2fs  (%0.2f of momp time)\n',...
                        tot_times(4), tot_times(4)/momp_time);
    if ~isnan(mpx_time), fprintf('MIST over MPx speedup: %0.2fX\n', mpx_time/momp_time); end

    fprintf('Top-1 motif pair is %0.2f in timeseries %d\n', momp_out, gbsf_Ti);
    T_orig = T_set(1:end, gbsf_Ti);
    figure;
    subplot(2,1,1); plot(T_orig, 'b', 'LineWidth', 1); hold on;
    plot((momp_loc(1):momp_loc(1)+m), T_orig(momp_loc(1):momp_loc(1)+m), 'r', 'LineWidth', 1);
    plot((momp_loc(2):momp_loc(2)+m), T_orig(momp_loc(2):momp_loc(2)+m), 'r', 'LineWidth', 1);
    title(sprintf('Time Series %d', tt), 'FontSize',14); box off;
    subplot(2,1,2);
    plot(zscore(T_orig(momp_loc(1):momp_loc(1)+m)), 'b', 'LineWidth', 2);
    hold on;
    plot(zscore(T_orig(momp_loc(2):momp_loc(2)+m)), 'r', 'LineWidth', 1);
    title(sprintf('MOMP : T_{%d} and T_{%d} (znorm)', momp_loc(1), momp_loc(2)),'FontSize',14);
    box off;
end

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


function [lbmp, camp] = compLB(n, m, amp, ktip, dd, plotting)
    
    subsequence_count = n - m +1;
    ip = ktip(1:dd:subsequence_count); 
    n_ip = length(ip);
    amp = amp(1:n_ip);
    camp_ds = (sqrt(dd)*amp) - ip;
    
    if plotting
        camp = repelem(camp_ds, dd);
        camp = camp(1:n-m+1);
    else
        camp = nan;
    end
    
    lbmp_ds = -inf(size(ip)); 
    
    for ii=1:n_ip
        temp_lb = camp_ds - ip(ii);
        temp_lb(ii) = nan;
        lbmp_ds = max(temp_lb, lbmp_ds);
    end
    
    lbmp = repelem(lbmp_ds, dd);
    lbmp = lbmp(1:subsequence_count);
        
end


function [lbmp, camp, uamp, absf, absf_loc, uktip, uamp_time, lbmp_time] = upsample_approximate_mp(T, m, dd, indices, ip, plotting)
    
    %Current assumption : T = k1*dd m = k2*dd 
    mask = indices <= length(ip);
    ktip = ip(indices(mask));
    
    
    uamp_tic = tic;
    n = length(T);
    pad = (ceil(n/dd)*dd) - n;
    T = [T ; randn(pad,1)];
    
    if dd > 1, [Tds, ~] = paa(T, floor(n/dd)); else, Tds = T; end
    mds = floor(m/dd);

    [amp, ~, loc, ~] = mpx_v2(Tds, floor(mds/2), mds);
    absf_loc = (loc(1:2,1) * dd) - dd + 1;
    absf_loc = indices(absf_loc);
    
    
    if plotting
        uamp = sqrt(dd) * repelem(amp, dd);
        uamp = uamp(1:n-m+1);
        ktip_ds = ktip(1:dd:length(T)-m +1);
        uktip = repelem(ktip_ds, dd);
        uktip = uktip(1:n-m+1);
    else
        uamp = sqrt(dd) * repelem(amp, dd);
        uamp = uamp(1:n-m+1);
        uktip = nan;
        
    end
    uamp_time = toc(uamp_tic); 

    lbmp_tic = tic;
    [lbmp, camp] = compLB(n, m, amp, ktip, dd, plotting);
    
    absf = min(lbmp);
    lbmp_time = toc(lbmp_tic); 

    
end



function [bsf, bsfloc] = refine(T, m, bsf, bsfloc, amloc, dd)


    ii = amloc(1); jj = amloc(2);
    
   
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




function [pruned_T, pruned_indices] = prune(T, m, indices, uamp, bsf)
    if min(uamp) > bsf
        pruned_T = []; pruned_indices = []; 
    else
        targets = find(uamp <= bsf);
        [pruned_T, pruned_indices] = idxFilter(T, targets, indices, m);
    end
    
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


function logInfo(n, dd, tt, bsf, bsfloc, ...
                    pruning, curr_time, profile_steps, uamp_time, lbmp_time,...
                                                    refine_time, prune_time)

    fprintf('MOMP - T%d : Tpaa1in%d | BSF: %0.2f {%d, %d}| pruning: %0.4f | Time: %0.2f\n',...
           tt, dd, bsf, bsfloc(1), bsfloc(2), pruning, curr_time);
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
    