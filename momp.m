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
            % paa.m
            
%% Main Functions

function [momp_out, momp_loc] = momp(T, m, verbose, run_mp, plotting)

    
    if nargin < 2
        error('incorrect number of input arguments');
    end


    if ~exist('verbose', 'var')
        verbose = 1;
    end
    if ~exist('run_mp', 'var')
        run_mp = 1;
    end
    if ~exist('plotting', 'var')
        plotting = 0;
    end
    



    T_orig = T;
    n_orig = length(T);
    dd = m / 16;
    bsf = inf; bsfloc = nan;
    indices = (1:n_orig)';
    profiling = 1;
    
    if verbose
        fprintf('T is length %d, and m is set to %d\n', n_orig, m);
        fprintf('Starting downsampling rate: %d\n', dd);
    end
    if run_mp
        tic;
        [mp, ~, loc, ~] = mpx_v2(T, m, m);
%         figure; plot(mp);
        mp_time = toc;
        fprintf('MP : %0.2f {%d, %d} | time: %0.2f\n', min(mp), loc(1,1), loc(2,1), mp_time);
    end
    
    ktip_tot_time = 0;
    uamp_tot_time = 0;
    refine_tot_time = 0;
    prune_tot_time = 0;
    curr_ktip_time = 0;
    curr_uamp_time = 0;
    curr_refine_time = 0;
    curr_prune_time = 0;
    
    momp_tot_tic = tic;

    while true
           
        [uamp_base, uamp, absf, local_bsf_loc, ip, ktip_time, uamp_time] = upsample_approximate_mp(T, m, dd, indices);
        
        refine_tic = tic;
        [bsf, bsfloc, local_bsf] = refine(T_orig, m, bsf, bsfloc, local_bsf_loc, dd);
        refine_time = toc(refine_tic);
        
        prune_tic = tic;
        [pruned_T, pruned_indices] = prune(T_orig, m, indices, uamp, absf, bsf);
        prune_time = toc(prune_tic);
        
        pruning = 1 - (length(pruned_T) / n_orig) ;
        
        % Summary
        fprintf('MOMP : Tpaa1in%d | BSF: %0.2f {%d, %d} | localBSF: %0.2f {%d, %d} | pruning: %0.4f\n',...
                  dd, bsf, bsfloc(1), bsfloc(2), local_bsf, local_bsf_loc(1), local_bsf_loc(2), pruning);
  
        if plotting
            
            figure;
            subplot(4,1,1); plot(T, 'Color', 'b', 'LineWidth', 1);
            title(sprintf('Tpaa1in%d Input T (No downsampling)', dd), 'FontSize',14)
            box off;
            subplot(4,1,2); plot(uamp_base, 'k', 'LineWidth', 1); 
            hold on; plot(ip, 'r', 'LineWidth', 1); 
            title('Initial UAMP', 'FontSize',14);  box off;
            legend('','KTIP');
            subplot(4,1,3); plot(uamp, 'b', 'LineWidth', 1); 
            title('Final UAMP', 'FontSize',14);  box off;
            hold on; yline(bsf, 'r', 'LineWidth', 1);
            legend('','BSF');
            subplot(4,1,4); plot(pruned_T, 'b', 'LineWidth', 1);
            title('Pruned T', 'FontSize',14); box off;
        end
        
        if profiling
            
            ktip_tot_time = ktip_tot_time + ktip_time;
            uamp_tot_time = uamp_tot_time + uamp_time;
            refine_tot_time = refine_tot_time + refine_time;
            prune_tot_time = prune_tot_time + prune_time;
            
            fprintf('Tpaa1in%d  : KTIP: %0.2fs \n', dd, ktip_time);
            fprintf('Tpaa1in%d  : UAMP: %0.2fs \n', dd, uamp_time);
            fprintf('Tpaa1in%d  : Refinement: %0.2fs \n', dd, refine_time);
            fprintf('Tpaa1in%d  : Prunning: %0.2fs \n', dd, prune_time);
        end
        
        dd = dd/2;
        T = pruned_T;
        indices = pruned_indices;
        
        if dd == 1
            [mp, ~, loc, ~] = mpx_v2(T, m, m); %for sanity check
            momp_out = min(mp); momp_loc = indices(loc(1:2,1));
            momp_time = toc(momp_tot_tic);
            fprintf('MOMP : Tpaa1in%d | BSF: %0.2f {%d, %d} | Tot Time: %0.2fs \n',...
                     dd, momp_out, momp_loc(1), momp_loc(2), momp_time);
            if run_mp
                fprintf('Speedup: %dX \n', floor(mp_time/momp_time));
            end
            
            if profiling
                fprintf('Tot KTIP time : %0.2fs  (%0.2f of momp time)\n',...
                                    ktip_tot_time, ktip_tot_time/momp_time);
                fprintf('Tot UAMP time : %0.2fs  (%0.2f of momp time)\n',...
                                    uamp_tot_time, uamp_tot_time/momp_time);
                fprintf('Tot Refine time : %0.2fs  (%0.2f of momp time)\n',...
                                    refine_tot_time, refine_tot_time/momp_time);  
                fprintf('Tot Prune time : %0.2fs  (%0.2f of momp time)\n',...
                                    prune_tot_time, prune_tot_time/momp_time);
            end
            
            break;
        end
       
    end

end



%%
function [uamp_base, uamp, absf, absf_loc, ip, ktip_time, uamp_time] = upsample_approximate_mp(T, m, dd, indices)
    %Current assumption : T = k1*dd m = k2*dd 
    ktip_tic = tic;
    [ip] = ktip_v1(T, m, dd);
    ktip_time = toc(ktip_tic); 
    
    uamp_tic = tic;
    n = length(T);
    [Tds, ~] = paa(T, floor(n/dd));
    mds = m/dd;
    [amp, ~, loc, ~] = mpx_v2(Tds, mds, mds);
    absf_loc = (loc(1:2,1) * dd) - dd + 1;
    absf_loc = indices(absf_loc);
    
    uamp_base = sqrt(dd) * repelem(amp, dd);
    uamp_base = uamp_base(1:length(T)-m+1);

    ipds = ip(1:dd:end);
    uip = repelem(ipds, dd);
    ip = uip;
    uamp = uamp_base - uip(1:length(uamp_base));
    absf = min(uamp);
    uamp_time = toc(uamp_tic); 
    
    
end



function [bsf, bsfloc, localbsf] = refine(T, m, bsf, bsfloc, amloc, dd)

    ii = amloc(1); jj = amloc(2);
    [localbsf] = MASS_s2(T(ii:ii+m), T(jj:jj+m));
    
    st1 = max([1, ii - dd + 1]); end1 = min([length(T), ii + m + dd - 1]);
    st2 = max([1, jj - dd + 1, end1]); end2 = min([length(T), jj + m + dd - 1]);
    indices = [(st1:end1) , (st2:end2)];
    
    [mp, ~, loc, ~] = mpx_v2(T(indices), (end1-st1), m);
    
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
    margin = m/4;
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
    