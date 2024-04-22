%% MOMP (Motif-Only Matrix Profile)
%% Inputs:  
            % T (Input time series)
            % m (Subsequence length) 
            % verbose (Output verbosity)
            
%% Output:
            % momp_out (motifs' distance)
            % momp_loc (motif locations)

%% Supporting functions
            % mpx_v2.m
            % SimMat.m (this will be updated in the next versions)
            % paa.m
            
%% Main Functions

function [momp_out, momp_loc] = momp(T, m, verbose)

    plotting = 0;
    run_mp = 1;
    
    T_orig = T;
    n_orig = length(T);
    dd = m / 4;
    bsf = inf; bsfloc = nan;
    indices = (1:n_orig)';
    
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
    
    tic;

    while true
           
        [uamp_base, uamp, absf, local_bsf_loc, ip] = upsample_approximate_mp(T, m, dd, indices);

        [bsf, bsfloc, local_bsf] = refine(T_orig, m, bsf, bsfloc, local_bsf_loc, dd);

        [pruned_T, pruned_indices] = prune(T_orig, m, indices, uamp, absf, bsf);
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
        
        dd = dd/2;
        T = pruned_T;
        indices = pruned_indices;
        
        if dd == 1
            [mp, ~, loc, ~] = mpx_v2(T, m, m); %for sanity check
            momp_time = toc;
            momp_out = min(mp); momp_loc = indices(loc(1:2,1));
            fprintf('MOMP : Tpaa1in%d | BSF: %0.2f {%d, %d} | Tot Time: %0.2f \n',...
                     dd, momp_out, momp_loc(1), momp_loc(2), momp_time);
            if run_mp
                fprintf('Speedup: %dX \n', floor(mp_time/momp_time));
            end
            
            break;
        end
       
    end

end



%%
function [uamp_base, uamp, absf, absf_loc, ip] = upsample_approximate_mp(T, m, dd, indices)
    %Current assumption : T = k1*dd m = k2*dd 
    tic;
    [ip] = ktip_v1(T, m, dd);
%     ktip_time = toc; fprintf('Tpaa1in%d  : ktip time: %0.2f \n', dd, ktip_time);
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
%     uamp = uamp_base - ip(1:length(uamp_base));
    absf = min(uamp);
    
end



function [bsf, bsfloc, localbsf] = refine(T, m, bsf, bsfloc, amloc, dd)

    ii = amloc(1); jj = amloc(2);
    localbsf = min(mpx_v2([T(ii:ii+m);T(jj: jj+m)], 4, m));
    st1 = max([1, ii - dd + 1]); end1 = min([length(T), ii + m + dd - 1]);
    st2 = max([1, jj - dd + 1, end1]); end2 = min([length(T), jj + m + dd - 1]);
    indices = [(st1:end1) , (st2:end2)];
    [mp, ~, loc, ~] = mpx_v2(T(indices), m, m);
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
    