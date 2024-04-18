function momp(T, m, verbose)

    T_orig = T;
    n_orig = length(T);
    dd = m / 4;
    bsf = inf; bsfloc = nan;
    run_mp = 1;
    plotting = 1;
    indices = (1:n_orig)';
    
    if verbose
        fprintf('T is length %d, and m is set to %d\n', n_orig, m);
        fprintf('Starting downsampling rate: %d\n', dd);
    end
    if run_mp
        [mp, ~, loc, ~] = mpx_v2(T, m, m);
        fprintf('MP : %0.2f {%d, %d}\n', min(mp), loc(1,1), loc(2,1));
    end
    
    [ip] = KTIP(T,m, dd);

    while true
            
        [uamp_base, uamp, absf, local_bsf_loc] = upsample_approximate_mp(T, m, dd, ip);
         
            
        [bsf, bsfloc, local_bsf] = refine(T_orig, m, bsf, bsfloc, local_bsf_loc, dd);
        [pruned_T, pruned_indices] = prune(T_orig, m, indices, uamp, absf, bsf);
        pruning = 1 - (length(T) / n_orig) ;
        
        % Summary
        fprintf('MOMP : Tpaa1in%d | BSF: %0.2f {%d, %d} | localBSF: %0.2f {%d, %d} | pruning: %0.2f\n',...
                  dd, bsf, bsfloc(1), bsfloc(2), local_bsf, local_bsf_loc(1), local_bsf_loc(2), pruning);
  
        if plotting
            
            figure;
            subplot(3,1,1); plot(T, 'Color', 'b', 'LineWidth', 1);
            title(sprintf('Tpaa1in%d Input T (No downsampling)', dd), 'FontSize',14)
            box off;
            subplot(3,1,2); plot(uamp_base, 'Color', 'k', 'LineWidth', 1); 
            hold on; plot(ip, 'b', 'LineWidth', 1)
            hold on; plot(uamp, 'Color', 'r', 'LineWidth', 1)
            hold on; yline(bsf, 'g', 'LineWidth', 1);
            title('UAMP', 'FontSize',14)
            legend('initial UAMP', 'KTIP','Corrected UAMP', 'bsf');
            box off;
            subplot(3,1,3); plot(pruned_T, 'LineWidth', 1);
            title('Pruned T', 'FontSize',14); box off;
        end
        
        dd = dd/2;
        T = pruned_T;
        indices = pruned_indices;
        ip = ip(indices);
        
        if dd == 1
            [mp, ~, loc, ~] = mpx_v2(T, m, m); %for sanity check
            fprintf('MOMP : Tpaa1in%d | BSF: %0.2f {%d, %d} | localBSF: %0.2f {%d, %d}\n',...
                            dd, min(mp), loc(1,1), loc(2,1),min(mp), loc(1,1), loc(2,1));
            break;
        end
       
    end

end



%% Computing Upsampled Approximate Matrix Profile
function [ip] = KTIP(T,m, dd)
    ip = nan(length(T),1);
    
    for ii = 1:length(T)-m
        
        query = zscore(T(ii:ii+m));
        jj = min([ii+m+dd-1 , length(T)]);
        curr_T = zscore(T(ii+1:jj));
        [dist_profile] = SimMat(curr_T, m, query);
        dist_profile = dist_profile(1,:);
        ip(ii) = max(dist_profile);
        
    end
    
end


function [uamp_base, uamp, absf, absf_loc] = upsample_approximate_mp(T, m, dd, ip)
    %Current assumption : T = k1*dd m = k2*dd 
    n = length(T);
    [Tds, ~] = paa(T, floor(n/dd));
    mds = m/dd;
    [amp, ~, loc, ~] = mpx_v2(Tds, mds, mds);
    absf_loc = (loc(1:2,1) * 2) - 1;
    
    uamp_base = sqrt(dd) * repelem(amp, dd);
    uamp_base = uamp_base(1:length(T)-m+1);

    uamp = uamp_base - ip(1:length(uamp_base));
    
    absf = min(uamp);
    
end

function [bsf, bsfloc, localbsf] = refine(T, m, bsf, bsfloc, amloc, dd)

    ii = amloc(1); jj = amloc(2);
    localbsf = min(mpx([T(ii:ii+m);T(jj: jj+m)], 4, m));
    st1 = max([1, ii - dd]); end1 = min([length(T), ii + m + dd]);
    st2 = max([1, jj - dd, end1]); end2 = min([length(T), jj + m + dd]);
    indices = [(st1:end1) , (st2:end2)];
    [mp, ~, loc, ~] = mpx_v2(T(indices), m/2, m);
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
    margin = m;
    gaps = find([2; diff(targets)] > 1);
    pruned_indices = [];
    for ii=1:length(gaps)-1
        set = targets(gaps(ii):gaps(ii+1)-1);
        ss = max([1, set(1)-margin]);
        ee = min([set(end)+margin , indices(end)]);
        pruned_indices = [pruned_indices ; indices(ss:ee)];
    end
    if isempty(pruned_indices)
        ss = max([1 , targets(gaps(end))-margin]);
    else
        ss = max([pruned_indices(end), targets(gaps(end))-margin]);
    end
    ee = min([targets(end)+margin , indices(end)]);
    pruned_indices = [pruned_indices ; indices(ss:ee)];
end
    