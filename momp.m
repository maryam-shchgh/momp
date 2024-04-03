function [momp_out, mp_out] = momp(T, m, verbose)
    n_orig = length(T);
    if verbose
        fprintf('T is length %d, and m is set to %d\n', n_orig, m);
        [mp, ~, motifsIdx, ~] = mpx_v2(T, m, m);
        mp_out = motifsIdx(1:2,1);
    end
    
    %% downsamping the T and m
    dsr = m / 4;
    bsf = inf; refinedLoc = nan;
    while dsr > 2
        [amp, absf, bsf, prunedT] = process(T, m, dsr, bsf, refinedLoc, [1:length(T)]);
        figure(1);
        subplot(3,1,1)
        plot(T)
        subplot(3,1,2)
        plot(amp, 'b');
        hold on;
        plot(mp, 'g');
        yline(absf, 'r');
        yline(bsf, 'g');
        subplot(3,1,3)
        plot(prunedT)
        fprintf('T is length %d, and pruned is %d\n', length(T), length(prunedT));
        input('Enter to continue')
        T = prunedT;
    end
    momp_out = 0;
%     Tds = paa(T, floor(length(T)/dsr));
%     mds = floor(m/dsr);
%     
%     %% Computing the AMP and calculating the aBSF motifs
%     [amp, ~, motifsIdx, ~] = mpx_v2(Tds, mds, mds);
%     absfIdx = motifsIdx(1:2,1); % Best motif 
%     absf = amp(absfIdx(1));
%     amp = sqrt(dsr) * repelem(amp, dsr);
%     absf = sqrt(dsr) * absf;
%     absfIdx = dsr * absfIdx - 1;
%     
%     momp_out = 0;
%     %% BSF local search
%     bsf = bsfLocalSearch(T, m , dsr, absfIdx);    
%     %% Prunnig
%     filteredIdx = filterIndices(amp, absf, bsf, m, n_orig);
%     prunedT = T(filteredIdx);  
   
%     if verbose
%         figure(1);
%         subplot(3,1,1)
%         plot(T)
%         subplot(3,1,2)
%         plot(amp, 'b');
%         hold on;
%         plot(mp, 'g');
%         yline(absf, 'r');
%         yline(bsf, 'g');
%         subplot(3,1,3)
%         plot(prunedT)
%         fprintf('T is length %d, and pruned is %d\n', length(T), length(prunedT));
%     end
    

end


function [amp, absf, bsf, prunedT] = process(T, m, dsr, bsf, refinedLoc, idxList)

    Tds = paa(T, floor(length(T)/dsr));
    mds = floor(m/dsr);
    
    %% Computing the AMP and calculating the aBSF motifs
    [amp, ~, motifsIdx, ~] = mpx_v2(Tds, mds, mds);
    absfIdx = motifsIdx(1:2,1); % Best motif 
    absf = amp(absfIdx(1));
    amp = sqrt(dsr) * repelem(amp, dsr);
    absf = sqrt(dsr) * absf;
    absfIdx = dsr * absfIdx - 1;
    
    %% BSF local search
    [bsf, refindLoc] = bsfLocalSearch(T, m , dsr, absfIdx, bsf, refinedLoc, idxList);    
    %% Prunnig
    idxList = filterIndices(T, amp, absf, bsf, m);
    prunedT = T(idxList);  

end




function [bsf, refinedLoc] = bsfLocalSearch(T, m , dsr, absfIdx, bsf, refinedLoc, idxList)

    offset = dsr  - 1;
    range1 = (max(1,absfIdx(1) - offset) : min(length(T),absfIdx(1) + offset));
    range2 = (max(1,absfIdx(2) - offset) : min(length(T),absfIdx(2) + offset));
    for ii = range1
        sub1 = T(ii:ii+m);
        for jj = range2
            sub2 = T(jj:jj+m);
            [curr_bsf] = MASS_s2(sub1, sub2);
            if curr_bsf <= bsf
                bsf = curr_bsf;
                refinedLoc = [idxList(ii), idxList(jj)];
            end
        end
    end
    
%     width = 2*offset + m;
%     temp_T = zeros(1, width * length(absfIdx));
%     for ii = 1:length(absfIdx)
%         ll = absfIdx(ii);
%         temp_T((ii-1)*width + 1 : ii*width) = T(ll - offset + 1 : ll + offset + m);
%     end
%     [mp, ~] = mpx(temp_T, m, m);
%     bsf = min(mp);
end


function filteredIdx = filterIndices(T, amp, absf, bsf, m)

    indices = find(amp >= absf & amp <= bsf);
    % Generate indices for right neighbors
    candidates = indices + (1:m);

    % Clip indices to ensure they stay within bounds
    candidates = min(length(T) - 1, candidates);

    % Flatten the array and remove duplicates
    filteredIdx = unique(candidates(:));

    % Sort the result
    filteredIdx = sort(filteredIdx);


end