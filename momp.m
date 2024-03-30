function [momp_out, mp_out] = momp(T, m, verbose)

    if verbose
        fprintf('T is length %d, and m is set to %d\n', length(T), m);
        [mp, ~, motifsIdx, ~] = mpx_v2(T, m, m);
        mp_out = motifsIdx(1:2,1);
    end
    
    %% downsamping the T and m
    dsr = 4;
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
    bsf = bsfLocalSearch(T, m , dsr, absfIdx);
    
    %% Prunnig
    filteredIdx = filterIndices(amp, absf, bsf, m);
    prunedT = T(filteredIdx);
    [prunedmp, ~, motifsIdx, ~] = mpx_v2(prunedT, m, m);
    momp_out = filteredIdx(motifsIdx(1:2,1));
    figure(2); 
    subplot(2,1,1); plot(T);
    subplot(2,1,2); plot(prunedT);hold on, plot(prunedmp);
    
    
    if verbose
        figure(1);
        plot(amp, 'b');
        hold on;
        plot(mp, 'g');
        yline(absf, 'r');
        yline(bsf, 'g');
    end
    
    
    

    
end


function bsf = bsfLocalSearch(T, m , dsr, absfIdx)

    offset = dsr  - 1;
    width = 2*offset + m;
    temp_T = zeros(1, width * length(absfIdx));
    for ii = 1:length(absfIdx)
        ll = absfIdx(ii);
        temp_T((ii-1)*width + 1 : ii*width) = T(ll - offset + 1 : ll + offset + m);
    end
    [mp, ~] = mpx(temp_T, m, m);
    bsf = min(mp);
end


function filteredIdx = filterIndices(amp, absf, bsf, m)

    firstCandidates = find(amp >= absf & amp <= bsf);
    filteredIdx = [];
    
    for ii = 1:length(firstCandidates)
        first = max(1, firstCandidates(ii) - m);
        last = min(firstCandidates(ii) + m, length(amp));
        filteredIdx(end:end+last-first + 1) = [filteredIdx, first:last];
    end

end