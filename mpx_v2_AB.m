function [matrixProfile, matrixProfileIdx, motifsIdx] = mpx_v2_AB(timeSeries, minlag, subseqlen, split)


subcount = length(timeSeries) - subseqlen + 1;

if nargin ~= 4
    error('incorrect number of input arguments');
elseif ~isvector(timeSeries)
    error('first argument must be a 1D vector');
elseif ~(isfinite(subseqlen) && floor(subseqlen) == subseqlen) || (subseqlen < 2) || (subcount < 2) 
    error('subsequence length must be an integer value between 2 and the length of the timeSeries');
end

transposed_ = isrow(timeSeries);
if transposed_
    timeSeries = transpose(timeSeries);
end


nanmap = find(~isfinite(movsum(timeSeries, [0 subseqlen-1], 'Endpoints', 'discard')));
nanIdx = find(isnan(timeSeries));
timeSeries(nanIdx) = 0;

% We need to remove any NaN or inf values before computing the moving mean,
% because it uses an accumulation based method. We add additional handling
% elsewhere as needed.


mu = moving_mean(timeSeries, subseqlen);
mus = moving_mean(timeSeries, subseqlen - 1);
invnorm = zeros(subcount, 1);


% movstd wasn't sufficient to prevent errors around constant or nearly constant regions
% This is safer, since we can reliably skip them
for i = 1 : subcount
    invnorm(i) = 1 ./ norm(timeSeries(i : i + subseqlen - 1) - mu(i));
end

invnorm(nanmap) = NaN;
invnorm(~isfinite(invnorm)) = NaN;


% This method requires more arrays, and it seems to be marginally worse
% at times. It isn't impacted as much by cancellation though.
dr_bwd = [0; timeSeries(1 : subcount - 1) - mu(1 : subcount - 1)];
dc_bwd = [0; timeSeries(1 : subcount - 1) - mus(2 : subcount)];
dr_fwd = timeSeries(subseqlen : end) - mu(1 : subcount);
dc_fwd = timeSeries(subseqlen : end) - mus(1 : subcount);
matrixProfile = repmat(-1, subcount, 1);
matrixProfile(nanmap) = NaN;
matrixProfileIdx = NaN(subcount, 1);


for diag = 2 : subcount
    cov_ = sum((timeSeries(diag : diag + subseqlen - 1) - mu(diag)) .* (timeSeries(1 : subseqlen) - mu(1)));
    for row = max(1, split-diag+1) : min(split-1, subcount - diag + 1)
        col = diag + row - 1;
        if row > 1
            cov_ = cov_ - dr_bwd(row) * dc_bwd(col) + dr_fwd(row) * dc_fwd(col);
        end
        corr_ = cov_ * invnorm(row) * invnorm(col);
        if corr_ > matrixProfile(row)
            matrixProfile(row) = corr_;
            matrixProfileIdx(row) = col;
        end
        if corr_ > matrixProfile(col)
            matrixProfile(col) = corr_;
            matrixProfileIdx(col) = row;
        end
    end
end

[motifsIdx] = findMotifs(timeSeries, mu, invnorm, matrixProfile, matrixProfileIdx, subseqlen, minlag);

matrixProfile = sqrt(max(0, 2 * subseqlen * (1 - matrixProfile), 'includenan'));
% We don't make a copy earlier or at this point, because this usually does nothing.
timeSeries(nanIdx) = NaN;


if transposed_  % matches the profile and profile index but not the motif or discord index to the input format
    matrixProfile = transpose(matrixProfile);
    matrixProfileIdx = transpose(matrixProfileIdx);
end

end


function [motifIdxs, matrixProfile] = findMotifs(timeSeries, mu, invnorm, matrixProfile, profileIndex, subseqLen, exclusionLen)
% This is adapted match the output of some inline code written by Michael Yeh
% to find the top k motifs in a time series. Due to some bug fixes I applied, the two
% may return slightly different results, although they almost always agree on the top 2 motif pairs. 


% +2 allows us to pack the original motif pair with its neighbors.
% Annoyingly if it's extended, matlab will pad with zeros, leading to some
% rather interesting issues later.
motifCount = 3;
radius = 2;
neighborCount = 10;

motifIdxs = NaN(neighborCount + 2, motifCount);
% This formula just happens to work pretty well in picking time vs frequency based
% methods for optimal execution time. Feel free to change. It's a little
% skewed in favor of ffts, because in the current version of Matlab (2019),
% the implementation of conv2 does not to take advantage of simd
% optimization, based on my informal timing tests and various comparisons.
% 
if subseqLen < maxNumCompThreads * 128
   crosscov = @(idx) conv2(timeSeries, (timeSeries(idx + subseqLen - 1 : -1 : idx) - mu(idx)) .* invnorm(idx), 'valid');
else
   padLen = 2^nextpow2(length(timeSeries));
   crosscov = @(idx) ifft(fft(timeSeries, padLen) .* conj(fft((timeSeries(idx : idx + subseqLen - 1) - mu(idx)) .* invnorm(idx), padLen)), 'symmetric');
end

for i = 1 : motifCount
    [corr_, motIdx] = max(matrixProfile);
    % -1 means this is maximally negatively correlated and was therefore
    % never updated
    if ~isfinite(corr_) || corr_ == -1 
        break;
    end
    % order subsequence motif pair as [time series index of 1st appearance, time series index of 2nd appearance]
    motifIdxs(1 : 2, i) = [min(motIdx, profileIndex(motIdx)), max(motIdx, profileIndex(motIdx))];
    [corrProfile] = crosscov(motIdx);
    corrProfile = min(1, corrProfile(1 : length(timeSeries) - subseqLen + 1) .* invnorm, 'includenan');
    % This uses correlation instead of normalized Euclidean distance, because it's easier to work with and involves fewer operations.
    
    corrProfile(isnan(matrixProfile)) = NaN;
    if exclusionLen > 0
        for j = 1 : 2
            exclRangeBegin = max(1, motifIdxs(j, i) - exclusionLen + 1); 
            exclRangeEnd = min(length(matrixProfile), motifIdxs(j, i) + exclusionLen - 1);
            corrProfile(exclRangeBegin : exclRangeEnd) = NaN;
        end
    end
    
    for j = 3 : neighborCount + 2
        [neighborCorr, neighbor] = max(corrProfile);
        % If you want to put a reasonable global bound on it, 
        % set the if statement to also skip anything where neighborCorr <= 0 as well.
        % This was reverted to match earlier code, which did not enforce such a restriction.
        % 
        if  ~isfinite(neighborCorr) || ((1 - neighborCorr) >= radius * (1 - corr_)) 
            break;
        end
        motifIdxs(j, i) = neighbor;
        if exclusionLen > 0
            exclRangeBegin = max(1, neighbor - exclusionLen + 1); 
            exclRangeEnd = min(length(matrixProfile), neighbor + exclusionLen - 1);
            corrProfile(exclRangeBegin : exclRangeEnd) = NaN;
        end
    end
    % Matlab allows NaNs to be ignored from min/max calculations by default, so this
    % is a Matlab specific way to iteratively add excluded ranges of elements. 
    matrixProfile(isnan(corrProfile)) = NaN;
end
end


function [ res ] = moving_mean(a,w)
% moving mean over sequence a with window length w
% based on Ogita et. al, Accurate Sum and Dot Product

% A major source of rounding error is accumulated error in the mean values, so we use this to compensate. 
% While the error bound is still a function of the conditioning of a very long dot product, we have observed 
% a reduction of 3 - 4 digits lost to numerical roundoff when compared to older solutions.

res = zeros(length(a) - w + 1, 1);
p = a(1);
s = 0;

for i = 2 : w
    x = p + a(i);
    z = x - p;
    s = s + ((p - (x - z)) + (a(i) - z));
    p = x;
end

res(1) = p + s;

for i = w + 1 : length(a)
    x = p - a(i - w);
    z = x - p;
    s = s + ((p - (x - z)) - (a(i - w) + z));
    p = x;
    
    x = p + a(i);
    z = x - p;
    s = s + ((p - (x - z)) + (a(i) - z));
    p = x;
    
    res(i - w + 1) = p + s;
end

res = res ./ w;

end
