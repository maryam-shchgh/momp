
function T = genData(n, type)

    if nargin < 2 || isempty(type)
        type = 'rwalk'; 
    end
    switch type
        case 'rwalk'
            T = generate_random_walk(n);
        case 'noise'
            T = randn(n,1);
        otherwise
            error('Unknown data type specified');
    end


end


function T = generate_random_walk(n)
    rng(0);
    mlen = 512;
    rwalk = zeros(n, 1);
    step_size = 0.5;  % Limit step size
    for i = 2:n
        step = randn * step_size;  % Generate a random step within the specified size
        rwalk(i) = rwalk(i-1) + step;
    end
    rwalk = smooth(rwalk, 20);
    x = (0:5:mlen);
    motif = (5*sin(x)-1);
    size(motif)
    figure; plot(motif);
    ii = floor(n/4); jj = 3*floor(n/4);
    rwalk(ii:ii+length(motif)-1) = motif;
    rwalk(jj:jj+length(motif)-1) = motif;
    T = zscore(rwalk + randn(length(rwalk),1)/10);
%     T = smooth(T,10);
end