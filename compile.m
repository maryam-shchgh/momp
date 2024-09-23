clc;
% Options: 
%           Synthetic (syn)
%           Corchorus Capsularis (corcap)
%           data_341_col_2 (data_2)
%           Pure Random walk (pure_rwalk)


dataset_name = 'test';
scores = [];

switch dataset_name
    
    case 'genome'
        mf = matfile('3-data/chimp17.mat');
        m = 1024;
        [momp_out, momp_loc, scores] = momp_v7(mf, m , 1, 0, 0);

    case 'ktip'
        rng(934894); % set random seed
        x = (1:0.05:80);
        motif = sin(x)*5; m = 512;
        length(motif)
        rwalk = smooth((cumsum(randn(1,2^14))),32);  % make data 
        [T, ~] = cat_timeseris(rwalk(1:4000), motif');
        [T, ~] = cat_timeseris(T, rwalk(4001:12000));
        [T, ~] = cat_timeseris(T, motif');
        [T, ~] = cat_timeseris(T, rwalk(12001:end));
        T = T + randn(length(T),1)/10;
        
        ktip_v3(T, m, m/4);
    
    case 'test'
        
        % rng(20242024); % set random seed
        % T1 = smooth((cumsum(randn(1,2^16))),128);  % make data 
        % rng(20242025); % set random seed
        % T2 = smooth((cumsum(randn(1,2^16))),128);  % make data 
        m = 2048;
        %FIX DSR in MOMP
        [momp_out, momp_loc] = mist_v02(data, m, 1, 0);
        

    case 'syn'
        rng(122);
        x = (1:0.01:20);
        motif = sin(x)*10; m = 512;
        rwalk = smooth((cumsum(randn(1,2^14))),32);
        [T, ~] = cat_timeseris(rwalk(1:4000), motif');
        [T, ~] = cat_timeseris(T, rwalk(4001:12000));
        [T, ~] = cat_timeseris(T, motif');
        [T, ~] = cat_timeseris(T, rwalk(12001:end));
        T = T + randn(length(T),1)/20;
        [momp_out, momp_loc] = momp_v8(T, m, 1, 0, 0, 8);
        momp_out(momp_out <0 ) = 0;
        
        
    case 'pure_rwalk'
        for ii=1:100
            display(ii)
            T = genData('pure_rwalk');
            [momp_out, momp_loc, scores(:, end+1:end+2)] = momp_v5(T, 4096, 0, 0, 0);
        end

    case 'corcap'

        T = load('3-data/Corchorus_capsularis.mat');
        T = T.Corchorus_capsularis;
        
        m = 16384;
        dd = 2;
        mem = 1;
        
        [momp_out, momp_loc] = momp_v8(T(1:dd:end)', m/dd, 1);
        %[momp_out, momp_loc] = momp_v9_global(T(1:dd:end)', m/dd, 1, inf, nan, 1);

    case 'data_2'
        
        T = load('3-data/3410222_0003m_125_1.mat');
        T = T.val;
        T = T(2,1:end);
        m = 256*8;
        [momp_out, momp_loc] = momp_v9(T, m, 1, 0, 0);
        
    case 'human'
        [T] = findFixSteps(T);
        m = 65536;
        [momp_out, momp_loc, dsr_final] = momp_v9_global(T, m, 0.000, 1, 0, 0);

        
    case 'ab-join-syn'
        rng(934894); % set random seed
        x = (1:0.05:20);
        motif = sin(x)*10; m = 256;
        rwalk = smooth((cumsum(randn(1,2^16))),32);  % make data 
        [T, ~] = cat_timeseris(rwalk(1:4000), motif');
        [T, ~] = cat_timeseris(T, rwalk(4001:12000));
        [T, ~] = cat_timeseris(T, motif');
        [T, ~] = cat_timeseris(T, rwalk(12001:40000));
        [T, ~] = cat_timeseris(T, motif');
        [T, ~] = cat_timeseris(T, rwalk(40001:end));
        T = T + randn(length(T),1)/20;
        split = 20000;
        [momp_out, momp_loc] = momp_v7_AB(T(1:split-1), 256, T(split:end), 1, 0, 0);
        
        
    case 'human-chimp'
        human = load('3-data/TEMP_DNA/dna/21.mat');
        human = human.Timeseries;
        [human] = findFixSteps(human); human21 = double(human); clear human;
        chimp = load('3-data/cchimp21.mat');
        chimp = chimp.Timeseries;
        [chimp] = findFixSteps(chimp); chimp21 = double(chimp); clear chimp;
        A = human21;
        B = chimp21;
        [T, split] = cat_timeseris(A(1:32:end), B(1:32:end));
        [momp_out, momp_loc] = momp_v9_AB(T, (4096*16), split, 1, 0, 0);
        
end