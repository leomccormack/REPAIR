% -------------------------------------------------------------------------
% This file is part of REPAIR
% Copyright (c) 2021 - Leo McCormack, Archontis Politis & Nils Meyer-Kahlen
%
% REPAIR is free software; you can redistribute it and/or modify it under 
% the terms of the GNU General Public License as published by the Free 
% Software Foundation; either version 2 of the License, or (at your option) 
% any later version.
%
% REPAIR is distributed in the hope that it will be useful, but WITHOUT 
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
% FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
% more details.
%
% See <http://www.gnu.org/licenses/> for a copy of the GNU General Public
% License.
% -------------------------------------------------------------------------
%
% Unit tests for Reproduction and Parameterisation of Array Impulse 
% Responses (REPAIR) [1]
% 
% DEPENDENCES
%   Spherical-Harmonic-Transform Matlab library
%       https://github.com/polarch/Spherical-Harmonic-Transform
%   Higher-Order-Ambisonics Matlab library
%       https://github.com/polarch/Higher-Order-Ambisonics
%   Vector-Base-Amplitude-Panning
%       https://github.com/polarch/Vector-Base-Amplitude-Panning
%   Spherical-Array-Processing
%       https://github.com/polarch/Spherical-Array-Processing
%   shoebox-roomsim
%       https://github.com/polarch/shoebox-roomsim
%   HO-SIRR
%       https://github.com/leomccormack/HO-SIRR
%   (Optional) SDM Toolbox
%       https://se.mathworks.com/matlabcentral/fileexchange/56663-sdm-toolbox
%
% REFERENCES
%   [1] Submitted for review
%
% -------------------------------------------------------------------------
%
%   Leo McCormack, 01/11/2021
%   leo.mccormack@aalto.fi 
%
% -------------------------------------------------------------------------
clear all, close all, dbstop if error %#ok
if ~exist('render_array_rirs', 'file'), error('shoebox-roomsim not in path'); end
addpath('utils')

fs = 48e3;
ENABLE_PLOTTING = 1;
path_for_plots = [];                             % set to '[]' to not save plots, or, e.g. './output/plots'
path_for_renders = [];%'./output/renders';       % set to '[]' to not save renders,
path_for_renders_lt = [];%'./output/renders_lt'; % set to '[]' to not save listening test renders,
ideal_SH_order = 4;
mic_arrays = {'eigenmike32','tetra','intensity-probe'};  % Options: {'eigenmike32', 'tetra', 'intensity-probe'}
t_designs = load('n_designs_1_124'); % Courtesy of Chris Hold (https://github.com/chris-hld/spaudiopy)
signalLength = fs/8;
sofa_file = 'D1_48K_24bit_256tap_FIR_SOFA_KU100.sofa'; % Can be obtained from the SADIE HRIR database, or any other SOFA file can be used

% Default REPAIR configuration 
pars.grid_svecs = []; % Defined for each test (can be e.g. SH weights, space-domain steering vectors, etc.)
pars.grid_dirs_xyz = t_designs.N060; 
pars.grid_dirs_rad = unitCart2sph(pars.grid_dirs_xyz);
pars.grid_weights = findGridWeights(pars.grid_dirs_rad(:,1), pi/2-pars.grid_dirs_rad(:,2))./(4*pi); %(1/size(pars.grid_dirs_rad,1)).*ones(size(pars.grid_dirs_rad,1),1);
pars.fs = fs;  
pars.SCMavgOption = 'recur';  % Options: {'block', 'recur', 'alltime'}
pars.SCMavg_coeff = 0.5; % Temporal averaging coefficient, [0..1], if SCMavgOption is set to "recur"
pars.SCMavg_Nframes = 1;  % Number of frames in each averaging block, if SCMavgOption is set to "block"
pars.Kestimator = 'RECON';  % Options: {'SORTE', 'SORTED', 'RECON', 'ORACLE'}
pars.DoAestimator = 'MUSIC'; % Options: {'MUSIC', 'SRP', 'ORACLE'}
pars.winsize = 256;     % Window size, in time-domain samples
pars.freqGrouping = 'octave'; % Options: {'broadband', 'octave', 'erb', 'fullres'}
pars.streamBalance = 1; % 0: only diffuse stream, 1: both streams are balanced, 2: only direct stream
pars.ENABLE_DIFF_WHITENING = 1;              % Applies an operation that diagonalises the SCMs when under diffuse conditions 
pars.ENABLE_COHERENT_FOCUSING = 1;           % Only used if the steering vectors are frequency dependent, and if there is some band grouping
pars.ENABLE_AMBIENT_ENERGY_PRESERVATION = 1; % Forces the beamformers used for the ambient stream to be energy-preserving over the sphere
pars.decorrelation = 'covMatch'; % Options: {'off', 'convNoise', 'shapedNoise', 'phaseRand', 'covMatch'}
pars.beamformerOption = 'SD';    % Options: {'pinv', 'MF', 'SD'}
pars.ENABLE_QUANTISE_TO_NEAREST_LS = 0;      % Quantise to nearest loudspeaker instead of using VBAP
pars.maxAnaFreq_Hz = fs/2;       % above this frequency, everything is treated as one band
pars.ls_dirs_deg = 180/pi.*unitCart2sph(t_designs.N008);
pars.vbapNorm = 1; % 0:reverberant room, ~0.5: dry listening room, 1: anechoic

% Create output folder for the renders
if ~isempty(path_for_renders), if ~exist(path_for_renders, 'dir'), mkdir(path_for_renders); end, end
if ~isempty(path_for_renders_lt), if ~exist(path_for_renders_lt, 'dir'), mkdir(path_for_renders_lt); end, end
path_for_renders_lt_bin = [path_for_renders_lt filesep 'bin'];
if ~isempty(path_for_renders_lt), if ~exist(path_for_renders_lt_bin, 'dir'), mkdir(path_for_renders_lt_bin); end, end

% Array specifications and steering vectors for the enabled microphone arrays
for mi = 1:length(mic_arrays)
    mic_array = mic_arrays{mi};
    
    % Simulate array responses for the scanning/beamforming grid
    switch mic_array
        case 'eigenmike32'
            mic_spec{mi}.name = 'eigenmike32'; %#ok
            mic_spec{mi}.radius = 0.042; mic_spec{mi}.type = 'rigid';  mic_spec{mi}.dir_coeff = []; mic_spec{mi}.sh_order = 4; %#ok
            mic_spec{mi}.dirs_rad = pi/180.*[0 21; 32 0; 0 -21; 328 0; 0 58; 45 35; 69 0; 45 -35; 0 -58;  315 -35; 291 0; 315 35; 91   69; 90   32; 90  -31; 89  -69; 180  21; 212  0; 180 -21; 148  0; 180 58; 225  35; 249   0; 225 -35; 180 -58; 135 -35; 111  0; 135  35; 269 69; 270 32; 270 -32; 271 -69;]; %#ok
        case 'tetra'
            mic_spec{mi}.name = 'tetra'; %#ok
            mic_spec{mi}.radius = 0.02; mic_spec{mi}.type = 'directional'; mic_spec{mi}.dir_coeff = 0.5; mic_spec{mi}.sh_order = 1; %#ok
            mic_spec{mi}.dirs_rad = [0.785398163397448, 0.615472907423280; -0.785398163397448, -0.615472907423280; 2.35619449019235, -0.615472907423280; -2.35619449019235, 0.615472907423280;]; %#ok
        case 'intensity-probe'
            mic_spec{mi}.name = 'intensity-probe'; %#ok
            mic_spec{mi}.radius = 0.025; mic_spec{mi}.type = 'open'; mic_spec{mi}.dir_coeff = 1; mic_spec{mi}.sh_order = 1; %#ok
            mic_spec{mi}.dirs_rad = unitCart2sph([0.025 0 0; -0.025 0 0; 0 0.025 0; 0 -0.025 0; 0 0 0.025; 0 0 -0.025]); %#ok
    end
    [mic_spec{mi}.h_grid mic_spec{mi}.H_grid] = simulateSphArray(pars.winsize, mic_spec{mi}.dirs_rad, pars.grid_dirs_rad, mic_spec{mi}.type, mic_spec{mi}.radius, 50, fs, mic_spec{mi}.dir_coeff); %#ok
    
    % Encoding filters/matrix
    [mic_spec{mi}.E_sh, mic_spec{mi}.e_sh] = arraySHTfiltersMeas_regLSHD(mic_spec{mi}.H_grid, mic_spec{mi}.sh_order, pars.grid_dirs_rad, [], pars.winsize, 15); %#ok
    
    % Array SH domain steering vectors 
    for nb=1:pars.winsize/2+1
        mic_spec{mi}.H_grid_sh(nb,:,:) = mic_spec{mi}.E_sh(:,:,nb) * squeeze(mic_spec{mi}.H_grid(nb,:,:));
    end  
    if ENABLE_PLOTTING==1
        figure,
        for ii=1:size(mic_spec{mi}.e_sh,1), for jj=1:size(mic_spec{mi}.e_sh,2), plot(squeeze(mic_spec{mi}.e_sh(ii,jj,:))), hold on, end, end
        title([ mic_spec{mi}.name  ' encoding filters']), grid on
        evaluateSHTfilters(mic_spec{mi}.E_sh, mic_spec{mi}.H_grid, fs, sqrt(4*pi) * getSH(mic_spec{mi}.sh_order, [pars.grid_dirs_rad(:,1) pi/2-pars.grid_dirs_rad(:,2)], 'real'),[]);
    end 
    
    % Permute into the format that REPAIR requires
    mic_spec{mi}.H_grid = permute(mic_spec{mi}.H_grid, [2 3 1]);  %#ok
    mic_spec{mi}.H_grid_sh = permute(mic_spec{mi}.H_grid_sh, [2 3 1]);  %#ok
end

% For binauralising the test rendered outputs
if ~isempty(path_for_renders)
    sofa = loadSofaFile(sofa_file);
    hrirs = sofa.IR;
    hrirs = 0.25.*hrirs./max(abs(hrirs(:)));
    hrir_dirs_rad = pi/180.*sofa.SourcePosition([1 2], :).'; 
    assert(sofa.IR_fs==fs)
    h_ls2bin = permute(hrirs(:,:,findClosestGridPoints(hrir_dirs_rad, pars.ls_dirs_deg*pi/180)), [1 3 2]);
    clear sofa
end

%% Plane-waves (pw) in a free-field
pw_tests = {'one pw', 'two inc pw', 'three inc pw'}; % Options: {'one pw', 'two pw', 'two coh pw', 'three coh pw'}
NFFT = 256; 

for ti = 1:length(pw_tests)
    pw_test = pw_tests{ti}; 
    
    % pw signals and directions:
    switch pw_test
        case 'one pw'
            src_dirs_rad = pi/180.*[58.8941020926610,53.6513674205284];  
            src_sigs = (1/3).*randn(signalLength,size(src_dirs_rad,1));
        case 'two inc pw'
            src_dirs_rad = pi/180.*[58.8941020926610,53.6513674205284; -130.226160991636,-14.0948465161701];  
            src_sigs = (1/3).*randn(signalLength,size(src_dirs_rad,1));
        case 'two coh pw'
            src_dirs_rad = pi/180.*[58.8941020926610,53.6513674205284; -130.226160991636,-14.0948465161701];  
            src_sigs = (1/3).*repmat(randn(signalLength,1), [1 size(src_dirs_rad,1)]);
        case 'three inc pw'
            src_dirs_rad = pi/180.*[58.8941020926610,53.6513674205284; -130.226160991636,-14.0948465161701; 0,-14.0948465161701];
            src_sigs = (1/3).*randn(signalLength,size(src_dirs_rad,1));
        case 'three coh pw'
            src_dirs_rad = pi/180.*[58.8941020926610,53.6513674205284; -130.226160991636,-14.0948465161701; 0,-14.0948465161701];
            src_sigs = (1/3).*repmat(randn(signalLength,1), [1 size(src_dirs_rad,1)]);
    end  
    if ~isempty(path_for_renders) 
        h_src2bin = permute(hrirs(:,:,findClosestGridPoints(hrir_dirs_rad, src_dirs_rad)), [1 3 2]);
        y_bin = matrixConvolver(src_sigs, h_src2bin, size(hrirs,1));
        audiowrite( [path_for_renders filesep 'TEST ' pw_test ' reference' '.wav'], y_bin, fs);
        if ENABLE_PLOTTING==1  
            y_bin_fft = []; for ch = 1:2, y_bin_fft(:,:,ch) = fft(buffer(y_bin(:,ch),NFFT)); end  %#ok
            figure, semilogx(0:fs/NFFT:fs/2, 20*log10(abs(squeeze(mean(abs(y_bin_fft(1:end/2+1,:,:)), 2)))))
            title(['TEST ' pw_test ' reference']), grid on, xlim([100 20e3])
        end 
    end
    
    % Ideal SH receiver
    test_name = ['TEST ' pw_test ' (ideal SH receiver, order ' num2str(ideal_SH_order) ')']; disp(test_name);
    pars.grid_svecs = getSH(ideal_SH_order, [pars.grid_dirs_rad(:,1) pi/2-pars.grid_dirs_rad(:,2)], 'real').'; 
    srir_sh = src_sigs * (getSH(ideal_SH_order, [src_dirs_rad(:,1) pi/2-src_dirs_rad(:,2)], 'real'));
    [h_ls, h_ls_dir, h_ls_diff, analysis] = REPAIR(srir_sh, pars);
    if ENABLE_PLOTTING, PLOT_REPAIR(test_name, h_ls, h_ls_dir, h_ls_diff, pars, analysis, 0, path_for_plots); end
    if ~isempty(path_for_renders)  
        y_bin = matrixConvolver(h_ls, h_ls2bin, size(hrirs,1)); 
        audiowrite( [path_for_renders filesep test_name '.wav'], y_bin, fs);
        if ENABLE_PLOTTING==1
            y_bin_fft = []; for ch = 1:2, y_bin_fft(:,:,ch) = fft(buffer(y_bin(:,ch),NFFT)); end  %#ok
            figure, semilogx(0:fs/NFFT:fs/2, 20*log10(abs(squeeze(mean(abs(y_bin_fft(1:end/2+1,:,:)), 2)))))
            title(test_name), grid on, xlim([100 20e3])
        end
    end

    % Now loop over the microphone arrays
    for mi = 1:length(mic_arrays)
        % Space domain:
        test_name = ['TEST ' pw_test ' (array ' mic_spec{mi}.name ', array steering vectors [space domain])']; disp(test_name); 
        pars.grid_svecs = mic_spec{mi}.H_grid;
        h_src = simulateSphArray(pars.winsize, mic_spec{mi}.dirs_rad, src_dirs_rad, mic_spec{mi}.type, mic_spec{mi}.radius, 25, fs, mic_spec{mi}.dir_coeff);  
        srir = matrixConvolver(src_sigs, permute(h_src, [1 3 2]), pars.winsize);
        [h_ls, h_ls_dir, h_ls_diff, analysis] = REPAIR(srir, pars);
        if ENABLE_PLOTTING, PLOT_REPAIR(test_name, h_ls, h_ls_dir, h_ls_diff, pars, analysis, 0, path_for_plots); end
        if ~isempty(path_for_renders)  
            y_bin = matrixConvolver(h_ls, h_ls2bin, size(hrirs,1)); 
            audiowrite( [path_for_renders filesep test_name '.wav'], y_bin, fs);
            if ENABLE_PLOTTING==1
                y_bin_fft = []; for ch = 1:2, y_bin_fft(:,:,ch) = fft(buffer(y_bin(:,ch),NFFT)); end  %#ok
                figure, semilogx(0:fs/NFFT:fs/2, 20*log10(abs(squeeze(mean(abs(y_bin_fft(1:end/2+1,:,:)), 2)))))
                title(test_name), grid on, xlim([100 20e3])
            end
        end

        if ~strcmp(mic_spec{mi}.type, 'open')
            % SH domain (using broad-band SHs as steering vectors):
            test_name = ['TEST ' pw_test ' (array ' mic_spec{mi}.name ', broad-band SH steering vectors)']; disp(test_name); 
            pars.grid_svecs = getSH(mic_spec{mi}.sh_order, [pars.grid_dirs_rad(:,1) pi/2-pars.grid_dirs_rad(:,2)], 'real').'; 
            srir_sh = matrixConvolver(srir, permute(mic_spec{mi}.e_sh, [3 2 1]), pars.winsize); 
            [h_ls, h_ls_dir, h_ls_diff, analysis] = REPAIR(srir_sh, pars);
            if ENABLE_PLOTTING, PLOT_REPAIR(test_name, h_ls, h_ls_dir, h_ls_diff, pars, analysis, 0, path_for_plots); end
            if ~isempty(path_for_renders)  
                audiowrite( [path_for_renders filesep test_name '.wav'], matrixConvolver(h_ls, h_ls2bin, size(hrirs,1)), fs);
                if ENABLE_PLOTTING==1
                    y_bin_fft = []; for ch = 1:2, y_bin_fft(:,:,ch) = fft(buffer(y_bin(:,ch),NFFT)); end  %#ok
                    figure, semilogx(0:fs/NFFT:fs/2, 20*log10(abs(squeeze(mean(abs(y_bin_fft(1:end/2+1,:,:)), 2)))))
                    title(test_name), grid on, xlim([100 20e3])
                end
            end

            % SH domain (using the encoded SMA steering vectors):
            test_name = ['TEST ' pw_test ' (array ' mic_spec{mi}.name ', SH encoded array steering vectors)']; disp(test_name); 
            pars.grid_svecs = mic_spec{mi}.H_grid_sh; 
            [h_ls, h_ls_dir, h_ls_diff, analysis] = REPAIR(srir_sh, pars);
            if ENABLE_PLOTTING, PLOT_REPAIR(test_name, h_ls, h_ls_dir, h_ls_diff, pars, analysis, 0, path_for_plots); end
            if ~isempty(path_for_renders)  
                audiowrite( [path_for_renders filesep test_name '.wav'], matrixConvolver(h_ls, h_ls2bin, size(hrirs,1)), fs);
                if ENABLE_PLOTTING==1
                    y_bin_fft = []; for ch = 1:2, y_bin_fft(:,:,ch) = fft(buffer(y_bin(:,ch),NFFT)); end  %#ok
                    figure, semilogx(0:fs/NFFT:fs/2, 20*log10(abs(squeeze(mean(abs(y_bin_fft(1:end/2+1,:,:)), 2)))))
                    title(test_name), grid on, xlim([100 20e3])
                end
            end
        end
    end
end


%% Diffuse-field input
NFFT = 4096;
diff_dirs_xyz = t_designs.N040; 
diff_dirs_rad = unitCart2sph(diff_dirs_xyz);

% Reference
diff_sigs = (1/3).*randn(signalLength,size(diff_dirs_rad,1));
diff_sigs = sqrt(1/size(diff_dirs_rad,1)) .* diff_sigs;
if ~isempty(path_for_renders) 
    h_diff2bin = permute(hrirs(:,:,findClosestGridPoints(hrir_dirs_rad, diff_dirs_rad)), [1 3 2]);
    y_bin=matrixConvolver(diff_sigs, h_diff2bin, size(hrirs,1)); 
    if ENABLE_PLOTTING==1 
        y_bin_fft = []; for ch = 1:2, y_bin_fft(:,:,ch) = fft(buffer(y_bin(:,ch),NFFT)); end  %#ok
        figure, semilogx(0:fs/NFFT:fs/2, 20*log10(abs(squeeze(mean(abs(y_bin_fft(1:end/2+1,:,:)), 2)))))
        title('TEST diffusefield reference'), grid on, xlim([100 20e3])
    end  
    audiowrite( [path_for_renders filesep 'TEST diffusefield reference' '.wav'], y_bin, fs);
end

% Ideal SH receiver
test_name = ['TEST diffusefield (ideal SH receiver, order ' num2str(ideal_SH_order) ')']; disp(test_name);
pars.grid_svecs = getSH(ideal_SH_order, [pars.grid_dirs_rad(:,1) pi/2-pars.grid_dirs_rad(:,2)], 'real').'; 
srir_sh = diff_sigs * (getSH(ideal_SH_order, [diff_dirs_rad(:,1) pi/2-diff_dirs_rad(:,2)], 'real'));
[h_ls, h_ls_dir, h_ls_diff, analysis] = REPAIR(srir_sh, pars);
if ENABLE_PLOTTING==1
    PLOT_REPAIR(test_name, h_ls, h_ls_dir, h_ls_diff, pars, analysis, 0, path_for_plots); 
    y_bin=matrixConvolver(h_ls, h_ls2bin, size(hrirs,1)); 
    y_bin_fft = []; for ch = 1:2, y_bin_fft(:,:,ch) = fft(buffer(y_bin(:,ch),NFFT)); end  %#ok
    figure, semilogx(0:fs/NFFT:fs/2, 20*log10(abs(squeeze(mean(abs(y_bin_fft(1:end/2+1,:,:)), 2)))))
    title(test_name), grid on, xlim([100 20e3])
    figure, imagesc(h_ls_diff'*h_ls_diff), title(pars.decorrelation)
end
if ~isempty(path_for_renders)  
    audiowrite( [path_for_renders filesep test_name '.wav'], matrixConvolver(h_ls, h_ls2bin, size(hrirs,1)), fs);
end

% Now loop over the microphone arrays
for mi = 1:length(mic_arrays)
    % Space domain:
    test_name = ['TEST diffusefield (array ' mic_spec{mi}.name ', array steering vectors [space domain])']; disp(test_name); 
    pars.grid_svecs = mic_spec{mi}.H_grid;
    h_diff = simulateSphArray(pars.winsize, mic_spec{mi}.dirs_rad, diff_dirs_rad, mic_spec{mi}.type, mic_spec{mi}.radius, 25, fs, mic_spec{mi}.dir_coeff);  
    srir = matrixConvolver(diff_sigs, permute(h_diff, [1 3 2]), pars.winsize);
    [h_ls, h_ls_dir, h_ls_diff, analysis] = REPAIR(srir, pars);
    if ENABLE_PLOTTING==1
        PLOT_REPAIR(test_name, h_ls, h_ls_dir, h_ls_diff, pars, analysis, 0, path_for_plots); 
        y_bin=matrixConvolver(h_ls, h_ls2bin, size(hrirs,1)); 
        y_bin_fft = []; for ch = 1:2, y_bin_fft(:,:,ch) = fft(buffer(y_bin(:,ch),NFFT)); end  %#ok
        figure, semilogx(0:fs/NFFT:fs/2, 20*log10(abs(squeeze(mean(abs(y_bin_fft(1:end/2+1,:,:)), 2)))))
        title(test_name), grid on, xlim([100 20e3])
    end
    if ~isempty(path_for_renders)
        audiowrite( [path_for_renders filesep test_name '.wav'], matrixConvolver(h_ls, h_ls2bin, size(hrirs,1)), fs);
    end
   
    if ~strcmp(mic_spec{mi}.type, 'open')
        % SH domain (using broad-band SHs as steering vectors):
        test_name = ['TEST diffusefield (array ' mic_spec{mi}.name ', broad-band SH steering vectors)']; disp(test_name); 
        pars.grid_svecs = getSH(mic_spec{mi}.sh_order, [pars.grid_dirs_rad(:,1) pi/2-pars.grid_dirs_rad(:,2)], 'real').'; 
        srir_sh = matrixConvolver(srir, permute(mic_spec{mi}.e_sh, [3 2 1]), pars.winsize); 
        [h_ls, h_ls_dir, h_ls_diff, analysis] = REPAIR(srir_sh, pars);
        if ENABLE_PLOTTING==1
            PLOT_REPAIR(test_name, h_ls, h_ls_dir, h_ls_diff, pars, analysis, 0, path_for_plots); 
            y_bin=matrixConvolver(h_ls, h_ls2bin, size(hrirs,1)); 
            y_bin_fft = []; for ch = 1:2, y_bin_fft(:,:,ch) = fft(buffer(y_bin(:,ch),NFFT)); end  %#ok
            figure, semilogx(0:fs/NFFT:fs/2, 20*log10(abs(squeeze(mean(abs(y_bin_fft(1:end/2+1,:,:)), 2)))))
            title(test_name), grid on, xlim([100 20e3])
        end
        if ~isempty(path_for_renders)  
            audiowrite( [path_for_renders filesep test_name '.wav'], matrixConvolver(h_ls, h_ls2bin, size(hrirs,1)), fs);
        end

        % SH domain (using the encoded SMA steering vectors):
        test_name = ['TEST diffusefield (array ' mic_spec{mi}.name ', SH encoded array steering vectors)']; disp(test_name); 
        pars.grid_svecs = mic_spec{mi}.H_grid_sh; 
        [h_ls, h_ls_dir, h_ls_diff, analysis] = REPAIR(srir_sh, pars);
        if ENABLE_PLOTTING==1
            PLOT_REPAIR(test_name, h_ls, h_ls_dir, h_ls_diff, pars, analysis, 0, path_for_plots); 
            y_bin=matrixConvolver(h_ls, h_ls2bin, size(hrirs,1)); 
            y_bin_fft = []; for ch = 1:2, y_bin_fft(:,:,ch) = fft(buffer(y_bin(:,ch),NFFT)); end  %#ok
            figure, semilogx(0:fs/NFFT:fs/2, 20*log10(abs(squeeze(mean(abs(y_bin_fft(1:end/2+1,:,:)), 2)))))
            title(test_name), grid on, xlim([100 20e3])
        end
        if ~isempty(path_for_renders)  
            audiowrite( [path_for_renders filesep test_name '.wav'], matrixConvolver(h_ls, h_ls2bin, size(hrirs,1)), fs);
        end
    end
end


%% Image-source based RIRs
rooms = {'small', 'medium','large'}; % Options: {'small', 'medium', 'large'}

% Optionally, output listening test files
if ~isempty(path_for_renders_lt) 
    % Options:
    ENABLE_HOSIRR_RENDERS = 1;
    ENABLE_SDM_RENDERS = 1; 
    
    % Load anechoic stimuli
    addpath('../HO-SIRR') 
    addpath('../HO-SIRR/_Stimuli_') 
    if ~exist('HOSIRR', 'file'), error('HOSIRR not in path'); end 
    addpath('_Stimuli_')
    stimuli = {'music__KickDrumClicky'
        'music__Trombone'
        'music__Castanets'
        'speech__HarvardMale'
        'music__snares'};
    for si=1:length(stimuli)
        [srcsig_tmp, srcsig_fs] = audioread([stimuli{si} '.wav']); assert(srcsig_fs==fs)
        srcsig_tmp = [srcsig_tmp; zeros(fs/2,1)];
        srcsig{si} = 0.7.*srcsig_tmp./max(abs(srcsig_tmp(:))); %#ok
        
        hpf_cutoff = 60;
        w_hh = hpf_cutoff/(fs/2);
        h_filt = fir1(10000, w_hh, 'high').';
        srcsig{si} = fftfilt(repmat(h_filt, [1 size(srcsig{si},2)]), [srcsig{si}; zeros(5001,size(srcsig{si},2))]); %#ok
        srcsig{si} = srcsig{si}(5000:end,:); %#ok
    end
    
    % For binaural rendering
    h_ls2bin = 3.*permute(hrirs(:,:,findClosestGridPoints(hrir_dirs_rad, pars.ls_dirs_deg*pi/180)), [1 3 2]);
    
    % Configurations for other methods:
    if ENABLE_HOSIRR_RENDERS 
        addpath('../HO-SIRR') 
        if ~exist('HOSIRR', 'file'), error('HOSIRR not in path'); end 
        % Default HOSIRR toolbox configuration
        hosirr_pars.chOrdering = 'ACN'; 
        hosirr_pars.normScheme = 'N3D'; 
        hosirr_pars.fs = fs;  
        hosirr_pars.multires_winsize = 128;  
        hosirr_pars.multires_xovers = [ ];   
        hosirr_pars.RENDER_DIFFUSE = 1;
        hosirr_pars.decorrelationType = 'noise';
        hosirr_pars.BROADBAND_FIRST_PEAK = 1;  
        hosirr_pars.BROADBAND_DIFFUSENESS = 1;
        hosirr_pars.maxDiffFreq_Hz = 3000;  
        hosirr_pars.alpha_diff = 0.5; 
        hosirr_pars.ls_dirs_deg = pars.ls_dirs_deg;
    end
    if ENABLE_SDM_RENDERS 
        if ~exist('createSDMStruct', 'file'), error('SDM not in path'); end
        % Default SDM toolbox space-domain configuration
        sdm_pars = createSDMStruct('micLocs', [0.025 0 0; -0.025 0 0; 0 0.025 0; 0 -0.025 0; 0 0 0.025; 0 0 -0.025], 'fs',fs);
        ls_dirs_deg_r = [pars.ls_dirs_deg ones(size(pars.ls_dirs_deg,1),1)]; % add radius
    end
end

ls_dirs_rad = pars.ls_dirs_deg*pi/180;   

for ri = 1:length(rooms)
    room = rooms{ri};

    % Configure room
    switch room
        case 'small'
            room_dims = [6 5 3.1];
            rt60 = [0.25 0.3 0.2 0.15 0.05 0.03].*1.3;
        case 'medium'
            room_dims = [10 6 3.5];
            rt60 = [0.4 0.45 0.3 0.15 0.12 0.1].*1.3;
        case 'large'
            room_dims = [12 7.4 3.5];
            rt60 = [0.4 0.47 0.35 0.18 0.13 0.12].*1.4;
    end  
    abs_wall_ratios = [0.75 0.86 0.56 0.95 0.88 1];
    nBands = length(rt60);
    band_centerfreqs = zeros(nBands,1);
    band_centerfreqs(1) = 125;
    for nb=2:nBands, band_centerfreqs(nb) = 2*band_centerfreqs(nb-1); end 
    abs_wall = findAbsCoeffsFromRT(room_dims, rt60,abs_wall_ratios);
    [RT60_sabine,  d_critical] = room_stats(room_dims, abs_wall);
    rec_pos_xyz = room_dims/2 + [-0.42 -0.44 0.18].*1.3; 
    src_pos_xyz = rec_pos_xyz + [2 0.24 -0.15].*1.1; 
    
    % Generate reference loudspeaker RIRs  
    type = 'maxTime'; % 'maxTime' 'maxOrder' 
    maxlim = max(rt60).*ones(size(rt60))*2; % just cut if it's longer than that ( or set to max(rt60) )  
    src_pos_xyz2 = [src_pos_xyz(:,1) room_dims(2)-src_pos_xyz(:,2) src_pos_xyz(:,3)]; % change y coord for src/rec due to convention inside the IMS function
    rec_pos_xyz2 = [rec_pos_xyz(:,1) room_dims(2)-rec_pos_xyz(:,2) rec_pos_xyz(:,3)]; % change y coord for src/rec due to convention inside the IMS function
    abs_echograms = compute_echograms_arrays(room_dims, src_pos_xyz2, rec_pos_xyz2, abs_wall, maxlim);
    grid_cell{1} = ls_dirs_rad;
    array_irs{1} = [permute(eye(size(ls_dirs_rad,1)), [3 1 2]); zeros(7,size(ls_dirs_rad,1),size(ls_dirs_rad,1));]; % zero pad
    x_rir = render_array_rirs(abs_echograms, band_centerfreqs, fs, grid_cell, array_irs);  x_rir = x_rir{1}; % remove rec cell 
    x_rir = x_rir((500+8):end,:); % remove FIR filter incurred delays
    for ii=1:size(x_rir,2), x_rir_air(:,ii) = applyAirAbsorption(x_rir(:,ii), fs); end %#ok
    x_rir = x_rir_air;clear x_rir_air
    if ENABLE_PLOTTING==1
        figure, mesh(max(20*log10(abs(x_rir)), -49.8))
        colormap('jet'), xlim([1, size(pars.ls_dirs_deg,1)]), zlim([-50,10 ])
        title(['reference ' room ' room loudspeaker RIR']), xlabel('loudspeaker #'), zlabel('energy, dB'), ylabel('time, samples')
    end
    if ~isempty(path_for_renders)  
        audiowrite( [path_for_renders filesep 'TEST ' room ' room reference' '.wav'], matrixConvolver(x_rir, h_ls2bin, size(hrirs,1)), fs, 'BitsPerSample', 32);
    end 
    if ~isempty(path_for_renders_lt) 
        for si=1:length(stimuli)
            yi = fftfilt(x_rir, repmat(srcsig{si}, [1 size(pars.ls_dirs_deg,1)]));      
            %audiowrite( [path_for_renders_lt filesep room ' ' stimuli{si} ' ref' '.wav'], yi, fs, 'BitsPerSample', 32);
            audiowrite( [path_for_renders_lt_bin filesep room ' ' stimuli{si} ' ref' '.wav'], matrixConvolver(yi, h_ls2bin, size(hrirs,1)), fs, 'BitsPerSample', 32);
        end
    end
    
    % Obtain Oracle parameters
    maxK_oracle = 8;
    time_s = abs_echograms(1,1,1).time; 
    doas_s = abs_echograms(1,1,1).coords;
    values = abs_echograms(1,1,1).value;
    time_hop_index = floor((time_s*fs)/(pars.winsize/2));
    for ii=1:time_hop_index(end)
        idxs = find(ii==time_hop_index); 
        if ~isempty(idxs)
            [~,sortdec] = sort(values(idxs),'descend');
            time_s_tmp = time_s(idxs);
            doas_s_tmp = doas_s(idxs,:);
            values_tmp = values(idxs); 
            time_s(idxs) = time_s_tmp(sortdec);
            doas_s(idxs,:) = doas_s_tmp(sortdec,:);
            values(idxs) = values_tmp(sortdec); 
        end 
    end 
    nFrames = time_hop_index(end)+40;
    analysis_oracle.K = zeros(1,nFrames);
    analysis_oracle.est_idx = nan(maxK_oracle,nFrames);
    for i=1:length(time_hop_index)
        index = time_hop_index(i);
        src_dir_rad = unitCart2sph(doas_s(i,:)); 
        est_idx = findClosestGridPoints(pars.grid_dirs_rad, src_dir_rad);
        if analysis_oracle.K(1,index)<maxK_oracle
            analysis_oracle.K(1,index) = analysis_oracle.K(1,index) + 1; 
            analysis_oracle.est_idx(analysis_oracle.K(1,index),index) = est_idx;
        end
    end
    %if ENABLE_PLOTTING, figure, plot(time_s*fs, abs_echograms(1,1,1).value), end
    %if ENABLE_PLOTTING, figure, plot(x_rir(:,1)), end
   
    % Ideal SH receiver  
    test_name = ['TEST ' room ' room (ideal SH receiver, order ' num2str(ideal_SH_order) ' ' pars.Kestimator ')']; disp(test_name);
    pars.grid_svecs = getSH(ideal_SH_order, [pars.grid_dirs_rad(:,1) pi/2-pars.grid_dirs_rad(:,2)], 'real').'; 
    x_rir_sh = x_rir * (getSH(ideal_SH_order, [ls_dirs_rad(:,1) pi/2-ls_dirs_rad(:,2)], 'real'));
    [h_ls, h_ls_dir, h_ls_diff, analysis] = REPAIR(x_rir_sh, pars, analysis_oracle); 
    if ENABLE_PLOTTING, PLOT_REPAIR(test_name, h_ls, h_ls_dir, h_ls_diff, pars, analysis, 0, path_for_plots); end
    if ~isempty(path_for_renders_lt) 
        for si=1:length(stimuli)
            yi = fftfilt(h_ls, repmat(srcsig{si},[1 size(pars.ls_dirs_deg,1)]));      
            %audiowrite( [path_for_renders_lt filesep room ' ' stimuli{si} ' idealSH o' num2str(ideal_SH_order) ' ' pars.Kestimator ' REPAIR SHD.wav'], yi, fs, 'BitsPerSample', 32);
            ybin = matrixConvolver(yi, h_ls2bin, size(hrirs,1));
            if max(abs(ybin(:)))>1, ybin = ybin./max(abs(ybin(:))); end
            audiowrite( [path_for_renders_lt_bin filesep  room ' ' stimuli{si} ' idealSH o' num2str(ideal_SH_order) ' ' pars.Kestimator ' REPAIR SHD.wav'], ybin, fs, 'BitsPerSample', 32);
        end
    end 

    % HOSIRR Ideal SH receiver 
    if ENABLE_HOSIRR_RENDERS
        case_name = ['HOSIRR ' room ' room (ideal SH receiver, order ' num2str(ideal_SH_order) ')']; disp(case_name);
        [h_ls, h_ls_dir, h_ls_diff, analysis] = HOSIRR(x_rir_sh, hosirr_pars);
        if ENABLE_PLOTTING==1, PLOT_REPAIR(case_name, h_ls, h_ls_dir, h_ls_diff, pars, [], 1, path_for_plots); end
        for si=1:length(stimuli)
            yi = fftfilt(h_ls, repmat(srcsig{si},[1 size(hosirr_pars.ls_dirs_deg,1)])) * 2;      
            %audiowrite( [path_for_renders_lt filesep room ' ' stimuli{si} ' idealSH o' num2str(ideal_SH_order)  ' HOSIRR SHD.wav'], yi, fs, 'BitsPerSample', 32);
            audiowrite( [path_for_renders_lt_bin filesep room ' ' stimuli{si} ' idealSH o' num2str(ideal_SH_order)  ' HOSIRR SHD.wav'], matrixConvolver(yi, h_ls2bin, size(hrirs,1)), fs, 'BitsPerSample', 32);
        end
        
        case_name = ['HOSIRR ' room ' room (ideal SH receiver, order ' num2str(1) ')']; disp(case_name);
        [h_ls, h_ls_dir, h_ls_diff, analysis] = HOSIRR(x_rir_sh(:,1:4), hosirr_pars);
        if ENABLE_PLOTTING==1, PLOT_REPAIR(case_name, h_ls, h_ls_dir, h_ls_diff, pars, [], 1, path_for_plots); end
        for si=1:length(stimuli)
            yi = fftfilt(h_ls, repmat(srcsig{si},[1 size(hosirr_pars.ls_dirs_deg,1)])) * 2;      
            %audiowrite( [path_for_renders_lt filesep room ' ' stimuli{si} ' idealSH o' num2str(ideal_SH_order)  ' HOSIRR SHD.wav'], yi, fs, 'BitsPerSample', 32);
            audiowrite( [path_for_renders_lt_bin filesep room ' ' stimuli{si} ' idealSH o' num2str(1)  ' HOSIRR SHD.wav'], matrixConvolver(yi, h_ls2bin, size(hrirs,1)), fs, 'BitsPerSample', 32);
        end
    end
    
    % Now loop over the microphone arrays
    for mi = 1:length(mic_arrays) 
        % Space domain:
        test_name = ['TEST ' room ' room (array ' mic_spec{mi}.name ' ' pars.Kestimator  ', array steering vectors [space domain])']; disp(test_name); 
        pars.grid_svecs = mic_spec{mi}.H_grid;  
        h_ls2mic = simulateSphArray(1024, mic_spec{mi}.dirs_rad, ls_dirs_rad, mic_spec{mi}.type, mic_spec{mi}.radius, 50, fs, mic_spec{mi}.dir_coeff);  
        srir = matrixConvolver(x_rir, permute(h_ls2mic, [1 3 2]), 1024);  
        srir = srir(512:512+size(x_rir,1),:);
        if ENABLE_PLOTTING, figure, plot(srir), end
        [h_ls, h_ls_dir, h_ls_diff, analysis] = REPAIR(srir, pars);
        if ENABLE_PLOTTING==1, PLOT_REPAIR(test_name, h_ls, h_ls_dir, h_ls_diff, pars, analysis, 0, path_for_plots); end
        if ~isempty(path_for_renders_lt) 
            for si=1:length(stimuli)
                yi = fftfilt(h_ls, repmat(srcsig{si},[1 size(pars.ls_dirs_deg,1)]));    
                %audiowrite( [path_for_renders_lt filesep room ' ' stimuli{si} ' ' mic_spec{mi}.name ' ' pars.Kestimator  ' REPAIR SD' '.wav'], yi, fs, 'BitsPerSample', 32);
                ybin = matrixConvolver(yi, h_ls2bin, size(hrirs,1));
                if max(abs(ybin(:)))>1, ybin = ybin./max(abs(ybin(:))); end
                audiowrite( [path_for_renders_lt_bin filesep room ' ' stimuli{si} ' ' mic_spec{mi}.name ' ' pars.Kestimator  ' REPAIR SD' '.wav'], ybin, fs, 'BitsPerSample', 32);
            end
        end 
        
        if ~strcmp(mic_spec{mi}.type, 'open')
            % SH domain (using broad-band SHs as steering vectors):
            test_name = ['TEST ' room ' room (array ' mic_spec{mi}.name ', broad-band SH steering vectors)']; disp(test_name); 
            pars.grid_svecs = getSH(mic_spec{mi}.sh_order, [pars.grid_dirs_rad(:,1) pi/2-pars.grid_dirs_rad(:,2)], 'real').'; 
            srir_sh = matrixConvolver(srir, permute(mic_spec{mi}.e_sh, [3 2 1]), pars.winsize); 
            srir_sh = srir_sh((pars.winsize/2):(pars.winsize/2)+size(x_rir,1),:); 
            if ENABLE_PLOTTING, figure, plot(srir_sh), end
            [h_ls, h_ls_dir, h_ls_diff, analysis] = REPAIR(srir_sh, pars);
            if ENABLE_PLOTTING==1, PLOT_REPAIR(test_name, h_ls, h_ls_dir, h_ls_diff, pars, analysis, 0, path_for_plots); end
            if ~isempty(path_for_renders_lt) 
                for si=1:length(stimuli)
                    yi = fftfilt(h_ls, repmat(srcsig{si},[1 size(pars.ls_dirs_deg,1)]));      
                    %audiowrite( [path_for_renders_lt filesep room ' ' stimuli{si} ' ' mic_spec{mi}.name ' REPAIR SHD BB' '.wav'], yi, fs, 'BitsPerSample', 32);
                    ybin = matrixConvolver(yi, h_ls2bin, size(hrirs,1));
                    if max(abs(ybin(:)))>1, ybin = ybin./max(abs(ybin(:))); end
                    audiowrite( [path_for_renders_lt_bin filesep room ' ' stimuli{si} ' ' mic_spec{mi}.name ' REPAIR SHD BB' '.wav'], ybin, fs, 'BitsPerSample', 32);
                end
            end

            % SH domain (using the encoded SMA steering vectors):
            test_name = ['TEST ' room ' room (array ' mic_spec{mi}.name ', SH encoded array steering vectors)']; disp(test_name); 
            pars.grid_svecs = mic_spec{mi}.H_grid_sh; 
            [h_ls, h_ls_dir, h_ls_diff, analysis] = REPAIR(srir_sh, pars);
            if ENABLE_PLOTTING==1, PLOT_REPAIR(test_name, h_ls, h_ls_dir, h_ls_diff, pars, analysis, 0, path_for_plots); end
            if ~isempty(path_for_renders_lt) 
                for si=1:length(stimuli)
                    yi = fftfilt(h_ls, repmat(srcsig{si},[1 size(pars.ls_dirs_deg,1)]));      
                    %audiowrite( [path_for_renders_lt filesep room ' ' stimuli{si} ' ' mic_spec{mi}.name ' REPAIR SHD E' '.wav'], yi, fs, 'BitsPerSample', 32);
                    ybin = matrixConvolver(yi, h_ls2bin, size(hrirs,1));
                    if max(abs(ybin(:)))>1, ybin = ybin./max(abs(ybin(:))); end
                    audiowrite( [path_for_renders_lt_bin filesep room ' ' stimuli{si} ' ' mic_spec{mi}.name ' REPAIR SHD E' '.wav'], ybin, fs, 'BitsPerSample', 32);
                end
            end
        
            % HOSIRR SH domain 
            if ENABLE_HOSIRR_RENDERS
                case_name = ['HOSIRR ' room ' room (' mic_spec{mi}.name ' broad-band SH steering vectors)']; disp(case_name); 
                [h_ls, h_ls_dir, h_ls_diff, analysis] = HOSIRR(srir_sh, hosirr_pars);
                if ENABLE_PLOTTING==1, PLOT_REPAIR(case_name, h_ls, h_ls_dir, h_ls_diff, pars, [], 1, path_for_plots); end
                for si=1:length(stimuli)
                    yi = fftfilt(h_ls, repmat(srcsig{si},[1 size(hosirr_pars.ls_dirs_deg,1)]));      
                    %audiowrite( [path_for_renders_lt filesep room ' ' stimuli{si} ' ' mic_spec{mi}.name ' HOSIRR SHD' '.wav'], yi, fs, 'BitsPerSample', 32);
                    audiowrite( [path_for_renders_lt_bin filesep room ' ' stimuli{si} ' ' mic_spec{mi}.name ' HOSIRR SHD'  '.wav'], matrixConvolver(yi, h_ls2bin, size(hrirs,1)), fs, 'BitsPerSample', 32);
                end
            end  
        end
        
        if strcmp(mic_spec{mi}.name, 'intensity-probe') && ENABLE_SDM_RENDERS
            case_name = ['SDM ' room ' room (' mic_spec{mi}.name ' Space-domain']; disp(case_name); 
            DOA = SDMPar(srir, sdm_pars);
            P = srir(:, 5);
            s = createSynthesisStruct('lspLocs',ls_dirs_deg_r,'snfft',length(P), 'ShowArray',false,'fs',fs,'c',343, 'LFEchannel',[]);
            h_ls = synthesizeSDMCoeffs(P, DOA, s);
            if ENABLE_PLOTTING==1, PLOT_REPAIR(case_name, h_ls, h_ls, zeros(size(h_ls)), pars, [], 1, path_for_plots); end
            for si=1:length(stimuli)
                yi = fftfilt(h_ls, repmat(srcsig{si},[1 size(hosirr_pars.ls_dirs_deg,1)]));      
                %audiowrite( [path_for_renders_lt filesep room ' ' stimuli{si} ' ' mic_spec{mi}.name ' SDM' '.wav'], yi, fs, 'BitsPerSample', 32);
                audiowrite( [path_for_renders_lt_bin filesep room ' ' stimuli{si} ' ' mic_spec{mi}.name ' SDM'  '.wav'], matrixConvolver(yi, h_ls2bin, size(hrirs,1)), fs, 'BitsPerSample', 32);
            end
        end
    end
end
