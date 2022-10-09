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

function [lsir, lsir_ndiff, lsir_diff, analysis] = REPAIR(srir, pars, analysis_oracle)
% Reproduction and Parameterisation of Array Impulse Responses (REPAIR) [1]
% -------------------------------------------------------------------------
% 
% DEPENDENCES
%   Spherical-Harmonic-Transform Matlab library
%       https://github.com/polarch/Spherical-Harmonic-Transform
%   Higher-Order-Ambisonics Matlab library
%       https://github.com/polarch/Higher-Order-Ambisonics
%   Vector-Base-Amplitude-Panning
%       https://github.com/polarch/Vector-Base-Amplitude-Panning
%
% INPUT ARGUMENTS 
%   srir                       : input microphone array/spatial RIR; signalLength x nCH
%   pars.grid_svecs            : array steering vectors; nCH x nDirs
%   pars.grid_dirs_xyz         : measurement grid for the array steering vectors (unit-length Cartesian vectors); nDirs x 3
%   pars.grid_dirs_rad         : measurement grid for the array steering vectors (in Spherical coordinates, radians); nDirs x 2
%   pars.grid_weights          : integration weights; nDirs x 1
%   pars.fs                    : sampling rate in Hz
%   pars.SCMavgOption          : options: {'block', 'recur', 'alltime'}
%   pars.SCMavg_coeff          : temporal averaging coefficient, [0..1], if SCMavgOption is set to "recur"
%   pars.SCMavg_Nframes        : number of frames in each averaging block, if SCMavgOption is set to "block"
%   pars.Kestimator            : options: {'SORTE', 'SORTED', 'RECON', 'ORACLE'}
%   pars.DoAestimator          : options: {'MUSIC', 'SRP', 'ORACLE'}
%   pars.winsize               : window size, in time-domain samples
%   pars.freqGrouping          : options: {'broadband', 'octave', 'erb', 'fullres'}
%   pars.streamBalance         : 0: only diffuse stream, 1: both streams are balanced, 2: only direct stream
%   pars.ENABLE_DIFF_WHITENING     : flag (0 or 1), applies an operation that diagonalises the SCMs when under diffuse conditions 
%   pars.ENABLE_COHERENT_FOCUSING  : flag (0 or 1), only used if the steering vectors are frequency dependent, and if there is some band grouping
%   pars.ENABLE_AMBIENT_ENERGY_PRESERVATION   : flag (0 or 1), forces the beamformers used for the ambient stream to be energy-preserving over the sphere
%   pars.decorrelation         : options: {'off', 'convNoise', 'shapedNoise', 'phaseRand', 'covMatch'}
%   pars.beamformerOption      : options: {'pinv', 'MF', 'SD'}
%   pars.ENABLE_QUANTISE_TO_NEAREST_LS  : flag (0 or 1), quantise to nearest loudspeaker instead of using VBAP
%   pars.maxAnaFreq_Hz         : above this frequency, everything is treated as one band
%   pars.ls_dirs_deg           : loudspeaker directions, in degrees; nLS x 2
%   pars.vbapNorm              : 0:reverberant room, ~0.5: dry listening room, 1: anechoic
%   analysis_oracle.K          : (optional) known number of reflections (otherwise this is estimated)
%   analysis_oracle.est_idx    : (optional) known DoA indices (otherwise they are estimated)
%
% OUTPUT ARGUMENTS
%   lsir                       : loudspeaker impulse responses; signalLength x nLS
%   lsir_ndiff                 : non-diffuse stream only; signalLength x nLS
%   lsir_diff                  : diffuse stream only; signalLength x nLS
%   analysis.K                 : (optional) estimated number of reflections
%   analysis.azim              : (optional) estimated reflection azimuths
%   analysis.elev              : (optional) estimated reflection elevations
%                                
% REFERENCES
%   [1] Submitted for review
%
% -------------------------------------------------------------------------
%
%   Leo McCormack, 28/10/2021, leo.mccormack@aalto.fi, with contributions 
%   from Archontis Politis and Nils Meyer-Kahlen
%
% -------------------------------------------------------------------------

maxNumReflections = 8; % An arbitrary limit on the maximum number of simultaneous reflections
nCH = size(srir,2);
nLS = size(pars.ls_dirs_deg,1);
nGrid = size(pars.grid_dirs_rad,1);
maxK = min(floor(nCH/2), maxNumReflections);  

%%% Defaults/Warnings/Errors etc. 
if strcmp(pars.DoAestimator, 'ORACLE') || strcmp(pars.Kestimator, 'ORACLE')
    if nargin<3, error('DoAestimator and/or Kestimator is set to ORACLE, but no Oracle analysis struct is given to REPAIR.'), end
end
if strcmp(pars.DoAestimator, 'ORACLE') && ~strcmp(pars.Kestimator, 'ORACLE')
    error('How would that work?')
end
if (abs(1-sum(pars.grid_weights))>0.001), error('grid_weights must sum to 1'), end
if (~isfield(pars, 'fs')), error('Please specify "fs"'); end
if (~isfield(pars, 'ls_dirs_deg')) 
    error('Please specify "ls_dirs_deg", in degrees')
end 

disp('REPAIR Configuration:'), pars %#ok

if length(size(pars.grid_svecs))==2, BROAD_BAND_SVECS=1; else, BROAD_BAND_SVECS=0; end

%%% Intialisations
fprintf('REPAIR: \n  - Initialising...\n');

lSig = size(srir,1);
winsize = pars.winsize;
fftsize = 2*winsize; % double the window size for FD convolution
hopsize = winsize/2; % half the window size time-resolution
nBins_anl = winsize/2 + 1; % nBins used for analysis 
nBins_syn = fftsize/2 + 1; % nBins used for analysis 
centrefreqs_anl = (0:winsize/2)'*pars.fs/winsize; 

% Frequency grouping 
switch pars.freqGrouping
    case 'broadband'  
        [~,ind_1kHz] = min(abs(centrefreqs_anl-1000)); 
        freqGrpInd = [1 nBins_anl];
    case 'octave'  
        [freqGrpInd, analysis.grpCentreFreq] = findOctavePartitions(centrefreqs_anl, 20e3);
    case 'erb'  
        [freqGrpInd] = findERBpartitions(centrefreqs_anl, 20e3);
    case 'fullres'
        freqGrpInd = [1:nBins_anl nBins_anl]; 
end
nFreqGrps = length(freqGrpInd)-1;

% Replicate grid_svecs per bin for coding convenience
if BROAD_BAND_SVECS, grid_svecs = repmat(pars.grid_svecs, [1 1 nBins_anl]);
else, grid_svecs = pars.grid_svecs;  end

% Diffuse coherence matrix (DCM)
DFCmtx = zeros(nCH, nCH, nBins_anl);
for nb = 1:nBins_anl
    tmp = grid_svecs(:,:,nb);
    DFCmtx(:,:,nb) = tmp*diag(pars.grid_weights)*tmp';  
end 

% Steering vector alignment over frequency
for grp=1:nFreqGrps
    if grp<nFreqGrps, grp_bins = freqGrpInd(grp):freqGrpInd(grp+1)-1; else, grp_bins = freqGrpInd(grp):freqGrpInd(grp+1); end
    if strcmp(pars.freqGrouping, 'broadband'), v0_ind(grp) = ind_1kHz; else, v0_ind(grp) = round(mean(grp_bins));  end
end
if pars.ENABLE_COHERENT_FOCUSING
    sh_order = 25;%floor(sqrt(nCH)-1);
    Y_grid = getSH(sh_order, [pars.grid_dirs_rad(:,1) pi/2-pars.grid_dirs_rad(:,2)], 'real');
    T_WINGS = zeros(nCH,nCH,nBins_anl,nFreqGrps);
    for grp=1:nFreqGrps
        if grp<nFreqGrps, grp_bins = freqGrpInd(grp):freqGrpInd(grp+1)-1; else, grp_bins = freqGrpInd(grp):freqGrpInd(grp+1); end
        % Band/bin to align to:
        for nb=grp_bins 
            if nb ==v0_ind(grp)
                T_WINGS(:,:,nb,grp) = eye(nCH);
            else 
                N_norm_0 = diag(nCH./diag((grid_svecs(:,:,v0_ind(grp))'*grid_svecs(:,:,v0_ind(grp)))+eps)); 
                N_norm_nb = diag(nCH./diag((grid_svecs(:,:,nb)'*grid_svecs(:,:,nb))+eps));  
                V_0 = grid_svecs(:,:,v0_ind(grp)) * N_norm_0 * diag(pars.grid_weights) * Y_grid; 
                V_nb = grid_svecs(:,:,nb) * N_norm_nb * diag(pars.grid_weights) * Y_grid;
                T_WINGS(:,:,nb,grp) =  V_0*pinv(V_nb); % T_0 set to identity in the eval described in Sec 6.D
            end 
       end 
    end
else
    T_WINGS = repmat(eye(nCH),[1 1 nBins_anl,nFreqGrps]);
end
% Note: If in the SH domain, or using full frequency resolution, T_WINGS = Identity (i.e. WINGS is bypassed)

% Frequency averaged diffuse coherence matrix
DFCmtx_grp = zeros(nCH, nCH, nFreqGrps);
for grp=1:nFreqGrps
    if grp<nFreqGrps, grp_bins = freqGrpInd(grp):freqGrpInd(grp+1)-1; else, grp_bins = freqGrpInd(grp):freqGrpInd(grp+1); end
    for nb=grp_bins
        DFCmtx_grp(:,:,grp) = DFCmtx_grp(:,:,grp) + T_WINGS(:,:,nb,grp)*DFCmtx(:,:,nb)*T_WINGS(:,:,nb,grp)';
    end  
end

% For spatial whitening to obtain identity-like SCM structure when under diffuse conditions
if pars.ENABLE_DIFF_WHITENING
    if ~BROAD_BAND_SVECS, reg_par = 0.01; else, reg_par = 0; end
    T_whiten = zeros(nCH, nCH, nFreqGrps);  
    T_unwhiten = zeros(nCH, nCH, nFreqGrps);  
    for grp=1:nFreqGrps
        DFCmtx_norm = DFCmtx_grp(:,:,grp).*(nCH./trace(DFCmtx_grp(:,:,grp)));
        [U,E] = svd(DFCmtx_norm); 
        T_whiten(:,:,grp) = sqrt(pinv(E + reg_par*eye(nCH)))*U';   
        T_unwhiten(:,:,grp) = U*sqrt(E+reg_par*eye(nCH)); 
    end 
else
    T_whiten = repmat(eye(nCH),[1 1 nFreqGrps]);
    T_unwhiten = repmat(eye(nCH),[1 1 nBins_anl]);
end
% Note: If in the SH domain, DCM - T_whiten*DCM*T_whiten' = 0, (i.e. whitening is bypassed)

% DoA estimator intialisations
switch pars.DoAestimator
    case 'MUSIC'
    case 'SRP'
        % Direction-dependent and frequency-dependent diffuse-field EQ to 
        % ensure that the computed pmap has the same energy for all
        % directions if the input is an isotropic diffuse-field, when
        % using arrays with irregular geometry and/or non-uniformly 
        % distributed sensors.
        pmap_diffEQ = zeros(nGrid,nFreqGrps);
        for grp = 1:nFreqGrps
            if grp<nFreqGrps, grp_bins = freqGrpInd(grp):freqGrpInd(grp+1)-1; else, grp_bins = freqGrpInd(grp):freqGrpInd(grp+1); end
            grid_svec_v0 = grid_svecs(:,:,v0_ind(grp));
            tmp = 0;
            for nb=grp_bins
                tmp = tmp + T_WINGS(:,:,nb,grp)*DFCmtx(:,:,nb)*T_WINGS(:,:,nb,grp)';
            end  
            pmap_diffEQ(:,grp) = 1./real(diag(grid_svec_v0'*tmp*grid_svec_v0)+eps);
        end
        % Note: For uniform spherical arrays, or ideal SH input, pmap_diffEQ will be the same for all directions (i.e. bypassed)
    otherwise, error('Unsupported DoA estimator');
end

% Loudspeaker panning gains / finding nearest loudspeakers
if pars.ENABLE_QUANTISE_TO_NEAREST_LS
    nrst_grid2ls = findClosestGridPoints(pars.ls_dirs_deg*pi/180, pars.grid_dirs_rad);
    spat_gains = zeros(size(pars.grid_dirs_rad,1), nLS);
    for i=1:size(pars.grid_dirs_rad,1), spat_gains(i,nrst_grid2ls(i))=1; end
else 
    ls_groups = findLsTriplets(pars.ls_dirs_deg);
    layoutInvMtx = invertLsMtx(pars.ls_dirs_deg, ls_groups);
    spat_gains = vbap(pars.grid_dirs_rad*180/pi, ls_groups, layoutInvMtx);
    % p-value for VBAP normalisation
    if pars.vbapNorm == 0, vbap_pValue = 2*ones(nBins_anl,1);
    else, vbap_pValue = getPvalue(pars.vbapNorm, centrefreqs_anl); end
end

% Diffuse stream rendering intialisations
% Find nearest grid points to the loudspeaker directions
diff_idx = findClosestGridPoints(pars.grid_dirs_rad, pars.ls_dirs_deg*pi/180);
assert(length(diff_idx)==length(unique(diff_idx)), 'No duplicates are permitted');
M_diff = zeros(nLS,nCH,nBins_anl);
diff_eq = zeros(nBins_anl,1);
for nb = 1:nBins_anl 
    Ad = grid_svecs(:,diff_idx,nb); 
    N_norm = (1/nCH).*eye(nLS); %N_norm = diag(1./diag((Ad'*Ad)+eps));
    if pars.ENABLE_AMBIENT_ENERGY_PRESERVATION
        [U,~,V] = svd(Ad);
        if nLS>nCH, A_p = V(:,1:nCH)*U';
        else,       A_p = V*U(:,1:nLS)'; end 
        M_diff(:,:,nb) = sqrt(N_norm)*A_p;
        % Note: If in the SH domain, this reverts to Energy-Preserving Ambisonic Decoding (EPAD)
    else 
        M_diff(:,:,nb) = sqrt(nLS./nCH).*(nCH/nLS).*(N_norm*Ad');
        % Note: If in the SH domain, this reverts to Sampling Ambisonic Decoding (SAD) 
    end
    
    % Diffuse-field equalise 
    diff_eq(nb,1) = sqrt(1./real(trace(DFCmtx(:,:,nb))./nCH+eps)); 
    if nb==2 && (diff_eq(nb,1)>4 || diff_eq(nb,1)<0.25), error('badly scaled steering vectors'), end
    M_diff(:,:,nb) = max(min(diff_eq(nb,1), 4), 0.25).*M_diff(:,:,nb);  % Max +/-12dB boost  
end 

% Covariance matching based decorrelation uses M_diff as the protoype
if strcmp(pars.decorrelation, 'covMatch')
    M_proto = interpolateFilters(M_diff, fftsize); 
    M = zeros(nLS, nLS, nBins_anl);
    Mr = zeros(nLS, nLS, nBins_anl);
end

%%% Pre-processing of SRIR
% transform window (hanning)
x = 0:(winsize-1);
win = sin(x.*(pi/winsize))'.^2;

% zero pad the signal's start and end for the 50% overlap STFT, and account for the temporal averaging:
srir = [zeros(pars.winsize*2, nCH); srir; zeros(fftsize*2, nCH)];
lSig_pad = size(srir,1);
nFrames = ceil((lSig_pad + pars.winsize)/hopsize)+1;

% Pre-compute input spectra
idx = 1;
framecount = 1;
inspec_frame = zeros(nBins_syn, nCH, nFrames);
while idx + fftsize <= lSig_pad 
    insig_win = win*ones(1,nCH) .* srir(idx+(0:winsize-1),:);
    inspec_tmp = fft(insig_win, fftsize);
    inspec_frame(:,:,framecount) = inspec_tmp(1:nBins_syn,:); % keep up to nyquist 
    idx = idx + hopsize;
    framecount = framecount + 1;
end

% Pre-compute time-averaged SCMs (any frequency averaging is done in main loop)
Cxx_frame = zeros(nCH, nCH, nBins_anl, nFrames);
switch pars.SCMavgOption
    case 'block'
        SCMavg_Nframes = pars.SCMavg_Nframes;
        if mod(SCMavg_Nframes,2)==0, SCMavg_Nframes = SCMavg_Nframes+1; end 
        cvxm_frameBuffer = zeros(nBins_anl, nCH, SCMavg_Nframes); 
        for framecount=1:nFrames 
            inspec_anl = inspec_frame(1:fftsize/winsize:end,:,framecount);  
            cvxm_frameBuffer(:,:,2:end) = cvxm_frameBuffer(:,:,1:end-1);
            cvxm_frameBuffer(:,:,1) = inspec_anl;
            cvxm_nFramesHalf = floor(SCMavg_Nframes/2);  
            if framecount>cvxm_nFramesHalf  
                for nb = 1:nBins_anl 
                    in = permute(cvxm_frameBuffer(nb,:,:), [ 2 3 1]);
                    Cxx_frame(:,:,nb,framecount) = in*in';
                end 
            end 
        end 

    case 'recur'
        for framecount=1:nFrames 
            inspec_anl = inspec_frame(1:fftsize/winsize:end,:,framecount); 
            for nb = 1:nBins_anl 
                in = permute(inspec_anl(nb,:,:), [2 3 1]);
                new_Cxx = in*in';
                if framecount==1
                    Cxx_frame(:,:,nb,1) = (1-pars.SCMavg_coeff).*new_Cxx;
                else
                    Cxx_frame(:,:,nb,framecount) = pars.SCMavg_coeff.*Cxx_frame(:,:,nb,framecount-1) + (1-pars.SCMavg_coeff).*new_Cxx;
                end
            end
        end
        
    case 'alltime'
        Cxx_avg = zeros(nCH, nCH, nBins_anl);
        inspec_anl = inspec_frame(1:fftsize/winsize:end,:,:); 
        for nb = 1:nBins_anl 
            in = permute(inspec_anl(nb,:,:), [2 3 1]);
            Cxx_avg(:,:,nb) = in*in';  
        end
        Cxx_frame = repmat(Cxx_avg, [1 1 1 nFrames]);
        
    otherwise
        error('unsuported pars.SCMavgOption')
end

% "Zero-pad" Oracle parameters
if strcmp(pars.DoAestimator, 'ORACLE') || strcmp(pars.Kestimator, 'ORACLE')
    delay_analysis = 2*winsize/hopsize; % Since srir is zero padded at the start by 2*winsize, we'll do the same here
    analysis_oracle.K = [zeros(1,delay_analysis) analysis_oracle.K ];
    analysis_oracle.est_idx = [nan(size(analysis_oracle.est_idx,1),delay_analysis) analysis_oracle.est_idx];
    assert(size(analysis_oracle.K,2) >= nFrames, "Oracle data does not cover at least the input length!")
    assert(size(analysis_oracle.est_idx,2) >= nFrames, "Oracle data does not cover at least the input length!")
end

%%% Main processing loop
idx = 1;
framecount = 1;
progress = 1; 
 
% storage for estimated parameters
if nargout==4
    analysis.K = zeros(ceil(lSig_pad/hopsize),nBins_anl);
    analysis.diffuseness = zeros(ceil(lSig_pad/hopsize),nBins_anl);
    analysis.azim = nan(ceil(lSig_pad/hopsize),maxK,nBins_anl);
    analysis.elev = nan(ceil(lSig_pad/hopsize),maxK,nBins_anl); 
    analysis.energy_in = zeros(ceil(lSig_pad/hopsize),nBins_anl);
    analysis.energy_ndiff = zeros(ceil(lSig_pad/hopsize),nBins_anl);
    analysis.energy_diff = zeros(ceil(lSig_pad/hopsize),nBins_anl);
    analysis.energy_total = zeros(ceil(lSig_pad/hopsize),nBins_anl);
end

% Signal buffers, mixing matrices etc. 
Ms = zeros(nLS, nCH, nBins_anl);
Md = zeros(nLS, nCH, nBins_anl);
lsir_ndiff = zeros(lSig_pad, nLS);
lsir_diff = zeros(lSig_pad, nLS); 
outspec_ndiff = zeros(nBins_syn, nLS);
outspec_diff = zeros(nBins_syn, nLS);

fprintf('  - Rendering: ')
while idx + fftsize <= lSig_pad
    % Load spectra and Cxx for this frame
    inspec_syn = inspec_frame(:,:,framecount); 
    Cxx = Cxx_frame(:,:,:,framecount);
 
    for grp = 1:nFreqGrps
        if grp<nFreqGrps, grp_bins = freqGrpInd(grp):freqGrpInd(grp+1)-1; else, grp_bins = freqGrpInd(grp):freqGrpInd(grp+1); end

        % Average SCMs across frequency  
        Cxx_WINGS = 0;
        for nb = grp_bins 
            Cxx_WINGS = Cxx_WINGS + T_WINGS(:,:,nb,grp)*Cxx(:,:,nb)*T_WINGS(:,:,nb,grp)';
        end
        Cxx_whiten = T_whiten(:,:,grp)*Cxx_WINGS*T_whiten(:,:,grp)';

        % Source number detection 
        [~,S] = sorted_eig(Cxx_whiten);  
        lambda = diag(real(S)); 
        diff_comedie = comedie(lambda);
        K_comedie = floor( (nCH-1)*(diff_comedie) + 1);
        switch pars.Kestimator
            case 'SORTE'
                if nCH>4,  KK = min([SORTE(lambda) maxK]);
                else,      KK = 1; end 
            case 'SORTED'
                if nCH>4,  KK = min([SORTE(lambda) K_comedie  maxK]);
                else,      KK = min([maxK K_comedie]); end 
            case 'RECON'
                KK=0:maxK; 
            case 'ORACLE'
                KK = analysis_oracle.K(1,framecount);  
        end  

        % Source DOA estimation 
        grid_svec_v0 = grid_svecs(:,:,v0_ind(grp));
        j=1; est_idx = {};
        switch pars.DoAestimator
            case 'MUSIC' 
                [V,~] = sorted_eig(Cxx_WINGS);
                for K=KK
                    if K>0 
                        Vn = V(:,K+1:end); 
                        pmap = sphMUSIC(grid_svec_v0, Vn); 
                        est_idx{j} = peakFind2d(pmap, pars.grid_dirs_rad, pars.grid_dirs_xyz, K); 
                    else
                        est_idx{j} = [];
                    end
                    j=j+1;
                end
            case 'SRP' 
                pmap = real(diag(grid_svec_v0'*Cxx_WINGS*grid_svec_v0)).*pmap_diffEQ(:,grp);
                for K=KK
                    if K>0 
                        est_idx{j} = peakFind2d(pmap, pars.grid_dirs_rad, pars.grid_dirs_xyz, K); 
                    else
                        est_idx{j} = []; 
                    end
                    j=j+1;
                end
            case 'ORACLE'
                est_idx{1} = analysis_oracle.est_idx(1:KK,framecount); 
        end 

        % If using the "reconnaissance" approach (i.e. the brute-force reconstruction error based determination of K)
        switch pars.Kestimator
            case {'SORTE', 'SORTED', 'ORACLE'}
                K = KK;
                est_idx = est_idx{1};
            case 'RECON'
                j=1;
                for K=KK 
                    % Obtain mixing matrices
                    if K==0
                        recon_Ms = zeros(nCH); 
                        recon_Md = eye(nCH); 
                    else  
                        Gs = grid_svecs(:,est_idx{j},v0_ind(grp)); 
                        Gd = eye(nCH); 
                        [recon_Ms, recon_Md] = constructMixingMatrices(pars.beamformerOption, grid_svecs(:,est_idx{j},v0_ind(grp)), DFCmtx_grp(:,:,grp), Gs, Gd, pars.streamBalance);
                    end   

                    % Use to reconstruct the input   
                    Cxx_norm = Cxx_WINGS./(trace(Cxx_WINGS)+eps);
                    Cxx_recon_s = recon_Ms*Cxx_norm*recon_Ms';
                    Cxx_recon_d = T_unwhiten(:,:,grp)*(recon_Md*Cxx_norm*recon_Md' .* eye(nCH))*T_unwhiten(:,:,grp)'; % Assuming perfect decorrelation, so more important that channel energies are OK
                    Cxx_recon = Cxx_recon_s + Cxx_recon_d; 
                    recon_err(j, 1) = trace((Cxx_norm-Cxx_recon)*(Cxx_norm-Cxx_recon)');  

                    j=j+1;
                end

                % Determine optimal K
                [~,min_ind] = min(recon_err);
                est_idx = est_idx{min_ind};
                K = KK(min_ind); 
        end  
        
        % CONSTRUCT REPAIR SYNTHESIS MATRICES  
        if K==0 % Bypass direct stream if there are no reflections
            for nb = grp_bins
                Ms(:,:,nb) = zeros(nLS,nCH);
                Md(:,:,nb) = M_diff(:,:,nb);
            end
        else  
            for nb = grp_bins
                % Construct mixing matrices using VBAP gains for target setup
                Gs = spat_gains(est_idx,:).';
                if ~pars.ENABLE_QUANTISE_TO_NEAREST_LS && pars.vbapNorm>0
                    Gs = Gs./(ones(nLS,1) * sum(Gs.^vbap_pValue(nb)).^(1/vbap_pValue(nb)));
                end 
                [Ms(:,:,nb), Md(:,:,nb)] = constructMixingMatrices(pars.beamformerOption, grid_svecs(:,est_idx,nb), DFCmtx(:,:,nb), Gs, M_diff(:,:,nb), pars.streamBalance);
            end
        end

        % for optional plotting
        if nargout==4
            for nb = grp_bins
                analysis.K(framecount,nb) = K;
                analysis.KperGroup(framecount,grp) = K;
                analysis.diffuseness(framecount,nb) = diff_comedie;
                if K>0
                    analysis.azim(framecount,1:K,nb) = pars.grid_dirs_rad(est_idx,1);
                    analysis.elev(framecount,1:K,nb) = pars.grid_dirs_rad(est_idx,2);
                    analysis.azimPerGroup(framecount,1:K,grp) = pars.grid_dirs_rad(est_idx,1);
                    analysis.elevPerGroup(framecount,1:K,grp) = pars.grid_dirs_rad(est_idx,2);
                    
                end
            end
        end
    end
   
    % for optional plotting
    if nargout==4
        for nb=1:nBins_anl
            analysis.energy_in(framecount,nb) = real(trace(Cxx(:,:,nb)));
            analysis.energy_ndiff(framecount,nb) = real(trace(Ms(:,:,nb)*Cxx(:,:,nb)*Ms(:,:,nb)'));
            analysis.energy_diff(framecount,nb) = real(trace( Md(:,:,nb)*Cxx(:,:,nb)*Md(:,:,nb)'  .* eye(nLS) ));
            analysis.energy_total(framecount,nb) = real(trace( Ms(:,:,nb)*Cxx(:,:,nb)*Ms(:,:,nb)' + (Md(:,:,nb)*Cxx(:,:,nb)*Md(:,:,nb)'  .* eye(nLS)) ));
        end
    end

    % REPAIR SYNTHESIS   
    Ms_int = interpolateFilters(Ms, fftsize); 
    Md_int = interpolateFilters(Md, fftsize);
    for nb=1:nBins_syn
        in = permute(inspec_syn(nb,:,:), [2 3 1]);
        outspec_ndiff(nb,:,:) = Ms_int(:,:,nb)*in;
        outspec_diff(nb,:,:) =  Md_int(:,:,nb)*in;
    end 
    
    % Decorrelation options
    switch pars.decorrelation
        case 'convNoise'   % Applied outside of main loop
        case 'shapedNoise' % Applied outside of main loop
        case 'phaseRand'
            randomPhi = rand(size(outspec_diff))*2*pi-pi;
            outspec_diff = abs(outspec_diff) .* exp(1i*randomPhi);
        case 'covMatch'
            % Compute mixing matrices
            for nb=1:nBins_anl 
                Cproto = M_diff(:,:,nb)*Cxx(:,:,nb)*M_diff(:,:,nb)';
                Ctarget = diag(diag(Md(:,:,nb)*Cxx(:,:,nb)*Md(:,:,nb)'));  
                [M(:,:,nb), Cr] = formulate_M_and_Cr(Cproto, Ctarget, eye(nLS), 0, 0.2);
                [Mr(:,:,nb), ~] = formulate_M_and_Cr(diag(real(diag(Cproto))), Cr, eye(nLS), 0, 0.2);  
            end
            M_int = interpolateFilters(M, fftsize); 
            Mr_int = interpolateFilters(Mr, fftsize);
            
            % Apply mixing matrices and override outspec_diff
            for nb=1:nBins_syn
                in = permute(inspec_syn(nb,:,:), [2 3 1]);
                protospec = M_proto(:,:,nb) * in;
                randomPhi = rand(size(protospec))*2*pi-pi;
                protospec_decor = abs(protospec) .* exp(1i*randomPhi); 
                outspec_diff(nb,:,:) = M_int(:,:,nb)*protospec + Mr_int(:,:,nb)*protospec_decor;
            end
    end

    % overlap-add
    lsir_win_ndiff = real(ifft([outspec_ndiff; conj(outspec_ndiff(end-1:-1:2,:))]));
    lsir_ndiff(idx+(0:fftsize-1),:) = lsir_ndiff(idx+(0:fftsize-1),:) + lsir_win_ndiff;
    lsir_win_diff = real(ifft([outspec_diff; conj(outspec_diff(end-1:-1:2,:))]));
    lsir_diff(idx+(0:fftsize-1),:) = lsir_diff(idx+(0:fftsize-1),:) + lsir_win_diff;

    % advance sample pointer
    idx = idx + hopsize;
    framecount = framecount + 1;
    if framecount >= floor(nFrames/10*progress)  
        fprintf('*');
        progress=progress+1; 
    end  
end
fprintf(' done!\n')
    
%%% Post-Processing of LSIR
fprintf('  - Post-Processing...\n');
% Remove delay caused by the filter interpolation of gains and circular shift
delay_sig = winsize*2+hopsize+1;
lsir_ndiff = lsir_ndiff(delay_sig:delay_sig+lSig-1,:);
lsir_diff =  lsir_diff (delay_sig:delay_sig+lSig-1,:);

% Remove the analysis frame results that occured before the first input sample actually encounted any analysis
if nargout==4
    delay_analysis = 2*winsize/hopsize;
    analysis_len = floor(lSig/(hopsize));
    analysis.K = analysis.K(delay_analysis:delay_analysis+analysis_len-1,:);
    analysis.diffuseness = analysis.diffuseness(delay_analysis:delay_analysis+analysis_len-1,:);
    analysis.azim = analysis.azim(delay_analysis:delay_analysis+analysis_len-1,:,:);
    analysis.elev = analysis.elev(delay_analysis:delay_analysis+analysis_len-1,:,:);
    analysis.energy_in = analysis.energy_in(delay_analysis:delay_analysis+analysis_len-1,:);
    analysis.energy_ndiff = analysis.energy_ndiff(delay_analysis:delay_analysis+analysis_len-1,:);
    analysis.energy_diff = analysis.energy_diff(delay_analysis:delay_analysis+analysis_len-1,:);
    analysis.energy_total = analysis.energy_total(delay_analysis:delay_analysis+analysis_len-1,:);
end

% Apply decorrelation
switch pars.decorrelation
    case 'convNoise' 
        % Apply convolution decorrelation to diffuse stream
        % we want to apply just enough noise-based reverberation as
        % to suitably decorrelate the signals, but not change the captured room
        % characteristics too much. T60s of a very, very dry room should suffice for
        % this task:
        t60 = [0.07 0.07 0.06 0.04 0.02 0.01];
        fcentre = [250 500 1e3 2e3 4e3];
        fcutoffs = fcentre/sqrt(2);
        randmat =  synthesizeNoiseReverb(nLS, pars.fs, t60, fcutoffs, 1);
        lsir_diff = fftfilt(randmat, lsir_diff);
        
    case 'shapedNoise' 
        fcentre = 10.^(0.1.*[19:43]);
        fcutoffs = fcentre * 10^0.05; 
        lsir_diff = synthesizeShapedNoise(lsir_diff, fcutoffs, pars.fs);
        
    case 'phaseRand' % Applied within the main loop
    case 'covMatch'  % Applied within the main loop
end

lsir = lsir_ndiff+lsir_diff;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Ms, Md] = constructMixingMatrices(beamformerOption, As, DFCmtx, Gs, Gd, streamBalance)
    nCH = size(As,1);
    K = size(As,2); 
    
    % Source stream beamforming weights, Ds
    switch beamformerOption
        case 'pinv'
            Ds = pinv(As); 
            % Trivia: if K==1, then this reverts to 'MF'
        case 'MF'
            As_n = As*diag(1./(diag(As'*As)+eps));
            Ds = As_n';
            % Trivia: if SHD, then this reverts to hyper-cardioid beamformers
        case 'SD'
            beta = 0.01;
            DFCmtx_dl = DFCmtx + (trace(DFCmtx)+0.000001).*eye(nCH).*beta;
            for k=1:K
                Ds(k,:) = (As(:,k)'/DFCmtx_dl*As(:,k))\(As(:,k)'/DFCmtx_dl);  
            end  
            % Trivia: if SHD, then this reverts to 'MF' (i.e. hyper-cardioid beamformers)
    end 
    
    % Ambient stream beamforming weights, Dd
    Dd = eye(nCH) -  As*Ds;  
    
    % Optional stream balance manipulations
    if streamBalance<1
        a = streamBalance;
        b = 1;
    elseif (streamBalance>=1) && (streamBalance<=2)
        a = 1;
        b = streamBalance;
    end  
    
    % Mixing matrices for direct and ambient streams (Ms & Md)
    Ms = a*Gs*Ds;
    Md = b*Gd*Dd; 
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rir_filt = synthesizeNoiseReverb(nCH, fs, t60, fcutoffs, FLATTEN)
%NOISEVERB Simulates a quick and dirty exponential decay reverb tail
%
% order:    HOA order
% fs:       sample rate
% t60:      reverberation times in different bands
% fc:       cutoff frequencies of reverberation time bands (1 more than t60)
%
%   Archontis Politis, 12/06/2018
%   archontis.politis@aalto.fi

    if nargin<5, FLATTEN = 0; end
    
    % number of HOA channels
    nSH = nCH;
    % number of frequency bands
    nBands = length(t60);
    % decay constants
    alpha = 3*log(10)./t60;
    % length of RIR
    %lFilt = ceil(max(t60)*fs);
    t = (0:1/fs:max(t60)-1/fs)';
    lFilt = length(t);
    % generate envelopes
    env = exp(-t*alpha);
    % generate RIRs
    rir = randn(lFilt, nSH, nBands);
    for k = 1:nBands
        rir(:, :, k) = rir(:,:,k).*(env(:,k)*ones(1,nSH));
    end
    % get filterbank IRs for each band
    filterOrder = 10000;
    h_filt = filterbank(fcutoffs, filterOrder, fs);
    % filter rirs
    rir_filt = zeros(lFilt+ceil(filterOrder/2), nSH);
    for n = 1:nSH
        h_temp = [squeeze(rir(:,n,:)); zeros(ceil(filterOrder/2), nBands)];
        rir_filt(:, n) = sum(fftfilt(h_filt, h_temp), 2);
    end
                            
    if FLATTEN, rir_filt = equalizeMinphase(rir_filt); end
    
    rir_filt = rir_filt(filterOrder/2+1:end,:); % remove delay
end

function ls_diff_shaped_noise = synthesizeShapedNoise(ls_diff, fcutoffs, fs)
    
nLS = size(ls_diff, 2);
nBands = length(fcutoffs)+1;

winSize_samp = 8000;
%win = hann(winSize_samp);
%win = win(winSize_samp/2:end);

% Exponential window
tau = winSize_samp / 6;
win = exp(-(1:winSize_samp)/tau);
win = win / sum(win);

% get filterbank IRs for each band
filterOrder = 70000;
h_filt = filterbank(fcutoffs, filterOrder, fs);
irLen = size(ls_diff, 1);

% zero padding
ls_diff = [ls_diff; zeros(filterOrder, nLS)];

% init
ls_diff_shaped_noise = zeros(irLen + filterOrder, nLS);
rir_filt = zeros(size(ls_diff, 1),  nBands, nLS);

for n = 1:nLS
    
    % analyse in freq bands
    rir_filt(:, :, n) = fftfilt(h_filt, ls_diff(:, n));
    
    % get the envelop in each band
    %env(:, :, n) = real(fftfilt(win, abs(rir_filt(:, :, n))));
    env(:, :, n) = real(fftfilt(win, (rir_filt(:, :, n))));
   
    % synthesize noise
    noise = randn(irLen+filterOrder, 1);
    
    % weight noise with envelope in each band
    shaped_noise = noise .* env(:, :, n);
    
    % filter
    shaped_noise_filt = fftfilt(h_filt, shaped_noise);
    
    ls_diff_shaped_noise(:, n) = sum(shaped_noise_filt, 2);
end

% shift back because filters are zero phase
ls_diff_shaped_noise = ls_diff_shaped_noise(filterOrder+1:end, :);

% normalize
ls_diff_shaped_noise = ls_diff_shaped_noise .* rms(ls_diff) ./ rms(ls_diff_shaped_noise);

%     figure
%     subplot(3, 1, 1)
%     plot(ls_diff)
%     subplot(3, 1, 2)
%     plot(ls_diff_shaped_noise)
%     subplot(3, 1, 3)
%     plot(squeeze(sum(env, 2)))
%     
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function h_filt = filterbank(fcutoffs, filterOrder, fs)
% fc:   the cutoff frequencies of the bands
% Nord: order of hte FIR filter
%
%   Archontis Politis, 12/06/2018
%   archontis.politis@aalto.fi

    if length(fcutoffs) == 0 %#ok
        h_filt = 1;

    elseif length(fcutoffs) == 1
        h_filt = zeros(filterOrder+1, 2);

        % lowpass
        f_ll = fcutoffs(1);
        w_ll = f_ll/(fs/2);
        h_filt(:, 1) = fir1(filterOrder, w_ll);
        % highpass
        f_hh = fcutoffs(2);
        w_hh = f_hh/(fs/2);
        h_filt(:, 2) = fir1(filterOrder, w_hh, 'high');

    else
        Nfilt = length(fcutoffs)+1;
        h_filt = zeros(filterOrder+1, Nfilt);

        % lowpass
        f_ll = fcutoffs(1);
        w_ll = f_ll/(fs/2);
        h_filt(:, 1) = fir1(filterOrder, w_ll);
        % highpass
        f_hh = fcutoffs(end);
        w_hh = f_hh/(fs/2);
        h_filt(:, end) = fir1(filterOrder, w_hh, 'high');
        % bandpass
        for k = 1:Nfilt-2
            fl = fcutoffs(k);
            fh = fcutoffs(k+1);
            wl = fl/(fs/2);
            wh = fh/(fs/2);
            w = [wl wh];
            h_filt(:, k+1) = fir1(filterOrder, w, 'bandpass');
        end
    end

end

function psi = comedie(lambda)
% Implementation based on the COMEDIE estimator as described in:
%   Epain, N. and Jin, C.T., 2016. Spherical Harmonic Signal Covariance and
%   Sound Field Diffuseness. IEEE/ACM Transactions on Audio, Speech, and 
%   Language Processing, 24(10), pp.1796-1807.
% 
nCH = length(lambda);
if all(lambda==0)
    psi = 0;
else
    g_0 = 2*(nCH-1);
    mean_ev = sum(lambda)/nCH;
    g = (1/mean_ev)*sum(abs(lambda-mean_ev));
    psi = 1-g/g_0;
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rir_filt_flat = equalizeMinphase(rir_filt)
%MAKEFLATVERB Makes the decaying noise spectrally flat
%
%   Archontis Politis, 12/06/2018
%   archontis.politis@aalto.fi

Nrir = size(rir_filt,2);
for n=1:Nrir
    % equalise TDI by its minimum phase form to unity magnitude response
    tdi_f = fft(rir_filt(:,n));
    tdi_min_f = exp(conj(hilbert(log(abs(tdi_f)))));
    tdi_eq = real(ifft(tdi_f./tdi_min_f));
    rir_filt_flat(:,n) = tdi_eq;
end

end

function [est_idx, est_dirs] = peakFind2d(pmap, grid_dirs, grid_xyz, nPeaks)

kappa  = 50;
P_minus_peak = pmap;
est_dirs = zeros(nPeaks, 2);
est_idx = zeros(nPeaks, 1);
for k = 1:nPeaks
    [~, peak_idx] = max(P_minus_peak);
    est_dirs(k,:) = grid_dirs(peak_idx,:);
    est_idx(k) = peak_idx;
    VM_mean = grid_xyz(peak_idx,:); % orientation of VM distribution
    VM_mask = kappa/(2*pi*exp(kappa)-exp(-kappa)) * exp(kappa*grid_xyz*VM_mean'); % VM distribution
    VM_mask = 1./(0.00001+VM_mask); % inverse VM distribution
    P_minus_peak = P_minus_peak.*VM_mask;
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function M_interp = interpolateFilters(M, fftsize)
% M filter matrix with y = M*x, size NxLxK,
%   N output channels, M input channels, K frequency bins
%
%   Archontis Politis, 12/06/2018
%   archontis.politis@aalto.fi

winsize = 2*(size(M,3)-1);
M_conj = conj(M(:,:,end-1:-1:2));
M_ifft = ifft(cat(3, M, M_conj), [], 3);
M_ifft = M_ifft(:,:, [(winsize/2+1:end) (1:winsize/2)]); % flip
M_interp = fft(M_ifft, fftsize, 3); % interpolate to fftsize
M_interp = M_interp(:,:,1:fftsize/2+1); % keep up to nyquist
end
 
function P_music = sphMUSIC(A_grid, Vn)
% SPHMUSIC DOA estimator based on MUltiple SIgnal Classification (MUSIC)
VnA = Vn'*A_grid;
P_music = 1./sum(conj(VnA).*VnA);
P_music = P_music.';

end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [K, obfunc] = SORTE(lambda)
% SORTE Second ORder sTatistic of Eigenvalues estimator of source number
% Implementation based on the SORTE estimator as described in:
%   He, Z., Cichocki, A., Xie, S., & Choi, K. (2010). 
%   Detecting the number of clusters in n-way probabilistic clustering. 
%   IEEE Transactions on Pattern Analysis and Machine Intelligence, 32(11), 
%   pp.2006-2021.
%
% INPUT ARGUMENTS
% lambda    % vector of eigenvalues of the spatial covariance matrix of
%             ambisonic or (diffuse-whitened) microphone array signals
%
% OUTPUT ARGUMENTS
% K         % the number of detected sources in the mixtures
% obfunc    % the values of the SORTE criteria for the different source
%             numbers - its minimum indicates theestimated source number
%
% Author:   Archontis Politis (archontis.politis@gmail.com)
% Copyright (C) 2021 - Archontis Politis
%
if all(lambda==0)
    K = 0;
else
    N = length(lambda);
    Dlambda = lambda(1:end-1) - lambda(2:end);
    for k=1:N-1
        meanDlambda = 1/(N-k)*sum(Dlambda(k:N-1));
        sigma2(k) = 1/(N-k)*sum( (Dlambda(k:N-1) - meanDlambda).^2 );
    end
    for k=1:N-2
        if sigma2(k)>0, obfunc(k) = sigma2(k+1)/sigma2(k);
        elseif sigma2(k) == 0, obfunc(k) = Inf;
        end
    end
    obfunc(end) = Inf;
    [~,K] = min(obfunc);    
end

end

function [oct_idx, oct_freqs] = findOctavePartitions(bandfreqs, maxFreqLim)
if (nargin<2), maxFreqLim = 20000; end
if (maxFreqLim>20000), maxFreqLim = 20000; end
oct_freqs = 250;
oct_idx = 1;

counter = 1;
while oct_freqs(counter)*2 < maxFreqLim
    oct_freqs(counter+1) = oct_freqs(counter) * 2;            % Upper frequency limit of band
    % find closest band frequency as upper partition limit
    [~, oct_idx(counter+1)] = min(abs(sqrt(1/2)*oct_freqs(counter+1) - bandfreqs));
    % force at least one bin bands (for low frequencies)
    if (oct_idx(counter+1) == oct_idx(counter)), oct_idx(counter+1) = oct_idx(counter+1)+1; end
    % erb_freqs(counter+1) = bandfreqs(erb_idx(counter+1));          % quantize new band limit to bins

    counter = counter+1;
end
% last limit set at last band
oct_freqs(counter+1) = bandfreqs(end);
oct_idx(counter+1) = length(bandfreqs);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [M, Cr] = formulate_M_and_Cr(Cx, Cy, Q, flag, reg)
% FORMULATE_M_AND_CR Computes optimal mixing matrix based on input-output
%                    covariance matrix
%
% J.Vilkamo's "Covariance Domain Framework for Spatial Audio Processing".
%
    if nargin == 4
        reg=0.2;
    end
    % flag = 0: Expect usage of residuals
    % flag = 1: Fix energies instead
    lambda=eye(length(Cy),length(Cx));

    % Decomposition of Cy
    [U_Cy,S_Cy]=svd(Cy);
    Ky=U_Cy*sqrt(S_Cy);

    % Decomposition of Cx
    [U_Cx,S_Cx]=svd(Cx);
    Kx=U_Cx*sqrt(S_Cx);

    %SVD of Kx
    Ux=U_Cx;
    Sx=sqrt(S_Cx);
    % Vx = identity matrix

    % Regularization Sx
    Sx_diag=diag(Sx);
    limit=max(Sx_diag)*reg+1e-20;
    Sx_reg_diag=max(Sx_diag,limit);

    % Formulate regularized Kx^-1
    Kx_reg_inverse=diag(1./Sx_reg_diag)*Ux';

    % Formulate normalization matrix G_hat
    Cy_hat_diag=diag(Q*Cx*Q');
    limit=max(Cy_hat_diag)*0.001+1e-20;
    Cy_hat_diag=real(max(Cy_hat_diag,limit));
    G_hat=diag(real(sqrt(diag(Cy)./Cy_hat_diag)));

    % Formulate optimal P
    [U,~,V]=svd(Kx'*Q'*G_hat'*Ky);
    P=V*lambda*U';

    % Formulate M
    M=Ky*P*Kx_reg_inverse;

    % Formulate residual covariance matrix
    Cy_tilde = M*Cx*M';
    Cr=Cy-Cy_tilde;

    % Use energy compensation instead of residuals
    if flag==1
        adjustment=diag(Cy)./diag(Cy_tilde + 1e-20);
        G=diag(sqrt(adjustment));
        M=G*M;
        Cr='unnecessary';
    end

end

function p = getPvalue(DTT,f)
% f:        frequency vector
% DTT = 0:  normal room (p=2)
% DTT = 1:  anechoic room (p<2) 
    a1 = 0.00045;
    a2 = 0.000085;
    p0 = 1.5 - 0.5 * cos(4.7*tanh(a1*f)).*max(0,1-a2*f);
    p = (p0-2)*sqrt(DTT)+2;

end
