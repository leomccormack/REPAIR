function PLOT_REPAIR(name, h_lspkr, h_lspkr_dir, h_lspkr_diff, pars, analysis, applyPROC, path_for_plots)

if nargin>7
    if ~isempty(path_for_plots)
        if ~exist(path_for_plots, 'dir'), mkdir(path_for_plots); end
    end
else
    path_for_plots = [];
end

f = 0:pars.fs/pars.winsize:pars.fs/2;
kkk = {'b', 'r', 'm', 'c', 'g', 'k', 'y', 'b', 'r', 'm', 'c', 'g', 'k', 'y', 'b', 'r', 'm', 'c', 'g', 'k', 'y', 'b', 'r', 'm', 'c', 'g', 'k', 'y'};
display_lpf_cutoff_freq_hz = 100;

%% SPATIAL ANALYSIS ESTIMATES
if ~isempty(analysis)
    figure, 
    h=subplot(4,1,1);
    s=mesh(1:size(analysis.K,1), f, analysis.K.');
    s.FaceColor = 'flat';
    set(h,'YScale','log'), 
    ylabel('frequency')
    colormap('jet'), title('number of sources detected'), xlabel('time, hops'), xlim([1 size(analysis.K,1)]), ylim([f(1) f(end)])
    zlabel('K')

    h=subplot(4,1,2);
    s=mesh(1:size(analysis.diffuseness,1), f, analysis.diffuseness.');
    s.FaceColor = 'flat';
    set(h,'YScale','log'), 
    ylabel('frequency')
    colormap('jet'), title('diffuseness'), xlabel('time, hops'), zlim([0 1]), xlim([1 size(analysis.K,1)]), ylim([f(1) f(end)])
    zlabel('diffuseness')

    h=subplot(4,1,3);
    for ii=1:size(analysis.azim,3)
        for jj=1:size(analysis.azim,2)  
            plot3(1:size(analysis.azim,1), ...
                 f(ii)*ones(size(analysis.azim,1),1), 180/pi.*analysis.azim(:,jj,ii), '.', 'Color', kkk{jj} ); hold on,
        end
    end
    set(h,'YScale','log'), 
    ylabel('frequency'), title('azimuth'), xlabel('time, hops'), xlim([1 size(analysis.K,1)]), ylim([f(1) f(end)])
    zlabel('degrees'), zlim([-180 180])
    grid on

    h=subplot(4,1,4);
    for ii=1:size(analysis.elev,3)
        for jj=1:size(analysis.elev,2)  
            plot3(1:size(analysis.elev,1), ...
                 f(ii)*ones(size(analysis.elev,1),1), 180/pi.*analysis.elev(:,jj,ii), '.', 'Color', kkk{jj} ); hold on,
        end
    end
    set(h,'YScale','log'), 
    ylabel('frequency'), title('elevation'), xlabel('time, hops'), xlim([1 size(analysis.K,1)]), ylim([f(1) f(end)])
    zlabel('degrees'), zlim([-90 90])
    grid on

    sgtitle(name)

    if ~isempty(path_for_plots)
        print([path_for_plots  filesep ['SpatialAnalysis ' name]], '-dpng', '-r300');
        %savefig([path_for_plots  filesep name]);
    end
    
    
    figure, 
    h=subplot(4,1,1);
    energy_dB_ref = 10*log10(abs(analysis.energy_in.'));
    energy_min = max(min(energy_dB_ref(:))-6, -60);
    energy_max = max(energy_dB_ref(:))+6;
    s=mesh(1:size(analysis.energy_in,1), f, max(energy_dB_ref,-59.8));
    s.FaceColor = 'flat';
    set(h,'YScale','log'), 
    ylabel('Energy')
    colormap('jet'), title('Input energy'), xlabel('time, hops'), xlim([1 size(analysis.K,1)]), ylim([f(1) f(end)]), zlim([energy_min energy_max])
    zlabel('Energy')

    h=subplot(4,1,2);
    s=mesh(1:size(analysis.energy_ndiff,1), f, max(10*log10(abs(analysis.energy_ndiff.')),-59.8));
    s.FaceColor = 'flat';
    set(h,'YScale','log'), 
    ylabel('frequency')
    colormap('jet'), title('Non-diffuse energy'), xlabel('time, hops'),  xlim([1 size(analysis.K,1)]), ylim([f(1) f(end)]), zlim([energy_min energy_max])
    zlabel('Energy')
    
    h=subplot(4,1,3);
    s=mesh(1:size(analysis.energy_diff,1), f, max(10*log10(abs(analysis.energy_diff.')),-59.8));
    s.FaceColor = 'flat';
    set(h,'YScale','log'), 
    ylabel('frequency')
    colormap('jet'), title('Diffuse energy'), xlabel('time, hops'),  xlim([1 size(analysis.K,1)]), ylim([f(1) f(end)]), zlim([energy_min energy_max])
    zlabel('Energy')
    
    h=subplot(4,1,4);
    s=mesh(1:size(analysis.energy_diff,1), f, max(10*log10(abs(analysis.energy_total.')),-59.8));
    s.FaceColor = 'flat';
    set(h,'YScale','log'), 
    ylabel('frequency')
    colormap('jet'), title('Total output energy'), xlabel('time, hops'),  xlim([1 size(analysis.K,1)]), ylim([f(1) f(end)]), zlim([energy_min energy_max])
    zlabel('Energy')
end

%% OUTPUT LOUSPEAKER RIR PLOTS
figure, 

h_plot = h_lspkr;
if applyPROC 
    [blp,alp]=butter(1, display_lpf_cutoff_freq_hz/(pars.fs/2)); % low-pass filter
    e_ref=filter(blp,alp,h_plot.^2); e_ref=10*log10(e_ref+eps);
    maxx=max(max(e_ref)); e_ref=e_ref-maxx;
    h_plot=max(-59.8,e_ref);
else
    h_plot = max(20*log10(abs(h_plot)), -59.8);
end
subplot(1,3,3), s=mesh(h_plot); 
s.FaceColor = 'flat';
colormap('jet'), xlim([1, size(pars.ls_dirs_deg,1)]), zlim([-60, 10 ])
title('combined'), xlabel('loudspeaker #'), zlabel('energy, dB'), ylabel('time, samples')

h_plot = h_lspkr_dir;
if applyPROC 
    [blp,alp]=butter(1, display_lpf_cutoff_freq_hz/(pars.fs/2)); % low-pass filter
    e_ref=filter(blp,alp,h_plot.^2); e_ref=10*log10(e_ref+eps);
    e_ref=e_ref-maxx;
    h_plot=max(-59.8,e_ref);
else
    h_plot = max(20*log10(abs(h_plot)), -59.8);
end
subplot(1,3,1), s=mesh(h_plot);
s.FaceColor = 'flat';
colormap('jet'), xlim([1, size(pars.ls_dirs_deg,1)]), zlim([-60,10 ])
title('direct stream'), xlabel('loudspeaker #'), zlabel('energy, dB'), ylabel('time, samples')

h_plot = h_lspkr_diff;
if applyPROC 
    [blp,alp]=butter(1, display_lpf_cutoff_freq_hz/(pars.fs/2)); % low-pass filter
    e_ref=filter(blp,alp,h_plot.^2); e_ref=10*log10(e_ref+eps);
    e_ref=e_ref-maxx;
    h_plot=max(-59.8,e_ref);
else
    h_plot = max(20*log10(abs(h_plot)), -59.8);
end
subplot(1,3,2), s=mesh(h_plot); 
s.FaceColor = 'flat';
colormap('jet'), xlim([1, size(pars.ls_dirs_deg,1)]), zlim([-60, 10 ])
title('ambient stream'), xlabel('loudspeaker #'), zlabel('energy, dB'), ylabel('time, samples')


sgtitle(name)

if ~isempty(path_for_plots)
    print([path_for_plots  filesep ['LRIR ' name]], '-dpng', '-r300');
    %savefig([path_for_plots  filesep name]); 
end

end

