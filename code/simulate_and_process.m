%
% simulate, post-process, plot simulation, extract features in segments
%
function [EEGsim, P_spatial, f_spatial, P_sim, f_sim] = simulate_and_process(SubjRslt, project_params, psd_scaling, sim_file_id, title_str)

    [EEGsim, SpatialSpectra, ~, isSimSuccess] = simulate_nft(SubjRslt.NFTparams, SubjRslt.Spectra, project_params, sim_file_id, 0);
    if ~isSimSuccess
        error([title_str ' SIMULATION FAILED!'])
    end
    EEGsim.data = sqrt(psd_scaling)*EEGsim.data; %scaling

    %post-processing
    EEGsim = pop_eegfiltnew(EEGsim, project_params.pipelineParams.passBandHz{1}, project_params.pipelineParams.passBandHz{2}); %high-pass (band pass) filter
%     order=4; cutoffFreq=0.3;
%     [b, a] = butter(order, 2*pi*cutoffFreq/(EEGsim.srate/2), 'high');
%     EEGsim.data = filter(b, a, EEGsim.data, [], 2);
%     EEGsim.data(:,1:250) = EEGsim.data(:,501:750); %remove transient

    %prepare spectra
    P_spatial = SpatialSpectra.P(SpatialSpectra.f>=SubjRslt.Spectra.freqBandHz(1) & SpatialSpectra.f<=SubjRslt.Spectra.freqBandHz(2))';
    f_spatial = SpatialSpectra.f(SpatialSpectra.f>=SubjRslt.Spectra.freqBandHz(1) & SpatialSpectra.f<=SubjRslt.Spectra.freqBandHz(2));
    P_spatial = psd_scaling*P_spatial; %scaling
    
    [P_sim, f_sim] = compute_psd(project_params.nftfit.psdMethod, EEGsim.data, EEGsim.srate, project_params.psd.window_sec, ...
        project_params.psd.overlap_percent, project_params.nftfit.freqBandHz, false);   

end
