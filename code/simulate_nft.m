%
% Fit EEG spectra to NFT model
% 
% inputs:
%       params - fited NFT parameters
%       Spectra - computed and fitted Spectra (for plotting)
%       project_params
%       plot_flg
% otputs:
%       EEG - EEGLAB structure
%       SpatialSpectra.P, SpatialSpectra.f
%       central_chan_data - Cz signal
%       isSimSuccess - does the spectrum looks good?
%
function [EEG, SpatialSpectra, central_chan_data, isSimSuccess] = simulate_nft(params, Spectra, project_params, file_id, plot_flg)

%nftsim params
configs_path = append(project_params.code_fp, "\nft\nftsim\configs\");
config_name = 'fit-braintrak-reproduce';
firemode = [0 0 0 0]+0;%1
int_time = project_params.minSectLenSec + 5; %first 5 seconds are transient
waves = [];
ranseed = [];
fprefix = configs_path+config_name;
nftsim_path = [project_params.code_fp '\nft\nftsim'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

EEG = [];
SpatialSpectra = [];
central_chan_data = [];
isSimSuccess = false;

%run NFTsim
nftsim(params, file_id, firemode, int_time, project_params.nftsim.grid_edge,... %generate config file
    project_params.nftsim.fs, waves, ranseed, fprefix, project_params.nftsim.out_dt, project_params.stim.area, project_params.stim.dendrites);
% pause; %edit fit-braintrak-reproduce_00.conf here
nf_struct = run_nftsim_pwd(sprintf('%s_%02d',config_name,file_id), nftsim_path, false);
[grid_data, longside_nodes, shortside_nodes] = nf.grid(nf_struct, 'Propagator.1.phi');
if project_params.nftsim.grid_edge ~= longside_nodes || project_params.nftsim.grid_edge ~= shortside_nodes
    error('different grid size');
end

if project_params.isMEG_flg==false && project_params.nftsim.grid_edge>1
    grid_data = nf.spatial_filter(nf_struct,'Propagator.1.phi');
%     [f_nftsim, P_nftsim] = nf.spectrum(nf_struct, 'Propagator.1.phi', [], true, true); %n_windows=ceil(size(grid_data,3)/project_params.psd.window_sec*nf_struct.deltat)
elseif project_params.isMEG_flg==true && project_params.nftsim.grid_edge==1
    error('can''t disable volume conduction for 1-node configuration in MEG')
end

if project_params.nftsim.grid_edge == 1
    central_chan_data = grid_data;
    [P,f] = pwelch(central_chan_data,[],[],[],1/nf_struct.deltat);
    isSimSuccess = is_sumulation_successfull(f,P,project_params);
    return;
elseif mod(project_params.nftsim.grid_edge,2) == 1
    mid_coord = ceil(project_params.nftsim.grid_edge/2);
else % mod(project_params.nftsim.grid_edge,2) == 0
    mid_coord = [project_params.nftsim.grid_edge/2, project_params.nftsim.grid_edge/2 + 1];
end
% mid_coord = ceil(project_params.nftsim.grid_edge/2); %don't interpolate
central_chan_data = squeeze(mean(grid_data(mid_coord,mid_coord,:),[1,2]));

[f_nftsim, P_nftsim] = nf.spatial_spectrum(nf_struct,'Propagator.1.phi',[],[],1); %n_windows=ceil(size(grid_data,3)/project_params.psd.window_sec*nf_struct.deltat)
isSimSuccess = is_sumulation_successfull(f_nftsim,P_nftsim,project_params);
SpatialSpectra.P = P_nftsim; SpatialSpectra.f = f_nftsim;

%plot
if plot_flg
    P_nftsim = P_nftsim(f_nftsim>=Spectra.freqBandHz(1) & f_nftsim<=Spectra.freqBandHz(2));
    f_nftsim = f_nftsim(f_nftsim>=Spectra.freqBandHz(1) & f_nftsim<=Spectra.freqBandHz(2));
    
    figure;
    subplot(2,1,1);loglog(Spectra.f,mean(Spectra.P,1), Spectra.f_fit,mean(Spectra.P_fit,1), f_nftsim,P_nftsim);xlim(Spectra.freqBandHz);xlabel('Hz');title('fitted and simulated spectra');
    legend('experimental','fitted','simulated');
    P_fit_adjust = Spectra.P_fit * trapz(Spectra.f,mean(Spectra.P,1)) / trapz(Spectra.f_fit,mean(Spectra.P_fit,1));
    P_nftsim_adjust = P_nftsim * trapz(Spectra.f,mean(Spectra.P,1)) / trapz(f_nftsim,P_nftsim);
    subplot(2,1,2);loglog(Spectra.f,mean(Spectra.P,1), Spectra.f_fit,mean(P_fit_adjust,1), f_nftsim,P_nftsim_adjust);xlim(Spectra.freqBandHz);xlabel('Hz');title('fitted and simulated spectra (power adjusted)');
    legend('experimental','fitted','simulated');
    nf.plot_timeseries(nf_struct,  {'Propagator.1.phi'}, {[1:longside_nodes*shortside_nodes]}, false, false); title('generated time series');
end

%EEG struct with simulated time series
EEG = [];
EEG.setname = config_name;
EEG.filename = [config_name '.set'];
EEG.filepath = '';
EEG.srate = round(1/nf_struct.deltat);
EEG.data = reshape(grid_data,[],size(grid_data,3)); %nf.extract(nf_struct, {'Propagator.1.phi'})';
EEG.times = nf_struct.time*1000;
EEG.xmin = nf_struct.time(1);
EEG.xmax = nf_struct.time(end);
EEG.nbchan = size(EEG.data,1);
EEG.pnts = size(EEG.data,2);
EEG.trials = 1;
EEG.ref        = [];
EEG.icawinv    = [];
EEG.icasphere  = [];
EEG.icaweights = [];
EEG.icaact     = [];
EEG.bad_channels = {'a'};

%set electrode locations
x_vals = (0.5:longside_nodes-0.5) * (params.Lx/longside_nodes); %the activity measured at the middle of the node
y_vals = (0.5:shortside_nodes-0.5) * (params.Ly/shortside_nodes);
[X_coord,Y_coord] = meshgrid(x_vals,y_vals);
coord = [X_coord(:) Y_coord(:)]';
theta = 90 + atan2d( (2*params.Lx*coord(2,:) - params.Lx*params.Ly), (2*params.Ly*coord(1,:) - params.Lx*params.Ly) ); %see reverse transition at spatial_gab.m
theta(theta>180) = theta(theta>180)-360;
theta(theta<=-180) = theta(theta<=-180)+360;
radius = sqrt( (2*coord(1,:)/params.Lx - 1).^2 + (2*coord(2,:)/params.Ly - 1).^2 );
for iElec = 1:length(theta)
    loc = struct('labels',['N' num2str(iElec)], 'theta',theta(iElec), 'radius',radius(iElec), 'sph_radius',10); % sph_radius will be updated later at teh code
    EEG.chanlocs(iElec) = loc;
end
EEG = eeg_checkset(EEG);
EEG = pop_select(EEG, 'nochannel', find([EEG.chanlocs.radius] > project_params.head_radius) ); %remove channels that out of head radius

%interpolate with experimental electrodes
EXPchanlocs = readlocs(project_params.electrodes_fn);
EXPchanlocs([EXPchanlocs.radius] > project_params.head_radius) = []; %remove channels that out of head radius
EXPchanlocs = rmfield(EXPchanlocs, setdiff(fieldnames(EXPchanlocs),fieldnames(EEG.chanlocs)));
nbchan = EEG.nbchan;
if isfield(EXPchanlocs,'sph_radius')
    sph_radius = num2cell(ones(1,nbchan)*mean([EXPchanlocs.sph_radius])); %use average sph_radius of the data
    [EEG.chanlocs.sph_radius] = sph_radius{:};
end
EEG = eeg_interp(EEG, EXPchanlocs);
EEG = pop_select( EEG, 'nochannel', 1:nbchan);
EEG = eeg_checkset(EEG);

end

%%%%%%%%%%%%%%%%%%%
function isSimSuccess = is_sumulation_successfull(f,P,project_params)
% does the spectrum looks good? check if the slope is negative

isSimSuccess = false;
logP = log10(P(f>=project_params.nftfit.freqBandHz(1) & f<=project_params.nftfit.freqBandHz(2)));
logF = log10(f(f>=project_params.nftfit.freqBandHz(1) & f<=project_params.nftfit.freqBandHz(2)));
p = polyfit(logF, logP, 1);
if p(1)<0
    isSimSuccess = true;
end

end
