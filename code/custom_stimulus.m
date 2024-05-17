%
clear all; close all;

project_params = doc_nft_params();
addpath(genpath([project_params.code_fp '\nft']));

in_fp = project_params.data_fp;
fn = 'dataset.mat';
% [fn, in_fp] = uigetfile([in_fp '*.mat'], 'Select Results file', 'MultiSelect','off');

%      resH = '0_20_HEALTHY resampled 624section pruned with ICA'; %m
%      resM = '131_mcs- resampled 1814section pruned with ICA'; %m
   resH = '0_02_HEALTHY resampled -1section pruned with ICA'; %m
   resM = '40_UWS resampled 1394section pruned with ICA'; %f
load([in_fp fn]);
HEALTHYres = Results(strcmp({Results.Name},resH));
MCSres = Results(strcmp({Results.Name},resM));

%setup
% stimArea = 'Cortical';
% phia_inx = 1;
% inhibitoryFactor = 0.7; %inhibitory cortical population has a different influence than the exitator. inhibitoryFactor = 1 however, gives the best spectrum
% stimArea = 'Reticular';
% phia_inx = 2;
stimArea = 'S_Relay';
phia_inx = 3;

%parameters%
phinScaling = 1e0; %phin scaling. adjust to put the system in the correct mode
no_noise_stop_band_flg = false;
use_healthy_P_EXP_flg = false;
instability_type = 'STABLE'; % STABLE X+Y=1 Z=1 Z=0
%
plot_flg = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% usefull definitions

%destabilize the MCS model
switch instability_type
    case 'X+Y=1' %critical point (zero-frequency instability)
        gee_rat = (1-MCSres.NFTparams.xyz(2))/MCSres.NFTparams.xyz(1);
        MCSres.NFTparams.gab(1)=MCSres.NFTparams.gab(1)*gee_rat; %Gee
%         ges_rat = (1-MCSres.NFTparams.xyz(1))/MCSres.NFTparams.xyz(2);
%         MCSres.NFTparams.gab(3)=MCSres.NFTparams.gab(3)*ges_rat; %Ges
    case 'Z=1' %spindle instability
        gsrs_rat = 1/MCSres.NFTparams.xyz(3);
%         MCSres.NFTparams.gab(5)=MCSres.NFTparams.gab(5)*gsrs_rat; %Gsr
        MCSres.NFTparams.gab(8)=MCSres.NFTparams.gab(8)*gsrs_rat; %Grs
    case 'Z=0' %theta instability
        % MCSres.NFTparams.gab(5)=0; %Gsr
        MCSres.NFTparams.gab(8)=0; %Grs  
        %also need to move Y towards -1 by modifying Gesre through MCSres.NFTparams.gab(5) and/or MCSres.NFTparams.gab(7)??
end

fs = project_params.nftsim.fs;
f = MCSres.Spectra.f;
df = f(2)-f(1);
N = round(fs/df/2)*2;
t = (0:N-1)/fs;
f_padded = 0:df:fs/2-df;
w = 2*pi*f_padded';
project_params.psd.window_sec = N/fs;

Gex=1;Giy=1;Grz=1;Gsw=1;
p = MCSres.NFTparams;
Gee=p.gab(1);Gei=p.gab(2);Ges=p.gab(3);Gse=p.gab(4);Gsr=p.gab(5);Gsn=p.gab(6);Gre=p.gab(7);Grs=p.gab(8);
if p.alpha(1)<=0 || ~isfinite(p.alpha(1)) || p.beta(1)<=0 || ~isfinite(p.beta(1)) || size(p.phia,1)~=1 || ...
   any(p.gab == 0) || any(~isfinite(p.gab)) || Gee<0 || Gei>0 || Ges<0 || Gse<0 || Gsr>0 || Gsn<0 || Gre<0 || Grs<0
    warning('PARAMETERS error!');
end

Lw = ((1-i*w/p.alpha(1)).^-1).*((1-i*w/p.beta(1)).^-1);
rho = p.phia(phia_inx)*(1-p.phia(phia_inx)/p.qmax)/p.sigma;
if rho<=0 || ~isfinite(rho)
    error('RHO error!');
end

if MCSres.NFTparams.phin ~= HEALTHYres.NFTparams.phin
    error('PHIN not equal!');
end    
phin = phinScaling*p.phin;
if project_params.nftsim.grid_edge > 1
    deltax = p.Lx/project_params.nftsim.grid_edge;
    sigma = sqrt(2*4*pi^3*fs*phin^2/deltax^2);
else
    sigma = sqrt(2*pi*fs*phin^2);
end
if sigma<=0 || ~isfinite(sigma)
    error('SIGMA error!');
end
HEALTHYres.NFTparams.phin = phin;
MCSres.NFTparams.phin = phin;


%% P_ratio
HEALTHY_P = mean(HEALTHYres.Spectra.P,1);
MCS_P = mean(MCSres.Spectra.P,1);
[~,HEALTHY_P_fit] = HEALTHYres.NFTparams.spectrum(0,0,f,0);
[~,MCS_P_fit] = MCSres.NFTparams.spectrum(0,0,f,0);
% HEALTHY_P_fit = full_gab_p_fit(HEALTHYres)';
% MCS_P_fit = full_gab_p_fit(MCSres)';
% HEALTHY_P_fit = mean(HEALTHYres.Spectra.P_fit,1);
% MCS_P_fit = mean(MCSres.Spectra.P_fit,1);
PHw = zeros(N/2,1); PMw = zeros(N/2,1);
if use_healthy_P_EXP_flg
    PHw(f_padded>=f(1) & f_padded<=f(end)) = HEALTHY_P;
else
    PHw(f_padded>=f(1) & f_padded<=f(end)) = HEALTHY_P_fit;
end
PMw(f_padded>=f(1) & f_padded<=f(end)) = MCS_P_fit;
P_ratio = PHw./PMw;
P_ratio(isnan(P_ratio) | isinf(P_ratio)) = 0;


%% theta, r, dendrites
mcs_phase_w = unifrnd(-pi,pi,N/2,1); %normal random phase
% mcs_phase_w = 0.2+0.41*t(1:N/2)'; %0*unifrnd(-pi,pi,N/2,1);
r = sigma*(1+sqrt(P_ratio));
switch stimArea
    case 'Cortical'
        Giy = inhibitoryFactor*Gex;
        Ew = Gei*(Giy-Gex) + Gex./Lw - Gex*Gsr*Grs.*Lw + (Gex-Giy)*Gei*Gsr*Grs.*Lw.^2;
        theta = mcs_phase_w - angle(Ew) + w*p.t0/2 - pi;
        r = r.*Ges*Gsn./abs(Ew);
        project_params.stim.dendrites.nu = [Gex Giy]/rho;
        project_params.stim.dendrites.alpha = [p.alpha(1) p.alpha(1)];
        project_params.stim.dendrites.beta = [p.beta(1) p.beta(1)];
    case 'Reticular'
        theta = mcs_phase_w - angle(Lw);
        r = r.*Gsn./(-Gsr*Grz*abs(Lw));
        project_params.stim.dendrites.nu = Grz/rho;
        project_params.stim.dendrites.alpha = p.alpha(1);
        project_params.stim.dendrites.beta = p.beta(1);
    case 'S_Relay'
        theta = mcs_phase_w - pi;
        r = r*Gsn/Gsw;
        project_params.stim.dendrites.nu = Gsw/rho;
        project_params.stim.dendrites.alpha = p.alpha(1);
        project_params.stim.dendrites.beta = p.beta(1);
end
theta(1) = 0;
if no_noise_stop_band_flg % disable noise at all the other frequencies. negative values affect only the phase
    r = r.*(P_ratio>0);
end


%% fft, ifft
r = [r; 0; flip(r(2:end))];
theta = [theta; 0; -flip(theta(2:end))];
Sw = r.*exp(i*theta); %pol2cart is also possible
St = ifft(Sw,'symmetric'); 
Sw_gag = fft(St); 
Sw_gag = Sw_gag(1:N/2); %reconstructed


%% represent ifft with trigonometric functions as Fourier series
a0 = Sw(1)/N;
A =  2*real(Sw(2:end))/N;
B = -2*imag(Sw(2:end))/N;
St_Fser = a0;
for k = 1:N/2
    St_Fser = St_Fser + A(k)*cos(2*pi*fs*k*t/N) + B(k)*sin(2*pi*fs*k*t/N);
end
St_Fser = St_Fser';
Sw_Fser_gag = fft(St_Fser); Sw_Fser_gag = Sw_Fser_gag(1:N/2);


%% reconstruct PHw
switch stimArea
    case 'Cortical'
        PMw_stim = PMw.*abs((Ew./(Ges*Gsn*exp(i*w*p.t0/2))) .* Sw_gag/sigma + 1*exp(i*mcs_phase_w)).^2;
        PMw_Fser_stim = PMw.*abs((Ew./(Ges*Gsn*exp(i*w*p.t0/2))) .* Sw_Fser_gag/sigma + 1*exp(i*mcs_phase_w)).^2;
    case 'Reticular'
        PMw_stim = PMw.*abs((Gsr*Grz*Lw/Gsn) .* Sw_gag/sigma + 1*exp(i*mcs_phase_w)).^2;
        PMw_Fser_stim = PMw.*abs((Gsr*Grz*Lw/Gsn) .* Sw_Fser_gag/sigma + 1*exp(i*mcs_phase_w)).^2;
    case 'S_Relay'
        PMw_stim = PMw.*abs((Gsw/Gsn) * Sw_gag/sigma + 1*exp(i*mcs_phase_w)).^2;
        PMw_Fser_stim = PMw.*abs((Gsw/Gsn) * Sw_Fser_gag/sigma + 1*exp(i*mcs_phase_w)).^2;
end

if plot_flg
    figure('Name','Stimulus Reconstruction');
    subplot(2,2,1);plot(t,St, t,St_Fser,'--');xlabel('sec');title('reconstructed stimulus signal');
    legend('ifft', 'fourier series');
    subplot(2,2,2);loglog(f_padded,P_ratio, f_padded,abs(Sw_gag).^2, f_padded,abs(Sw_Fser_gag(1:N/2)).^2,'--');xlabel('hz');title('reconstructed stimulus PSD');
    legend('PSD ratio', 'ifft', 'fourier series'); %ylim([1e-1 1e1]);
    subplot(2,2,3);loglog(f_padded,PHw, f_padded,PMw, f_padded,PMw_stim,'--');xlabel('hz');title('HEALTHY and MCS PSD');
    legend('PHw', 'PMw', 'PMw stim');
    subplot(2,2,4);loglog(f_padded,PHw, f_padded,PMw, f_padded,PMw_Fser_stim,'--');xlabel('hz');title('HEALTHY and MCS PSD - fourier series');
    legend('PHw', 'PMw', 'PMw stim');    
end
%
figure('Name','publish stim time series');
plot(t,St,'linewidth',project_params.grapics.linewidth-2);
ylim([-0.05,0.05]);
ax = gca;
ax.XGrid = "on"; ax.YGrid = "on"; ax.GridColor = project_params.grapics.GridColor; ax.GridAlpha = project_params.grapics.GridAlpha;
ax.XMinorGrid = "on"; ax.YMinorGrid = "on"; ax.MinorGridColor = project_params.grapics.GridColor; ax.MinorGridAlpha = project_params.grapics.GridAlpha;
ax.Box = "off";
ax.FontSize = project_params.grapics.axisTickFntSz;
ax.FontName = project_params.grapics.fontName;
xlabel('sec','FontSize',project_params.grapics.axisLabelFntSz, 'FontName',project_params.grapics.fontName);ylabel('V       ', 'Rotation',0, 'FontSize',project_params.grapics.axisLabelFntSz, 'FontName',project_params.grapics.fontName); 
% title('Stimulus Time Series', 'FontSize',project_params.grapics.titleFntSz, 'FontName',project_params.grapics.fontName);

figure('Name','Publish Spectra Comparison - Analytical');
loglog(f_padded,PHw,'b', f_padded,PMw,'r', f_padded,PMw_stim,'g-.', 'linewidth',project_params.grapics.linewidth);
xlim([0,f(end)+1]);
ax = gca;
ax.FontSize = project_params.grapics.axisTickFntSz;
xlabel('Hz', 'FontSize',project_params.grapics.axisLabelFntSz); ylabel('V^2', 'FontSize',project_params.grapics.axisLabelFntSz); title('Analytical Power Spectra Comparison', 'FontSize',project_params.grapics.titleFntSz);
legend({' Healthy',' DOC',' DOC+Stim'}, 'FontSize',project_params.grapics.axisLabelFntSz);
%

%% save file
fileID = fopen('..\nft\nftsim\stimulus.bin','w');
fwrite(fileID,[N; a0; A; B],'double');
fclose(fileID);


%% simulate

%for plotting purppses only
[~,CZinx] = min([MCSres.chanlocs.radius]);
if size(MCSres.Spectra.P,1) == 1
    HEALTHYres.Spectra.P = repmat(HEALTHYres.Spectra.P,length(HEALTHYres.chanlocs),1);
    MCSres.Spectra.P = repmat(MCSres.Spectra.P,length(MCSres.chanlocs),1);    
end
HEALTHYres.Spectra.P_fit = HEALTHY_P_fit;
MCSres.Spectra.P_fit = MCS_P_fit;
HEALTHYres.Spectra.f_fit = f; 
MCSres.Spectra.f_fit = f;

psd_scaling = 1/sigma^2;
if use_healthy_P_EXP_flg
    HEALTHY_P_plot = HEALTHY_P;
else
    [~, HP_spatial, f_spatial, HP_sim, f_sim] = simulate_and_process(HEALTHYres, project_params, psd_scaling, 0, 'HEALTHY');
    HEALTHY_P_plot = mean(HP_sim,1);
end
[~, MP_spatial,f_spatial, MP_sim, f_sim] = simulate_and_process(MCSres, project_params, psd_scaling, 1, 'MCS');
if use_healthy_P_EXP_flg
    MCS_P_plot = MCS_P;
else
    MCS_P_plot = mean(MP_sim,1);
end

project_params.stim.area = stimArea;
[~, Mstim_P_spatial, f_spatial, Mstim_P_sim, f_sim] = simulate_and_process(MCSres, project_params, psd_scaling, 2, 'MCSstim');
[Mstim_P_spatial,~] = envelope(Mstim_P_spatial,ceil(df/(f_spatial(2)-f_spatial(1))),'peak'); %envelope is needed due to oversampling
Mstim_P_spatial = max(0,Mstim_P_spatial);   

mse_spatial = (mean(Mstim_P_spatial)-mean(HEALTHY_P_plot))^2;
mse_cz = mean((Mstim_P_sim(CZinx,:)-HEALTHYres.Spectra.P(CZinx,:)).^2);
mse_av = mean((Mstim_P_sim-HEALTHYres.Spectra.P).^2,"all");

% plot results
figure('Name','Spectra Comparison');
subplot(3,1,1); loglog(f, HEALTHY_P_plot,'k', f_spatial, MP_spatial,'g', f_spatial, Mstim_P_spatial,'r');
xlabel('Hz'); ylabel('V^2'); title(['Spatial Spectra, mse=' num2str(mse_spatial)]); 
legend('PHw', 'PMw', 'PMw stim'); 
subplot(3,1,2); loglog(f_sim, HEALTHYres.Spectra.P(CZinx,:),'k', f_sim, MP_sim(CZinx,:),'g', f_sim, Mstim_P_sim(CZinx,:),'r');
xlabel('Hz'); ylabel('V^2'); title(['Cz  Spectra, mse=' num2str(mse_cz)]);  
subplot(3,1,3); loglog(f_sim, HEALTHY_P_plot,'k', f_sim, MCS_P_plot,'g', f_sim, mean(Mstim_P_sim,1),'r');
xlabel('Hz'); ylabel('V^2'); title(['Average  Spectra, mse=' num2str(mse_av)]);          
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%
% function P_fit = full_gab_p_fit(SubjRslt)
%     bt_model = bt.model.full_gab();
%     bt_model.p = SubjRslt.NFTparams;
%     bt_model.set_weights_freq(SubjRslt.Spectra.freqBandHz,true);
%     [~,P_fit] = bt_model.objective_with_spectrum(bt_model.params_from_p(SubjRslt.NFTparams), SubjRslt.Spectra.f);
% end
