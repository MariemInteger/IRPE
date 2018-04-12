clear all
close all
clc
warning('off')

% Chemin
cheminScript = mfilename('fullpath');
filenameScript = mfilename;
chemin = strsplit(cheminScript,filenameScript);
chemin = chemin{1};

%% Loading spectra
% Loading experimental spectra of RIS5
[B_exp, spc_exp, Par_exp] = eprload([chemin,'Data/Simu_RIS5_noisy']);
spc_exp = real(spc_exp);
% Loading experimental spectra of empty tube
[B_noised, spc_noised, Par_noised] = eprload([chemin,'Data/Simu_tube_noisy']);
spc_noised = real(spc_noised);

%% 2nd step : extract noise from empty tube
S_noised = basecorr(spc_noised);
[mu, sigma] = normfit(S_noised);

%% Simulation spectre RPE
Sys0.S = 1/2;
Sys0.g = [2.0060 2.0049 2.0022];
Sys0.lw = [0.5956 0.8106];

freqexp1=Par_exp.MWFQ/1e9; 
deltaB1=28;
B_exp=(B_exp+deltaB1)/10; 
Exp0=struct('Range',[B_exp(1) B_exp(end)],'mwFreq',freqexp1);
Exp0.ModAmp = 0.5;
Exp0.nPoints = 1760;
Exp0.Harmonic = 1;

[B,spc]=pepper(Sys0,Exp0);
B = B(:);
spc = spc(:);
spc=rescale(spc, 'minmax'); 
spc = spc .* abs(spc_exp(1110) - spc_exp(926));
spc = basecorr(spc);

%% Génération de N spectres bruités
N = 1;
spc_N = normrnd(mu,sigma,N,length(spc))+spc';
spc_N = spc_N';

%% Initialisation de parametres
focus = 660 : 1310; % zone du signal qui nous interesse
a=0.68; % pourcentage de validation d'appartenance residus dans le bruit

%% Finding the best seuil and Denoising
tic()
for i=1:N
    [seuil(i),Best_denoised(:,i)] = best_seuil(focus,a,spc,spc_N(:,i),S_noised);
end

%% Finding the best chebfun's degree of filtered spectra
for i=1:N
    [Best_poly(i),Best_denoised_polyN(:,i)] = best_degree(spc_N(:,i),Best_denoised(:,i),B,a,S_noised);
end
toc()