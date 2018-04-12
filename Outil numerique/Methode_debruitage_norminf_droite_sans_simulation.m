clear all
close all
clc

%% Obtention du bruit typique spectre
cheminScript = mfilename('fullpath');
filenameScript = mfilename;
chemin = strsplit(cheminScript,filenameScript);
chemin = chemin{1};

% [B_exp, spc_exp, Par_exp] = eprload([chemin,'Ongle_keratine']);
% [B_exp, spc_exp, Par_exp] = eprload([chemin,'Simu_RIS5_noisy']);

[B_exp_i, spc_exp_i, Par_exp_i] = eprload([chemin,'imaging_CNRS3']);
tic()
for h=424:size(spc_exp_i,2)
%     clear B_exp spc_exp y_max x_max y_min x_min S_noised mu sigma noise spc_N focus Frequence_coupure E Spc_D droite_debruite Error_pct_fft

%     droite_debruite = zeros(size(spc_exp_i,1),1);
%     Error_pct_fft = zeros(1,1);
    
    B_exp = [];
    spc_exp = [];
    y_max = [];
    x_max = [];
    y_min = []; 
    x_min = []; 
    S_noised = [];
    mu = [];
    sigma = [];
    noise = [];
    spc_N = [];
    focus = [];
    Frequence_coupure = [];
    E = [];
    Spc_D = [];
    droite_debruite = [];
    Error_pct_fft = [];
    
    
    h
    B_exp = B_exp_i{1,1};
    spc_exp = spc_exp_i(:,h);
%     spc_exp = spc_exp./max(spc_exp);

    spc_exp = real(spc_exp);
    spc_exp = basecorr(spc_exp);
    spc_exp = spc_exp - mean(spc_exp(1 : 180,:));

    [y_max, x_max] = max(spc_exp);
    [y_min, x_min] = min(spc_exp);

    %%%% 2nd step : extract noise
    zone= 100:180; %; 
    S_noised = spc_exp(zone);
    % S_noised = (spc_exp(zone));
    % Parametres du bruit
    [mu, sigma] = normfit(S_noised);



    %% Simulation spectre RPE

    x = 1:1:length(spc_exp);
    x = x(:);

    %%%% Droite

    % coeff = 0;
    % spc = @(x) coeff * x + mu ;

    %%%% Gaussienne

    spc = gaussmf(x, [x_max - x_min length(spc_exp)/2]);
    spc = spc .* (y_max - y_min);

    spc = spc(x);

    noise = normrnd(mu,sigma,1,length(x));
    noise = noise(:);

    spc_N = spc + noise;

    focus = 1 : length(spc_N); 



%%% Low Pass Filter %%%

%     a=0.68;
    a = 0.35;
    [Frequence_coupure,E,Spc_D]=denoise_residu_droite(spc_exp,S_noised,a); % l'erreure est tjrs ici
    
    
    for k = 1 : length(Frequence_coupure)
        
        l = Frequence_coupure(k);  
        li = l;
        lf = length(spc_N) - l;
        t = fft(spc_N);
        t(li:lf)=0;

        droite_debruite(:,k) = real(ifft(t));

        Error_pct_fft(k) = 100 * norm(spc-droite_debruite(:,k),inf);
    
    end

    
    [Y_Best, X_Best] = min(Error_pct_fft);

    Best_denoised = droite_debruite(:,X_Best);

    if length(Frequence_coupure)==1
        Best_seuil_fft = Frequence_coupure;
    else
        Best_seuil_fft = Frequence_coupure(X_Best);
    end
    
    spc_exp_iD(:,h) = Spc_D{Best_seuil_fft};
    
    
    
    %%% Polynomial smoothing %%%

    
%     [Degree_smooth,E_smooth,Spc_Smooth]=smooth_residu_droite(spc_exp,S_noised,a);
%     
%     
%     for k = 1 : length(Degree_smooth)
%         
%         
%         l = Degree_smooth(k); 
% 
%         droite_smooth(:,k) = datasmooth(spc_N, l, 'binom'); 
% 
%         Error_pct_smooth(k) = 100 * norm(spc-droite_smooth(:,k),inf);
% 
%     end
% 
%     
%     [Y_Best_smooth, X_Best_smooth] = min(Error_pct_smooth);
% 
%     Best_smooth = droite_smooth(:,X_Best_smooth);
% 
%     Best_seuil_smooth = Degree_smooth(X_Best);

end

toc()

save('denoisedCNRS');
%% Figures %%%%
   
    
% figure()   
% subplot(2,1,1)
% plot(x, [spc_N Best_denoised], 'Linewidth', 1.5)
% xlabel('Magnetic field, G')
% ylabel('dy / dB, a.u')
% legend('Noisy spectrum', 'Denoised spectrum')
% 
% subplot(2,1,2)
% plot(x, [spc Best_denoised], 'Linewidth', 1.5)
% xlabel('Magnetic field, G')
% ylabel('dy / dB, a.u')
% title(['Erreur = ',num2str(Y_Best),'%'])
% legend('Original spectrum', 'Denoised spectrum')
% 
% 
% 
% figure()   
% subplot(2,1,1)
% plot(x, [spc_N Best_smooth], 'Linewidth', 1.5)
% xlabel('a.u')
% ylabel('dy / dB, a.u')
% legend('Noisy spectrum', 'Smooth spectrum')
% 
% subplot(2,1,2)
% plot(x, [spc Best_smooth], 'Linewidth', 1.5)
% xlabel('a.u')
% ylabel('dy / dB, a.u')
% title(['Erreur = ',num2str(Y_Best_smooth),'%'])
% legend('Original spectrum', 'Smooth spectrum')
% 
% 
% figure()   
% subplot(2,1,1)
% plot(B_exp, [spc_exp Spc_D{Best_seuil_fft}], 'Linewidth', 1.5)
% xlabel('Magnetic field, G')
% ylabel('dy / dB, a.u')
% legend('Noisy spectrum', 'Low pass filter')
% 
% subplot(2,1,2)
% plot(B_exp, [spc_exp Spc_Smooth{Best_seuil_smooth}], 'Linewidth', 1.5)
% xlabel('Magnetic field, G')
% ylabel('dy / dB, a.u')
% legend('Original spectrum', 'Smooth spectrum')
