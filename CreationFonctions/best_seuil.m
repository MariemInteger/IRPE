function [seuil,Best_denoised] = best_seuil(focus,a,spc,spc_N,S_noised)
    warning('off')

    [seuil,E,Spc_D]=denoise_residu_last(spc_N,S_noised,a);
    Extraction = (find(E>=a));
    
    for k = 1 : length(find(E>=a))
        spc_debruite(:,k) = Spc_D{Extraction(k)};
        Error_pct_fft(k) = 100 * norm(spc(focus)-spc_debruite(focus,k),inf);
    end
    
    [Y_Best, X_Best] = min(Error_pct_fft);
    Best_denoised = spc_debruite(:,X_Best);
    Best_seuil_fft = seuil(Extraction(X_Best));
    
    seuil = Best_seuil_fft;
end