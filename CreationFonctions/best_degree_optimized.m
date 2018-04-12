function [Best_poly,Best_denoised_polyN] = best_degree( spc,Best_denoised,B,a,S_noised)
    
    warning('off')

    [mu, sigma] = normfit(S_noised);

    % 4th step : RSP, changing referential : equi -> cheb2 (transformation en chebfun
    % sans dire que c'est un echantillonnage equi)
    f_cheb = chebfun(spc,[B(1) B(end)]); % spectre theorique
    f_cheb_ftt_N = chebfun(Best_denoised,[B(1) B(end)]); % spectre débruite
    
    % calculate Chebyshev points
    N = length(f_cheb_ftt_N);
    Original_Chebyshev_points = chebpts(N,[B(1) B(end)]); 

    parfor i = 1 : N
        x_cheb{i} = chebpts(i,[B(1) B(end)]); % [-1; +1]
    end
    
    % calculate y(N') = f_cheb_ftt(X_cheb defined by degree N')
    parfor ii = 1 : N
        y_Nprim{ii} = f_cheb_ftt_N(x_cheb{1,ii});
    end

    % Transform y_Nprim into a chebfun
    parfor iii = 1 : N
        f_cheb_fft_Nprim{iii} = chebfun(y_Nprim{iii},[B(1) B(end)]);
    end
    
    %%%% 8th step : ask chebfun function f_cheb_ftt_Nprim on originals chebyshev points
    parfor iiii = 1 : N
        y_N_new_fft{iiii} =  f_cheb_fft_Nprim{iiii}(Original_Chebyshev_points);
    end


    %%%% 9th step : calculate residuals
    y_N_original = f_cheb(Original_Chebyshev_points);

    parfor iiiii = 1 : N
        residuals{iiiii} = y_N_new_fft{1,iiiii} - y_N_original;
        [Mu_residuals{iiiii}, Sigma_residuals{iiiii}] = normfit(residuals{iiiii});
    end
    
    %%% Study of the residuals

    parfor i = 1:N
        E(i) = ettest_last(residuals{i},[mu, sigma],a);
    end

    focus = 1 : length(B);

    Extraction = (find(E>=a));

    for k = 1 : length(find(E>=a))
        spc_debruite_poly(:,k) = y_N_new_fft{Extraction(k)};
        Error_pct(k) = 100 * norm(y_N_original(focus)-spc_debruite_poly(focus,k),inf);   
    end
    
    [Y_Best, X_Best] = min(Error_pct);

    Best_poly = Extraction(X_Best);

    Best_denoised_polyN = spc_debruite_poly(:,X_Best);
end

