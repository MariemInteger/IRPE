function [seuil,E,Spc_D]=denoise_residu_last(spc_N,S_noised,a)

N = length(spc_N);
spc_N_fft = fft(spc_N);
%seuil = 200:-1:1; 
seuil = 1:1:200; 
Test=0;
[mu, sigma] = normfit(S_noised);

j=1;

while seuil(j) < 200
    j
    l = seuil(j);  
    li = ceil(l);
    lf = N-l;
    t = spc_N_fft;
    t(li:lf)=0;
    
    spc_denoised = real(ifft(t));
    
    R = spc_denoised - spc_N;
    E(j) = ettest_last(R,[mu, sigma],a);

    %seuil(j) = seuil(j) + 1;
    
    Spc_D{j} = spc_denoised;
    
    j=j+1;
    end
end