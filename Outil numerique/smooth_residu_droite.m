function [Degree_smooth,E_smooth,Spc_Smooth]=smooth_residu_droite(spc_exp,S_noised,a)

N = length(spc_exp);

Degree = 1:1:200; 

[mu, sigma] = normfit(S_noised);

j=1;

while Degree(j) < 200
    
    l = Degree(j); 
    spc_smooth = datasmooth(spc_exp, l, 'binom');   
    
    R_smooth = spc_smooth - spc_exp;
    E_smooth(j) = ettest_last(R_smooth,[mu, sigma],a);
    
    Spc_Smooth{j} = spc_smooth;
    
    j=j+1;
end

Degree_smooth = (find(E_smooth>=a));

end