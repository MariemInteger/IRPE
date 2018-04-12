function [App_cheb,Lpp_cheb,Intensity_cheb] = Cheb_processing(spc_cheb,B_exp)

    [Max,I_max] = max(spc_cheb);
    [Min,I_min] = min(spc_cheb);
    App_cheb = Max-Min;
    Lpp_cheb = (I_min)- (I_max);  
    Intensity_cheb = sum(cumsum(spc_cheb));

    
    
