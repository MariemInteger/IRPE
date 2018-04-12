function [App_man,Lpp_man,Intensity_man] = manual_processing(spc_man,B)


    [Max,I_max] = max(spc_man);
    [Min,I_min] = min(spc_man);
    App_man = Max-Min;
    Lpp_man = B(I_min)-B(I_max);  
	S1 = ((B(end)-B(1))/length(spc_man))*cumsum(spc_man);
    Intensity_man = ((B(end)-B(1))/length(spc_man'))* sum(S1);
    
