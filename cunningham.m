function [ cc ] = cunningham( dp,temp )

    lambda = vapaa_matka(temp);
    
    expo = -0.4985.*dp./lambda;
    kerroin = 2.33+0.966.*exp(expo);
    
    cc = 1+lambda./dp.*kerroin;
    
end

