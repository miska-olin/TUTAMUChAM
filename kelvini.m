function kelvin=kelvini(surftens, molecmass, temp, rho, dp) 
    boltz = 1.381e-23;
    kelvin = exp(4*surftens*molecmass/(boltz*temp*rho*dp));

end

