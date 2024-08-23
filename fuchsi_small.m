function b=fuchsi_small(temp,dp, dg, mp, mg, diff_p, diff_g, accommodationCoeff) 
	
kn = knudseni_small(temp,dp,dg,mp,mg,diff_p,diff_g);
up = 0.75*accommodationCoeff*(1.0+kn);
down = kn*kn+kn+0.283*kn*accommodationCoeff+0.75*accommodationCoeff;
b = up/down;

b=rajat(0,b,1);

	
end

