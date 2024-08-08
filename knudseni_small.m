function out=knudseni_small(temp,dp, dg, mp, mg, diff_p, diff_g) 
		
	v_ave_p = velo_ave(temp,mp);
	v_ave_g = velo_ave(temp,mg);
	
	out= 6.0/(dp+dg)*(diff_p+diff_g)/sqrt(v_ave_p*v_ave_p+v_ave_g*v_ave_g);
end

