function [ la ] = vapaa_matka( temp )

    su = 110.4;
	temp_r = 296.15;
	l_r = 67.3e-9;
	
	eka = temp./temp_r;
	tokayla = 1.0+su/temp_r;
	tokaala = 1.0+su./temp;
	toka = tokayla./tokaala;
	
	la= l_r.*eka.*toka;
    
end

