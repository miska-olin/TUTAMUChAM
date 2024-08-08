


function out=laskeMomentit(NPL,alpha,D1,D2,NLN,cmd,sigma)

out(1)=NPL;

if abs(alpha)<0.001 || abs(alpha+2)<0.001 || abs(alpha+3)<0.001 
    alpha=alpha+0.001;
end

out(2)=D1^2*NPL*(alpha/(alpha+2))*((D2/D1)^(alpha+2)-1)/((D2/D1)^alpha-1);

out(3)=D1^3*NPL*(alpha/(alpha+3))*((D2/D1)^(alpha+3)-1)/((D2/D1)^alpha-1);


out(4)=NLN;


out(5)=  NLN*cmd^2*exp(2*ln2(sigma));

out(6)= NLN*cmd^3*exp(4.5*ln2(sigma));


end