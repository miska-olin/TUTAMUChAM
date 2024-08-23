function out = massGrowthRate(dp,dMolec,temp,visc,diffMolec,CInf,CSat,rho,molarmass,surfaceTension,accommodationCoeff)

    diffP = diff_p(dp,temp,visc);
    activityVapor = 1;
    K = kelvini(surfaceTension,molarmass/6.022e23,temp,rho,dp);
    fuchs = fuchsi_small(temp,dp,dMolec,pi/6*rho*dp^3,molarmass/6.022e23,diffP,diffMolec,accommodationCoeff);

    out=2*pi*(dp+dMolec)*(diffP+diffMolec)*(CInf-activityVapor*CSat*K)*fuchs;

    out=max(out,0);

end

