function out =  simulationFS(p)
%SIMULATIONFS simulates the FS model using input parameter set 'p'    

    nSec = n_sec(p.model);
    Dp_edges = logspace(log10(p.dCluster),log10(p.highestDiameter),nSec+1);
    logEdges = log10(Dp_edges);
    logCenters = logEdges(1:end-1) + diff(logEdges)./2;
    Dp_centers = 10.^(logCenters); 
    p.DeltaDp = diff(Dp_edges);
    p.Dp_centers = Dp_centers;

    p.kk = zeros(nSec,nSec); % Preallocate coagulation coefficients
    if p.coag
        for i = 1:nSec
            p.kk(i,:) = 2*pi*coag_kernel_dp(Dp_centers(i),Dp_centers,p.T,p.rho)*1e6;
        end
    end  

    % force number of vapors to 0 if GR is not modelled
    if p.GRModel == 0
        p.nVapor = 0;
    end

    % Calculating coagulation sink factor
    p.coagSinkFactor = coag_kernel_dp(5e-9,p.coagSinkCMD,p.T,p.rho)*2*pi*5e-9^(-p.coagSinkExponent);

    % Total simulation time
    p.totalTime = p.timeVec(end)-p.timeVec(1);

    % The method how initial particle size distribution is given
    switch(p.initialPSDMethod)
        case 'parameters'
            initialN = (n_j_powerlaw(p.initialPSDParameters(1),p.initialPSDParameters(2),p.dCluster,p.initialPSDParameters(3),Dp_centers) +...
                n_j(Dp_centers,p.initialPSDParameters(4),p.initialPSDParameters(5),p.initialPSDParameters(6)))*ln(Dp_edges(2)/Dp_edges(1));
        case 'vector'
            initialN = p.initialMomentVec;
        otherwise
            error('Initial PSD method needs to be ''parameters'' or ''vector'' for FS models')
    end

    % Correct initial particle size distribution based on the sampling line penetration efficiency
    initialN = initialN/p.particleNumberPenetration; 

    % Include vapor concentrations in initialN if applicable
    if p.GRModel == 1
        initialN = [initialN p.vaporConc0'];
    end

    % Options to ODE solver
    options = odeset('RelTol',p.relativeTolerance,'nonnegative',1:(nSec+p.nVapor),'InitialStep',1,'MaxStep',1000,'stats','off');
    
    if p.plotWaitbarDuringSim
        hWait = waitbar(0,'Simulating...');
    end
    
    % Solving ODE
    [t,Y] = feval(p.solverName,@modelFunc,p.timeVec,initialN,options,p);

    if p.plotWaitbarDuringSim
        close(hWait);
    end
    
    disp('Simulation performed')
    disp('Making output struct...')

    % extract vapor concentrations from Y matrix
    vaporConc = Y(:,(end-p.nVapor+1):end);

    % remove vapor concentrations from Y matrix
    Y = Y(:,1:end-p.nVapor);

    % Correct the modelled particle size distributions based on the sampling line penetration efficiency
    Y = Y*p.particleNumberPenetration;

    % Making the output struct 'out'
    out.t = t; 
    out.Y = Y;
    out.p = p;

    % Zeros
    out.N = zeros(length(t),1);
    out.M_2 = out.N;
    out.M_3 = out.N;
    out.GMD = out.N;
    out.GSD = out.N;

    % Update values
    for i=1:length(t)
        out.N(i) = sum(Y(i,:));
        out.M_2(i) = sum(Y(i,:).*(p.Dp_centers*1e9).^2)*pi; % Surface area (nm2/cm3)
        out.M_3(i) = sum(Y(i,:).*p.Dp_centers.^3)*pi/6*p.rho*1e9*1e6; % Mass (ug/m3)
        lngmd = sum(Y(i,:).*log(p.Dp_centers))/sum(Y(i,:));
        ln2gsd = sum(Y(i,:).*(log(p.Dp_centers)-lngmd).^2)/sum(Y(i,:));
        out.GMD(i) = exp(lngmd);
        out.GSD(i) = exp(sqrt(ln2gsd));
    end

    % Removing a large interpolation table from the output struct
    out.p=rmfield(out.p,'alphaTaulukko');
    out.p=rmfield(out.p,'dTaulukko');

    out.vaporConc = vaporConc;

end



function dy = modelFunc(t,y,param)
%MODELFUNC is the internal function containing the FS model equations

    if param.plotWaitbarDuringSim
        waitbar(t/param.totalTime);
    end

    % extract vapor concentrations from y vector
    vaporConc = y((end-param.nVapor+1):end); % ug/m3 

    % remove vapor concentrations from y vector
    y = y(1:end-param.nVapor);

    % Find dilution rate, J, and GR from the matrices
    param.dilu = param.dilutionRateMatrix(2,find(t<=param.dilutionRateMatrix(1,:),1));
    param.J = param.JMatrix(2,find(t<=param.JMatrix(1,:),1));
    param.GR = param.GRMatrix(2,find(t<=param.GRMatrix(1,:),1));

    % if GR is modelled, J is also modelled
    if param.GRModel == 1
        param.J = nan(param.nVapor,1); % param.J becomes a vector
        
        for iVapor = 1:param.nVapor
            param.J(iVapor) = param.nucleationCoefficient(iVapor)*vaporConc(iVapor)^param.nucleationExponent(iVapor); % (#/cm3 s)
        end
    end
    
    % Plotting the distr.
    if param.plotDistrDuringSim
        figure(1)
        clf

        stairs(param.Dp_centers*1e9,y/(log10(param.Dp_centers(3))-log10(param.Dp_centers(2)))*param.particleNumberPenetration,'r')
        set(gca,'xscale','log','yscale','log')
        ylim([100 1e6])
        xlim([1 10000])
        title(strcat('t= ',num2str(t,'%1.3f'),' s'))
        xlabel('Dp (nm)')
        ylabel('dN/dlogDp (cm^{-3})')
        grid

        drawnow
    end
    
    % initialize vectors
    dy = zeros(size(y));
    C = zeros(size(y));
    MGR = nan(param.nVapor,1); % kg/s
    dVaporConc = zeros(param.nVapor,1); % ug/m3 s


    if param.losses == 1
        if param.GRModel == 1
            for iVapor = 1:param.nVapor
                % wall losses of vapors
                dVaporConc(iVapor) = dVaporConc(iVapor) - param.vaporWallLossRates(iVapor)*vaporConc(iVapor);
                % dilution of vapors
                dVaporConc(iVapor) = dVaporConc(iVapor) - vaporConc(iVapor)*param.dilu;
            end
        end
    end

    % transport equations for vapors
    if param.GRModel == 1
        dVaporConc = feval(param.vaporsTransportEquation,vaporConc,dVaporConc,param);
    end

    for i = 1:length(dy) % covers only particles, not vapors
        dp = param.Dp_centers(i);

        % dilution
        dy(i) = dy(i) - y(i)*param.dilu;

        % new particle formation (nucleation)
        if i==1 % particles are added only to the first size section (dCluster)
            if param.GRModel == 1 % nucleation reduces vapor concentrations if GR is modelled, and J is also modelled then
                for iVapor = 1:param.nVapor
                    dy(i) = dy(i) + param.J(iVapor); % #/cm3 s
                    dVaporConc(iVapor) = dVaporConc(iVapor) - param.J(iVapor)*pi/6*param.rho*dp^3*1e6*1e9; % ug/m3 s
                end
            else
                dy(i) = dy(i) + param.J; % #/cm3 s
            end
        end

        % condensation
        if param.GRModel == 1
            for iVapor = 1:param.nVapor
                MGR(iVapor) = ...
                    massGrowthRate(dp,param.dMolec(iVapor),param.T,param.visc,param.diffMolec(iVapor),vaporConc(iVapor)*1e-9,...
                    param.CSat(iVapor)*1e-9,param.rho,param.molarMass(iVapor),param.surfaceTension(iVapor),param.accommodationCoefficient(iVapor)); % kg/s

                C(i) = C(i) + MGR(iVapor)/(pi/2*param.rho*dp^2)/param.DeltaDp(i);

                dVaporConc(iVapor) = dVaporConc(iVapor) - MGR(iVapor)*y(i)*1e6*1e9; % ug/m3 s
            end
        else
            C(i) = param.GR/param.DeltaDp(i);
        end

        if i>1
            dy(i) = dy(i) - y(i)*C(i) + y(i-1)*C(i-1);
        else
            dy(i) = dy(i) - y(i)*C(i);
        end

        % depositional losses
        if param.losses == 1
            switch(param.lossesMethod)
                case 'exponent'
                    dy(i) = dy(i) - param.lossesCoeff*dp^param.lossesExponent*y(i);
                case 'vector'
                    dy(i) = dy(i) - param.lossesCoeff(i)*y(i);
                case 'equation'
                    dy(i) = dy(i) - eval(param.lossesCoeff)*y(i);
                otherwise
                    error('Losses method needs to be either ''exponent'', ''vector'', or ''equation''')
            end
        end

        % coagulation
        if param.coag
            cM = coagulationMatrix(param.Dp_centers,i);
            
            for j = 1:i
                NNkk = y(i)*y(j)*param.kk(i,j);
                
                if i == j
                    dy(j) = dy(j)-NNkk; % loss      
                    dy = dy + cM(j,:)'.*0.5.*NNkk; % gain
                else
                    dy(i) = dy(i)-NNkk; %loss
                    dy(j) = dy(j)-NNkk; %loss           
                    dy = dy + cM(j,:)'.*NNkk; % gain
                end
            end
        end

        % coagulational losses
        if param.coagSink
            dy(i) = dy(i) - param.coagSinkFactor*param.coagSinkN*y(i)*dp^param.coagSinkExponent*1e6;    
        end

    end

    % append vapor d vector to dy because vapor variables were removed from the y vector in the beginning
    dy = [dy;dVaporConc];

end
