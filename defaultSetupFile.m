% default input file

% Remove previous parameters except interpolation table variables
clear p2
p2.alphaTaulukko = p.alphaTaulukko;
p2.dTaulukko = p.dTaulukko;
p = p2;
clear p2

% How to represent particle size distribution?
% PL:    power law
% LN:    log-normal
% PLLN:  combined power law and log-normal
% FSn:   fixed-sectional, where n is the number of sections 
p.model = 'FS40';

% Method for solving a system of equations in PL and PLLN models
% 1: interpolation table
% 2: Levenberg-Marquardt iterative algorithm
p.PLEquations = 2;

% new particle formation rates, J
% doesn't have an effect if GR is modelled because J is also modelled then
% first row:    times (s)
% second row:   J (#/cm3 s)
p.JMatrix = [0 3600; 0 100];

% Diameter of a newly formed particle (m)
% It is also the lowest diameter in the FS model.
p.dCluster = 1e-9;

% The highest diameter in the FS model
p.highestDiameter = 1000e-9;

% Are growth rates modelled? Also J becomes modelled
p.GRModel = 0;
% If not, growth rates are given as constants...

    % condensational growth rates, GR
    % first row:    times (s)
    % second row:   GR (m/s)
    p.GRMatrix = [0 3600; 0 3e-12];

% If they are...

    % number of vapors
    p.nVapor = 2;

    % names of the vapors
    p.vaporName = {'IVOC';'SVOC'};

    % initial concentrations of the vapors (ug/m3) 
    p.vaporConc0 = [400; 50];

    % molar masses of the vapors (kg/mol)
    p.molarMass = [0.2; 0.3];

    % diameters of the vapor molecules (m)
    p.dMolec = [0.8e-9; 0.9e-9];

    % diffusion coefficients of the vapors in air (m2/s)
    p.diffMolec = [0.1e-4; 0.1e-4];

    % saturation concentrations of the vapors (ug/m3) 
    p.CSat = [1e5; 1e2];

    % surface tensions of the vapors (N/m) 
    p.surfaceTension = [0.03; 0.03];

    % accommodation coefficients of the vapors (dimensionless)
    p.accommodationCoefficient = [0.1; 1];

    % nucleation exponents for the vapors (dimensionless)
    p.nucleationExponent = [0; 2];

    % nucleation coefficients for the vapors ((ug/m3)^(-nucleation exponent)/cm3 s)
    p.nucleationCoefficient = [0; 0.01];

    % wall loss rates of the vapors (1/s)
    p.vaporWallLossRates = [3e-4; 3e-4];

    % name of the transport equation for the vapors
    p.vaporsTransportEquation = 'defaultVaporsTransportEquation';



% Is condensational transfer on?
% Only applicable with the PLLN model
p.condensationalTransfer = 1;

% Condensational transfer factor
p.condensationalTransferFactor = 0.5;

% Is coagulation on?
p.coag = 1;

% Is coagulational transfer from the PL distr. to LN distr. on?
% Only applicable with the PLLN model
p.coagulationalTransfer = 1;

% Size bins in numerical integration of coagulation terms
p.binsInCoagulation = 20;

% Are integrals of coagulation terms always calculated numerically?
p.numericCoagulation = 0;

% Temperature (K)
p.T = 300;

% Air viscosity (Pa s)
p.visc = viscosity(p.T);

% Particle density (kg/m3)
p.rho = 1000;

% Is coagulation to background mode modelled?
p.coagSink = 0;

% If it is...
    % CMD of background distribution (m)
    p.coagSinkCMD = 100e-9;
    % Exponent
    p.coagSinkExponent = -1.6;
    % Number conc of background distribution (#/cm3)
    p.coagSinkN = 1e3;

% Are particle losses modelled?
p.losses = 1;
% If they are...
    % Loss method ('exponent', 'vector', or 'equation')
    p.lossesMethod = 'exponent';
    % Loss coefficient (m/s)
    p.lossesCoeff = 5e-12;
    % Exponent
    p.lossesExponent = -1;

% Penetration of total particle number in the sampling lines
% If less than 1, it means that the initial particle concentrations are
% actually higher in the model than given. Also the outputted
% concentrations will be lower that what model internally calculates.
p.particleNumberPenetration = 1;

% dilution rates
% first row:    times (s)
% second row:   dilution rate (1/s), nonnegative number
p.dilutionRateMatrix = [0 3600; 0 0];
    
% Time vector (s)
p.timeVec = 0:3600;

% Method for giving initial PSD 
p.initialPSDMethod = 'parameters';

% Initial PSD parameters if parameters are given
% [NPL, alpha, D2, NLN, GSD, CMD]
p.initialPSDParameters = [0, 0, 0, 1e5, 2, 100e-9];

% Initial moment vector
p.initialMomentVec = zeros(1,6);

% ODE solver
p.solverName = 'ode45';

% Relative tolerance for ODE solver
p.relativeTolerance = 1e-3;

% Plot distribution during simulation
p.plotDistrDuringSim = 1;

% Plot waitbar during simulation
p.plotWaitbarDuringSim = 1;

% Plot output after simulation
p.plotOutputAfterSim = 1;
