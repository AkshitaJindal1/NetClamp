function [modelpars, varpars, sysPars,fun] = setDefaultParameters_opsinmodel_ODE()

    fun = @opsinmodel_ODE;

    % Parameters for 4-state opsin model (Williams et al., 2013) 
    modelpars.Gd1   = 0.08384246453338509;   % 0.01 to 0.1S    %0.09
    modelpars.Gd2 = 0.27859819229447524; % 0.01 to 0.1S
    modelpars.Gr   = -0.0033678682368642673;% 1e-5 to 1e-3       %1e-5
    modelpars.e12   = 0.05164699251693067;     %Gb: 0.001 to 0.1 per second
    modelpars.e21   = 0.3573660258332245;     %Gf : 0.001 to 0.1 per second
    modelpars.epsilon1   = 3.2477173862362223;    %Ga1: 0.1 to 1
    modelpars.epsilon2   = 0.3209832217446309;       %Ga2: 0.1 to 1
    modelpars.I  = 0;                  % depends on the experimental setup
    modelpars.lambda  = 470;         %450 to 500 nm
    modelpars.sigma = 6e-4; 
    modelpars.w  = 0.77;         %stays the same
    
    varpars = [modelpars.Gd1; modelpars.Gd2; ...
      modelpars.Gr; modelpars.e12; modelpars.e21;...
      modelpars.epsilon1; modelpars.epsilon2; modelpars.I; ...
      modelpars.lambda; modelpars.sigma; modelpars.w]; 
  
    sysPars.tstart = 0;
    sysPars.tfinal = 100;

end
