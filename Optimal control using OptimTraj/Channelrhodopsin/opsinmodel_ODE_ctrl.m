function dydt = opsinmodel_ODE_ctrl(~, y, vp, I)

% t

% Extract variables 
O1 = y(1,:);
C2 = y(2,:);
O2 = y(3,:);
C1 = 1-(O1+O2+C2);

%asigning numbers to the parameters - 

Gd1 = vp(1);       %decay rate of O1
Gd2 = vp(2);       %decay rate of O2
Gr = vp(3);        %recovery rate of C2
e12 = vp(4);       %transition rate from O1 to O2
e21 = vp(5);        %transition rate from O2 to O21
epsilon1 = vp(6); %sensitivity factor 1
epsilon2 = vp(7);  %sensitivity factor 2
lambda = vp(9);    %Wavelength
sigma = vp(10); %absorption cross-section for retinal
w = vp(11); %Scaling factr for losses of phtons during to scattering/absorptions * Planck's constant times speed of light in vacuum 

% Calculate ChR2 current (I_opsin)
Phi = sigma * I *lambda/w; % Calculate photon flux
k1 = epsilon1 * Phi; % Calculate rate constant k1
k2 = epsilon2 * Phi; % Calculate rate constant k2

% define ODEs
dydt = [(k1 .* C1) + (e21 * O2) - O1 * (e12 + Gd1);  %dynamics of O1
    (Gd2 * O2) - C2 .* (Gr + k2);   %dynamics of C2
    (k2 .* C2) + (e12 * O1) - O2 * (e21 + Gd2)]; %dynamics of O2

end
