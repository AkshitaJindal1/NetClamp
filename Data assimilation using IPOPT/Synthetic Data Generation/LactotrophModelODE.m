% Lactotroph model with hard-coded parameters
% Akshita Jindal
% 11.10.2022

% NOTE: varPars(vp) 
% vp = [ C, gCa, VCa, Vm, Sm, gk, Vk, Vn, Sn, Taun, gSK, Ks, gBK, Vf, Sf, TauBK, gL, VL, fc, Alpha, Kc]
% I_app is used for stimulation and is a time dependent variable.


function dydt = LactotrophModelODE(t, y, vp, I_app)

%local/state variables - 4 - time dependent variables and can be plotted with t 
  V  = y(1,1);  %bounds -80, 20
  n  = y(2,1);   % 0 , 1
  f  = y(3,1);    % 0, 1
  Ca = y(4,1);     % 0, 1
 
% assigning numbers to the parameters - 21 parameters
  C =vp(1);
  gCa = vp(2);
  VCa =vp(3);
  Vm =vp(4);
  Sm =vp(5);
  gk =vp(6);
  Vk =vp(7);
  Vn =vp(8);
  Sn =vp(9);
  Taun =vp(10);
  gSK =vp(11);
  Ks =vp(12);
  gBK =vp(13);
  Vf =vp(14);
  Sf =vp(15);
  TauBK =vp(16);
  gL =vp(17);
  VL =vp(18);
  fc =vp(19);
  Alpha =vp(20);
  Kc =vp(21);

  minf = (1 + exp ((Vm - V) / Sm))^(-1);
  ICa  = gCa * minf * (V - VCa);
   
  IK =  gk * n * (V - Vk);
  ninf = (1 + exp((Vn - V) / Sn))^(-1);
  sinf = Ca^2 / (Ca^2 + Ks^2);
  ISK  = gSK * sinf * (V - Vk);
  
  finf= (1 + exp((Vf - V) / Sf))^(-1);
  IBK = gBK * f * (V - Vk);
  Ileak  = gL * (V - VL);

  dydt = [-(ICa + IK + ISK + IBK + Ileak - I_app(t)) / C;
         (ninf - n) / Taun;
         (finf - f) / TauBK;
         -fc * (Alpha * ICa + Kc * Ca)];

end
