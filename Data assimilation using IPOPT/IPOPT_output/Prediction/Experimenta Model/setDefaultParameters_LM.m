function [modelPars,varPars,sysPars,fun] = setDefaultParameters_LM()

    % Parameters for a lactotroph model
    modelPars.C  = 11;
    modelPars.gCa  = 4;
    modelPars.VCa   = 50;
    modelPars.Vm   = -20.844206;
    modelPars.Sm = 8.688076;
    modelPars.gk   = 0;
    modelPars.Vk   = -82.505433;
    modelPars.Vn   = -14.997337;
    modelPars.Sn   = 29.996314;
    modelPars.Taun   = 49.851608;
    modelPars.gSK = 2;
    modelPars.Ks = 1;
    modelPars.gBK = 2;
    modelPars.Vf  = -35;
    modelPars.Sf  = 4.227308;
    modelPars.TauBK  = 2;
    modelPars.gL  = 1;
    modelPars.VL   = -40;
    modelPars.fc  = 0.020634;
    modelPars.Alpha   = 0.0015;
    modelPars.Kc   = 0.16;
    
    varPars = [modelPars.C, modelPars.gCa, modelPars.VCa, modelPars.Vm, modelPars.Sm, modelPars.gk, ...
        modelPars.Vk, modelPars.Vn, modelPars.Sn, modelPars.Taun, modelPars.gSK, modelPars.Ks, ...
        modelPars.gBK, modelPars.Vf, modelPars.Sf, modelPars.TauBK, modelPars.gL, modelPars.VL, ...
        modelPars.fc, modelPars.Alpha, modelPars.Kc];
    
    sysPars.tstart = 0.0;
    sysPars.tfinal = 500.0;
    
    fun = @LactotrophModelODE;
    
end