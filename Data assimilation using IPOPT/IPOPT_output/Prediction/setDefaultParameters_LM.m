function [modelPars,varPars,sysPars,fun] = setDefaultParameters_LM()

    % Parameters for a lactotroph model
    modelPars.C  = 10.696926;
    modelPars.gCa  = 3.235192;
    modelPars.VCa   = 69.999904;
    modelPars.Vm   = -15.659561;
    modelPars.Sm =  14.311551;
    modelPars.gk   = 0.709192;
    modelPars.Vk   = -89.999739;
    modelPars.Vn   = -14.997238;
    modelPars.Sn   = 29.916643;
    modelPars.Taun   = 42.906489;
    modelPars.gSK = 0.998482;
    modelPars.Ks = 0.261894;
    modelPars.gBK = 1.323580;
    modelPars.Vf  = -17.595142;
    modelPars.Sf  = 9.999982;
    modelPars.TauBK  = 9.999984;
    modelPars.gL  = 0.409886;
    modelPars.VL   = -40.003664;
    modelPars.fc  = 0.012834;
    modelPars.Alpha   = 0.001500;
    modelPars.Kc   = 0.142569;

    varPars = [modelPars.C, modelPars.gCa, modelPars.VCa, modelPars.Vm, modelPars.Sm, modelPars.gk, ...
        modelPars.Vk, modelPars.Vn, modelPars.Sn, modelPars.Taun, modelPars.gSK, modelPars.Ks, ...
        modelPars.gBK, modelPars.Vf, modelPars.Sf, modelPars.TauBK, modelPars.gL, modelPars.VL, ...
        modelPars.fc, modelPars.Alpha, modelPars.Kc];
    
    sysPars.tstart = 0.0;
    sysPars.tfinal = 500.0;
    
    fun = @LactotrophModelODE;
    
end
