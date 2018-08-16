function [headertext,statetext,paramtext,vartext,reactiontext,functiontext,footertext] = initializeHeaders(general)
%%% define headers
headertext={'********** MODEL NAME \n';...
[general.Name,' \n'];...
'********** MODEL NOTES \n';...
'********** MODEL STATE INFORMATION \n'};

statetext={};
paramtext={'********** MODEL PARAMETERS \n'};
vartext={'********** MODEL VARIABLES \n'};
reactiontext={'********** MODEL REACTIONS \n'};

functiontext={
'********** MODEL FUNCTIONS \n';...
'hillfun(V,h,K,u) = V * (u^h)/(u^h+K^h) \n'
'unitstep(t) = max(sign(t),0) \n'};

footertext={'********** MODEL EVENTS \n';...
'********** MODEL MATLAB FUNCTIONS \n';...
'\n'};

end