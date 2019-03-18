function [parameters,variables,functions]=getInputData(parameters,variables,functions,inputs)
% inputs(i).Identifier= 'u1';  % Input identifier
% inputs(i).Target= 'P1'; 
% inputs(i).Mode= 'tx'; % 'tx' or 'tl'
% inputs(i).ParamNames={'P1','P2'};
% inputs(i).ParamValues=[1e-3,0];
% inputs(i).Fun='P1*cos(P2*t)';


%% get targets of inputs, incorporate new parameters and add input functions
intargets={};
newinputparamnames={};
newinputparamvals=[];
for i = 1:length(inputs)
    % get targets of inputs
    intargets=[intargets,inputs(i).Target];
    
    % check interaction parameters
    if length(inputs(i).ParamNames) ~= length(inputs(i).ParamValues)
       error('Number of input parameter names and values does not match')
    end
    inputparamnames=strcat([inputs(i).Identifier,'_'],inputs(i).ParamNames);
    % incorporate new parameters
    newinputparamnames=[newinputparamnames,inputparamnames];
    newinputparamvals=[newinputparamvals,inputs(i).ParamValues];
    
    % add input functions
    funpnames=strjoin(inputs(i).ParamNames,',');
    functions.names=[functions.names;inputs(i).Identifier];
    functions.arguments=[functions.arguments;strjoin({'t',funpnames},',')];
    functions.vals=[functions.vals;inputs(i).Fun];
    
    % adapt variables
    if strcmp(inputs(i).Mode,'tx')
        targetName=[inputs(i).Target,'_init_transcr'];
    elseif strcmp(inputs(i).Mode,'tl')
        targetName=[inputs(i).Target,'_init_transl'];
    else
        error('wrong input');
    end
    
    varIndex=strcmp(variables.names,targetName);
    stringToVary=variables.vals{varIndex};
    oldStr=regexp(stringToVary,'max\((.*),0\)','tokens','once');
    oldStr=oldStr{1};
    multsplits=regexp(oldStr,'*','split');
    oldaddStr=regexp(multsplits{end},'\((.*)\)','tokens','once');
    if isempty(oldaddStr) % meaning no brackets are found, i.e. only constant init rate
        oldaddStr=multsplits{end};
    else
        oldaddStr=oldaddStr{1};
    end    
    newaddStr=[oldaddStr,' + ',inputs(i).Identifier,'(',strjoin({'time',inputparamnames},','),')'];
    variedString=replace(stringToVary,oldaddStr,newaddStr);
    variables.vals{varIndex}=variedString;
end
parameters.names=[parameters.names;newinputparamnames'];
parameters.vals=[parameters.vals;newinputparamvals'];