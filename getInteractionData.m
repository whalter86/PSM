function [parameters,variables,functions]=getInteractionData(parameters,variables,functions,interactions,genes)
% interactions(i).Identifier= 'i1';  % interaction identifier
% interactions(i).Source= 'g1';  % source gene identifier
% interactions(i).SourceType= 'Protein'; % Protein or mRNA
% interactions(i).Target= 'g2'; % target gene identifier
% interactions(i).Mode= 'tx'; % 'tx' or 'tl' 
% interactions(i).ParamNames={'k1'}; 
% interactions(i).ParamValues=[.000001,1];
% interactions(i).Fun='-k1 * u1';

%% get targets of interactions, incorporate new parameters and add interaction functions
intertargets={};
newinteractionparamnames={};
newinteractionparamvals=[];
for i = 1:length(interactions)
    % get targets of interactions
    intertargets=[intertargets,interactions(i).Target];
    
    % check interaction parameters
    if length(interactions(i).ParamNames) ~= length(interactions(i).ParamValues)
       error('Number of interaction parameter names and values does not match')
    end
    
    interactionparamnames=strcat([interactions(i).Identifier,'_'],interactions(i).ParamNames);
    % incorporate new parameters
    newinteractionparamnames=[newinteractionparamnames,interactionparamnames];
    newinteractionparamvals=[newinteractionparamvals,interactions(i).ParamValues];
    
    % add interaction functions
    funpnames=strjoin(interactions(i).ParamNames,',');
    funuwords=regexp(interactions(i).Fun,'\W?u\d*\W?','match');
    funuwords=unique( cellfun(@(x) regexp(x,'u\d*','match'),funuwords));
    fununames=strjoin(funuwords,',');
    functions.names=[functions.names;interactions(i).Identifier];
    functions.arguments=[functions.arguments;strjoin({fununames,'t',funpnames},',')];
    functions.vals=[functions.vals;interactions(i).Fun];
    
    % adapt variables
    if strcmp(interactions(i).Mode,'tx')
        targetName=[interactions(i).Target,'_init_transcr'];
    elseif strcmp(interactions(i).Mode,'tl')
        targetName=[interactions(i).Target,'_init_transl'];
    else
        error('wrong input');
    end
    if strcmp(interactions(i).SourceType,'Protein')
        sourceName=genes(strcmp({genes.ID},interactions(i).Source)).product;
    elseif strcmp(interactions(i).SourceType,'mRNA')
        sourceName= ['mRNA_',interactions(i).Source];
    else
        error('wrong input');
    end
    varIndex=strcmp(variables.names,targetName);
    stringToVary=variables.vals{varIndex};
    oldStr=regexp(stringToVary,'max\((.*?),0\)','tokens','once');
    oldStr=oldStr{1};
    newStr=[oldStr,' + ',interactions(i).Identifier,'(',strjoin({sourceName,'time',interactionparamnames},','),')'];
    variedString=replace(stringToVary,oldStr,newStr);
    variables.vals{varIndex}=variedString;
end
parameters.names=[parameters.names;newinteractionparamnames'];
parameters.vals=[parameters.vals;newinteractionparamvals'];

