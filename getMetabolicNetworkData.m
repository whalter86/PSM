function [states,parameters,reactions]=getMetabolicNetworkData(states,parameters,reactions,MN)
% example: A -> B
% MN(i).S           = [-1;1];  % Stoichometric matrix (# metabolites x # reactions)
% MN(i).Names       = {'A','B'};  % Name of metabolites 
% MN(i).IC          = [ 2,1]; % Initial condition of metabolites
% MN(i).ParamNames  = {'kcat','KM'};  % Parameter Names 
% MN(i).ParamValues = [1e-2,0.0017];  % Parameter values 
% MN(i).Fun         = {'kcat*cEnzyme*A / (KM+A)'};  % reaction function


if ~isempty(MN)
    if size(MN.ParamNames,1)>size(MN.ParamNames,2)
        newMNparamnames=MN.ParamNames;
    else
        newMNparamnames=MN.ParamNames';
    end
    
    if size(MN.ParamValues,1)>size(MN.ParamValues,2)
        newMNparamvals=MN.ParamValues;
    else
        newMNparamvals=MN.ParamValues';
    end    
        
    if size(MN.Names,1)>size(MN.Names,2)
        newMNstateNames=MN.Names;
    else
        newMNstateNames=MN.Names';
    end
    
    if size(MN.IC,1)>size(MN.IC,2)
        newMNstateIC=MN.IC;
    else
        newMNstateIC=MN.IC';
    end    
else
    newMNparamnames={};
    newMNparamvals=[];
    newMNstateNames={};
    newMNstateIC=[];
end


states.names=[states.names;newMNstateNames];
states.IC=[states.IC;newMNstateIC];
states.code=[states.code;repmat([0 5],length(newMNstateNames),1)];
parameters.names=[parameters.names;newMNparamnames];
parameters.vals=[parameters.vals;newMNparamvals];


numofreactions=size(MN.S,2);
MNreactionindex=1;
for i = 1: numofreactions
    substrateindices=MN.S(:,i)<0;
    productindices=MN.S(:,i)>0;
    substrates=MN.Names(substrateindices');
    subsSfakt=-MN.S(substrateindices,i);
    products=MN.Names(productindices');
    prodSfakt=MN.S(productindices,i);
    
    if isempty(substrates)
        substratestring=' ';
    else
        substratestring=[num2str(subsSfakt(1)),'*',substrates{1}];
        for j =2:length(substrates)
            substratestring=[substratestring,' + ',num2str(subsSfakt(j)),'*',substrates{j}];
        end
    end
    
    if isempty(products)
        productstring=' ';
    else
        productstring=[num2str(prodSfakt(1)),'*',products{1}];
        for j =2:length(products)
            productstring=[productstring,' + ',num2str(prodSfakt(j)),'*',products{j}];
        end
    end
    
    reactions.names=[reactions.names;['r_MN',num2str(MNreactionindex)]];
    reactions.educts=[reactions.educts;substratestring];
    reactions.operators=[reactions.operators;'=>'];
    reactions.products=[reactions.products;productstring];
    reactions.forwardrate=[reactions.forwardrate;MN.Fun{i}];
    reactions.backwardrate=[reactions.backwardrate;' '];
    MNreactionindex=MNreactionindex+1;
end
% 
% 
% %% get targets of inputs, incorporate new parameters and add input functions
% intargets={};
% newinputparamnames={};
% newinputparamvals=[];
% for i = 1:length(inputs)
%     % get targets of inputs
%     intargets=[intargets,inputs(i).Target];
%     
%     % check interaction parameters
%     if length(inputs(i).ParamNames) ~= length(inputs(i).ParamValues)
%        error('Number of input parameter names and values does not match')
%     end
%     inputparamnames=strcat([inputs(i).Identifier,'_'],inputs(i).ParamNames);
%     % incorporate new parameters
%     newinputparamnames=[newinputparamnames,inputparamnames];
%     newinputparamvals=[newinputparamvals,inputs(i).ParamValues];
%     
%     % add input functions
%     funpnames=strjoin(inputs(i).ParamNames,',');
%     functions.names=[functions.names;inputs(i).Identifier];
%     functions.arguments=[functions.arguments;strjoin({'t',funpnames},',')];
%     functions.vals=[functions.vals;inputs(i).Fun];
%     
%     % adapt variables
%     if strcmp(inputs(i).Mode,'tx')
%         targetName=[inputs(i).Target,'_init_transcr'];
%     elseif strcmp(inputs(i).Mode,'tl')
%         targetName=[inputs(i).Target,'_init_transl'];
%     else
%         error('wrong input');
%     end
%     
%     varIndex=strcmp(variables.names,targetName);
%     stringToVary=variables.vals{varIndex};
%     oldStr=regexp(stringToVary,'max\((.*?),0\)','tokens','once');
%     oldStr=oldStr{1};
%     newStr=[oldStr,' + ',inputs(i).Identifier,'(',strjoin({'time',inputparamnames},','),')'];
%     variedString=replace(stringToVary,oldStr,newStr);
%     variables.vals{varIndex}=variedString;
% end
% parameters.names=[parameters.names;newinputparamnames'];
% parameters.vals=[parameters.vals;newinputparamvals'];