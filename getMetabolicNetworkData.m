function [states,parameters,reactions]=getMetabolicNetworkData(states,parameters,reactions,MN)
% example: A -> B
% MN(i).S           = [-1;1];  % Stoichometric matrix (# metabolites x # reactions)
% MN(i).Names       = {'A','B'};  % Name of metabolites 
% MN(i).IC          = [ 2,1]; % Initial condition of metabolites
% MN(i).ParamNames  = {'kcat','KM'};  % Parameter Names 
% MN(i).ParamValues = [1e-2,0.0017];  % Parameter values 
% MN(i).Fun         = {'kcat*cEnzyme*A / (KM+A)'};  % reaction function

newMNparamnames={};
newMNparamvals=[];
newMNstateNames={};
newMNstateIC=[];


for m=1:length(MN)
    if size(MN(m).ParamNames,1)>size(MN(m).ParamNames,2)
        newMNparamnames=[newMNparamnames;MN(m).ParamNames];
    else
        newMNparamnames=[newMNparamnames;MN(m).ParamNames'];
    end

    if size(MN(m).ParamValues,1)>size(MN(m).ParamValues,2)
        newMNparamvals=[newMNparamvals;MN(m).ParamValues];
    else
        newMNparamvals=[newMNparamvals;MN(m).ParamValues'];
    end    

    if size(MN(m).Names,1)>size(MN(m).Names,2)
        newMNstateNames=[newMNstateNames;MN(m).Names];
    else
        newMNstateNames=[newMNstateNames;MN(m).Names'];
    end

    if size(MN(m).IC,1)>size(MN(m).IC,2)
        newMNstateIC=[newMNstateIC;MN(m).IC];
    else
        newMNstateIC=[newMNstateIC;MN(m).IC'];
    end
end

%% check for duplicates

[newMNparamnames,ind_unique_param]=unique(newMNparamnames,'stable');
if length(newMNparamnames) ~= length(newMNparamvals) 
    warning('Duplicates found in parameter names. Duplicate entries will be neglected!');
    newMNparamvals=newMNparamvals(ind_unique_param);
end

[newMNstateNames,ind_unique_states]=unique(newMNstateNames,'stable');
if length(newMNstateNames) ~= length(newMNstateIC) 
    warning('Duplicates found in state names. Duplicate entries will be neglected!');
    newMNstateIC=newMNstateIC(ind_unique_states);
end

ind_exist_states=ismember(newMNstateNames,states.names);
if sum(ind_exist_states)~=0
    warning('Some states already exist, initial conditions will be neglected')
    newMNstateNames(ind_exist_states)=[];
    newMNstateIC(ind_exist_states)=[];
end

ind_exist_param=ismember(newMNparamnames,parameters.names);
if sum(ind_exist_param)~=0
    warning('Duplicates found in parameter names. Duplicate entries will be neglected!');
    newMNparamnames(ind_exist_param)=[];
    newMNparamvals(ind_exist_param)=[];
end


states.names=[states.names;newMNstateNames];
states.IC=[states.IC;newMNstateIC];
states.code=[states.code;repmat([0 5],length(newMNstateNames),1)];
parameters.names=[parameters.names;newMNparamnames];
parameters.vals=[parameters.vals;newMNparamvals];


MNreactionindex=1;
for m=1:length(MN)
    numofreactions=size(MN(m).S,2);
    for i = 1: numofreactions
        substrateindices=MN(m).S(:,i)<0;
        productindices=MN(m).S(:,i)>0;
        substrates=MN(m).Names(substrateindices');
        subsSfakt=-MN(m).S(substrateindices,i);
        products=MN(m).Names(productindices');
        prodSfakt=MN(m).S(productindices,i);

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
        reactions.forwardrate=[reactions.forwardrate;MN(m).Fun{i}];
        reactions.backwardrate=[reactions.backwardrate;' '];
        MNreactionindex=MNreactionindex+1;
    end
end
