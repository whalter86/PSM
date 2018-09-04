function [statenames,statecodes]=createPSM(general,genes,varargin)
%statecoding: first column stands for geneID, second column for type of
%state!
% statecoding = 1: dna-slot
% statecoding = 2: mRNA
% statecoding = 3: rna-slot
% statecoding = 4: product
% statecoding = 5: metabolite
% usage
% [modelstates,statecoding]=createPSM(general,genes)
% [modelstates,statecoding]=createPSM(general,genes,interactions)
% [modelstates,statecoding]=createPSM(general,genes,interactions,inputs)
% [modelstates,statecoding]=createPSM(general,genes,interactions,inputs,MN)

switch nargin
    case 1
        error('Not enough inputs, at least general and genetic information is required')
    case 2
        interactions=[];
        inputs=[];
        MN=[];
    case 3  
        interactions=varargin{1};
        inputs=[];
        MN=[];
    case 4
        interactions=varargin{1};
        inputs=varargin{2};
        MN=[]; 
    case 5
        interactions=varargin{1};
        inputs=varargin{2};
        MN=varargin{3};
    otherwise
        warning('Too many inputs, all inputs after number 5 will be ignored')
end

%% check for standard parameters
generalparamnames= fieldnames(general);         
if ~ismember('Name',generalparamnames)
    general.Name ='PSM';
end
if ~ismember('RNAP_width',generalparamnames)
    general.RNAP_width =40;
end
if ~ismember('R_width',generalparamnames)
    general.R_width =76;
end
if ~ismember('transcr_speed',generalparamnames)
    general.transcr_speed =3300;
end
if ~ismember('transl_speed',generalparamnames)
    general.transl_speed =2970;
end
if ~ismember('NA',generalparamnames)
    general.NA =6.0221e+23;
end
if ~ismember('Vcell',generalparamnames)
    general.Vcell =6.7e-16;
end
if ~ismember('dilution',generalparamnames)
    general.dilution =0;
end
if ~ismember('RNAP_IC',generalparamnames)
    general.RNAP_IC =4600;
end
if ~ismember('R_IC',generalparamnames)
    general.R_IC =39400;
end


%% get basic gene data
[geneproduction,states,parameters,variables,reactions,functions]=getGeneData(general,genes);
% todo: energy/resource dependenc

%% incorporate genetic interactions
[parameters,variables,functions]=getInteractionData(parameters,variables,functions,interactions,genes);

%% incorporate inputs
[parameters,variables,functions]=getInputData(parameters,variables,functions,inputs);

%% incorporate metabolic network
[states,parameters,reactions]=getMetabolicNetworkData(states,parameters,reactions,MN);

%% write data to file
writeToFile(general,geneproduction,states,parameters,variables,reactions,functions);

%% outputs
statenames=states.names;
statecodes=states.code;
