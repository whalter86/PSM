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
% todo
            


%% get basic gene data
[geneproduction,states,parameters,variables,reactions,functions]=getGeneData(general,genes);
% todo: energy/resource dependenc

%% incorporate genetic interactions
[parameters,variables,functions]=getInteractionData(parameters,variables,functions,interactions,genes);

%% incorporate inputs
[parameters,variables,functions]=getInputData(parameters,variables,functions,inputs);

%% incorporate metabolic network
% todo

%% write data to file
writeToFile(general,geneproduction,states,parameters,variables,reactions,functions);

%% outputs
statenames=states.names;
statecodes=states.code;
