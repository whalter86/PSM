function [geneproduction,states,parameters,variables,reactions,functions]=getGeneData(general,genes)
%statecoding: first column stands for geneID, second column for type of
%state!
% statecoding = 1: dna-slot
% statecoding = 2: mRNA
% statecoding = 3: rna-slot
% statecoding = 4: product
% statecoding = 5: metabolite
 
%% initialize empty output variables
geneproduction.names={};
geneproduction.rhs={};
states.names={};
states.code=[];
states.IC=[];
parameters.names={};
parameters.vals=[];
variables.names={};
variables.vals={};
reactions.names={};
reactions.products={};
reactions.operators={};
reactions.educts={};
reactions.forwardrate={};
reactions.backwardrate={};
functions.names={};
functions.arguments={};
functions.vals={};



%% get general parameters
generalparamnames= fieldnames(general);
% remove Name and IC parameters
generalparamnames(ismember(generalparamnames,{'Name','RNAP_IC','R_IC','AGTP_IC','CUTP_IC','AA_IC'}))=[]; 
% set default dilution parameter if not exists
if ~ismember('dilution',generalparamnames)
    general.dilution=0;
    generalparamnames=[generalparamnames;'dilution'];
end
for i = 1:length(generalparamnames)
    generalparamvals(i) = getfield(general,generalparamnames{i});
end
parameters.names=[parameters.names;generalparamnames];
parameters.vals=[parameters.vals;generalparamvals'];

%% include RNAP and R as proteins if not exist
geneproducts={genes.product};
if ~ismember('RNAP',geneproducts)
    geneproduction.names=[geneproduction.names;'RNAP'];
    geneproduction.rhs=[geneproduction.rhs;'0'];
    states.names=[states.names;'RNAP'];
    states.code=[states.code;[0 4]];
end
if ~ismember('R',geneproducts)
    geneproduction.names=[geneproduction.names;'R'];
    geneproduction.rhs=[geneproduction.rhs;'0'];
    states.names=[states.names;'R'];
    states.code=[states.code;[0 4]];
end

%% include energy resources
states.names=[states.names;'AGTP'];
states.code=[states.code;[0 4]];
states.names=[states.names;'CUTP'];
states.code=[states.code;[0 4]];
states.names=[states.names;'AA'];
states.code=[states.code;[0 4]];
% if ~general.Consumption_flag
%     geneproduction.names=[geneproduction.names;'AGTP'];
%     geneproduction.rhs=[geneproduction.rhs;'0'];
%     geneproduction.names=[geneproduction.names;'CUTP'];
%     geneproduction.rhs=[geneproduction.rhs;'0'];
%     geneproduction.names=[geneproduction.names;'AA'];
%     geneproduction.rhs=[geneproduction.rhs;'0'];
% end


%% get gene parameternames
geneparams=fieldnames(genes);
geneparams(ismember(geneparams,{'product','ID'}))=[]; %remove product parameter

rnap_o_str='';
r_o_str='';

reactionindex=1;
for i = 1:length(genes)
    %%% states
    geneslots(i)=round(genes(i).genelength/general.RNAP_width)+1;
    RNAslots(i)=round(genes(i).genelength/general.R_width)+1;
    rnap_o_str=[rnap_o_str,genes(i).ID,'_numgenes *('];
    for j =1:geneslots(i)
        states.names=[states.names;['x',num2str(i),'_',num2str(j)]];
        states.code=[states.code;i,1];
        rnap_o_str=[rnap_o_str,states.names{end},'+'];
    end
    rnap_o_str(end)=')';
    states.names=[states.names;['mRNA_',genes(i).ID]];
    states.code=[states.code;i,2];
    r_o_str=[r_o_str,states.names{end},' *('];
    for j =1:RNAslots(i)
        states.names=[states.names;['y',num2str(i),'_',num2str(j)]];
        states.code=[states.code;i,3];
        r_o_str=[r_o_str,states.names{end},'+'];
    end 
    r_o_str(end)=')';
    if i ~=length(genes)
        rnap_o_str=[rnap_o_str,'+'];
        r_o_str=[r_o_str,'+'];
    end
    if ~ismember(genes(i).product,states.names)
        states.names=[states.names;genes(i).product];
        states.code=[states.code;i,4];
    else
        tmplog=logical(ismember(states.names,genes(i).product));
        states.code(tmplog,:)=[i,4];
%         geneprodtext(find(ismember(modelstates,genes(i).product)))=[];
    end

    %%% gene parameters
    for j = 1:length(geneparams)
        newparamstring=[genes(i).ID,'_',geneparams{j}];
        newparamval=genes(i).(geneparams{j});
        parameters.names=[parameters.names;newparamstring];
        parameters.vals=[parameters.vals;newparamval];
    end

    %%% gene variables
    variables.names=[variables.names;[genes(i).ID,'_init_transcr']];       
    variables.names=[variables.names;[genes(i).ID,'_init_transl']];
    variables.names=[variables.names;['c',genes(i).product]];
    variables.vals=[variables.vals;['max(',genes(i).ID,'_initrate_transcr',',0)']];
    variables.vals=[variables.vals;['max(',genes(i).ID,'_initrate_transl',',0)']];
    variables.vals=[variables.vals;['( ',genes(i).product ,'/(Vcell*NA))*1000']];
    
    
    %%% set reactions
    %%% mRNA degradation
    reactions.names=[reactions.names;['r',num2str(reactionindex)]];
    reactions.educts=[reactions.educts;[' mRNA_',genes(i).ID]];
    reactions.operators=[reactions.operators;'=>'];
    reactions.products=[reactions.products;' '];
    reactions.forwardrate=[reactions.forwardrate;['(dilution + ',genes(i).ID,'_decay_RNA) * mRNA_',genes(i).ID]];
    reactions.backwardrate=[reactions.backwardrate;' '];
    reactionindex=reactionindex+1;
    
    %%% transcription initiation
    reactions.names=[reactions.names;['r',num2str(reactionindex)]];
    reactions.educts=[reactions.educts;' '];
    reactions.operators=[reactions.operators;'=>'];
    reactions.products=[reactions.products;['x',num2str(i),'_1']];
    reactions.forwardrate=[reactions.forwardrate;[genes(i).ID,'_init_transcr * RNAP_f * (1-x',num2str(i),'_1)']];
    reactions.backwardrate=[reactions.backwardrate;' '];
    reactionindex=reactionindex+1;
    
    
    %%% transcription elongation 
    for j =1:geneslots(i)-1
        reactions.names=[reactions.names;['r',num2str(reactionindex)]];
        reactions.educts=[reactions.educts;['x',num2str(i),'_',num2str(j)]];
        reactions.operators=[reactions.operators;'=>'];
        reactions.products=[reactions.products;['x',num2str(i),'_',num2str(j+1)]];
        reactions.forwardrate=[reactions.forwardrate;['elongrate_transcr * x',num2str(i),'_',num2str(j),' * (1-x',num2str(i),'_',num2str(j+1),')']];
        reactions.backwardrate=[reactions.backwardrate;' '];
        reactionindex=reactionindex+1;
        %%% (including energy consumption)
        reactions.names=[reactions.names;['r',num2str(reactionindex)]];
        reactions.educts=[reactions.educts;['AGTP + CUTP']];
        reactions.operators=[reactions.operators;'=>'];
        reactions.products=[reactions.products;' '];
        reactions.forwardrate=[reactions.forwardrate;[num2str(genes(i).numgenes),'*Consumption_flag*1/2*RNAP_width*elongrate_transcr * x',num2str(i),'_',num2str(j),' * (1-x',num2str(i),'_',num2str(j+1),')']];
        reactions.backwardrate=[reactions.backwardrate;' '];
        reactionindex=reactionindex+1;
    end
    %%% transcription termination 
    reactions.names=[reactions.names;['r',num2str(reactionindex)]];
    reactions.educts=[reactions.educts;['x',num2str(i),'_',num2str(geneslots(i))]];
    reactions.operators=[reactions.operators;'=>'];
    reactions.products=[reactions.products;[num2str(genes(i).numgenes),'*mRNA_',genes(i).ID]];
    reactions.forwardrate=[reactions.forwardrate;['elongrate_transcr * x',num2str(i),'_',num2str(geneslots(i))]];
    reactions.backwardrate=[reactions.backwardrate;' '];
    reactionindex=reactionindex+1;
    %%% (including energy consumption)
    reactions.names=[reactions.names;['r',num2str(reactionindex)]];
    reactions.educts=[reactions.educts;'AGTP + CUTP'];
    reactions.operators=[reactions.operators;'=>'];
    reactions.products=[reactions.products;' '];
    reactions.forwardrate=[reactions.forwardrate;[num2str(genes(i).numgenes),'*Consumption_flag*1/2*RNAP_width*elongrate_transcr * x',num2str(i),'_',num2str(geneslots(i))]];
    reactions.backwardrate=[reactions.backwardrate;' '];
    reactionindex=reactionindex+1;
    
    %%% translation initiation
    reactions.names=[reactions.names;['r',num2str(reactionindex)]];
    reactions.educts=[reactions.educts;' '];
    reactions.operators=[reactions.operators;'=>'];
    reactions.products=[reactions.products;['y',num2str(i),'_1']];
    reactions.forwardrate=[reactions.forwardrate;[genes(i).ID,'_init_transl * R_f * (1-y',num2str(i),'_1)']];
    reactions.backwardrate=[reactions.backwardrate;' '];
    reactionindex=reactionindex+1;
    
    %%% translation elongation
    for j =1:RNAslots(i)-1
        reactions.names=[reactions.names;['r',num2str(reactionindex)]];
        reactions.educts=[reactions.educts;['y',num2str(i),'_',num2str(j)]];
        reactions.operators=[reactions.operators;'=>'];
        reactions.products=[reactions.products;['y',num2str(i),'_',num2str(j+1)]];
        reactions.forwardrate=[reactions.forwardrate;['elongrate_transl * y',num2str(i),'_',num2str(j),' * (1-y',num2str(i),'_',num2str(j+1),')']];
        reactions.backwardrate=[reactions.backwardrate;' '];
        reactionindex=reactionindex+1;
    %%% (including energy consumption)
        reactions.names=[reactions.names;['r',num2str(reactionindex)]];
        reactions.educts=[reactions.educts;'AGTP + AA'];
        reactions.operators=[reactions.operators;'=>'];
        reactions.products=[reactions.products;' '];
        reactions.forwardrate=[reactions.forwardrate;['mRNA_',genes(i).ID,'*Consumption_flag*1/3*R_width*elongrate_transl * y',num2str(i),'_',num2str(j),' * (1-y',num2str(i),'_',num2str(j+1),')']];
        reactions.backwardrate=[reactions.backwardrate;' '];
        reactionindex=reactionindex+1;
    end 
    
    %%% translation termination
    reactions.names=[reactions.names;['r',num2str(reactionindex)]];
    reactions.educts=[reactions.educts;['y',num2str(i),'_',num2str(RNAslots(i))]];
    reactions.operators=[reactions.operators;'=>'];
    reactions.products=[reactions.products;' '];
    reactions.forwardrate=[reactions.forwardrate;['elongrate_transl * y',num2str(i),'_',num2str(RNAslots(i))]];
    reactions.backwardrate=[reactions.backwardrate;' '];
    reactionindex=reactionindex+1;
    %%% (including energy consumption)
    reactions.names=[reactions.names;['r',num2str(reactionindex)]];
    reactions.educts=[reactions.educts;'AGTP + AA'];
    reactions.operators=[reactions.operators;'=>'];
    reactions.products=[reactions.products;' '];
    reactions.forwardrate=[reactions.forwardrate;['mRNA_',genes(i).ID,'*Consumption_flag*1/3*R_width*elongrate_transl * y',num2str(i),'_',num2str(RNAslots(i))]];
    reactions.backwardrate=[reactions.backwardrate;' '];
    reactionindex=reactionindex+1;
    
%     geneproduction.names=[geneproduction.names;genes(i).product];
%     geneproduction.rhs=[geneproduction.rhs;[' -(dilution + ',genes(i).ID,'_decay_prot) *',genes(i).product,'+ mRNA_',genes(i).ID,'*r',num2str(reactionindex)]];
    
    
    %%% protein production
    reactions.names=[reactions.names;['r',num2str(reactionindex)]];
    reactions.educts=[reactions.educts;' '];
    reactions.operators=[reactions.operators;'=>'];
    reactions.products=[reactions.products;genes(i).product];
    reactions.forwardrate=[reactions.forwardrate;['mRNA_',genes(i).ID,'*elongrate_transl * y',num2str(i),'_',num2str(RNAslots(i))]];
    reactions.backwardrate=[reactions.backwardrate;' '];
    reactionindex=reactionindex+1;
    
    %%% degradation and dilution
    reactions.names=[reactions.names;['r',num2str(reactionindex)]];
    reactions.educts=[reactions.educts;genes(i).product];
    reactions.operators=[reactions.operators;'=>'];
    reactions.products=[reactions.products;' '];
    reactions.forwardrate=[reactions.forwardrate;['(dilution + ',genes(i).ID,'_decay_prot)*',genes(i).product]];
    reactions.backwardrate=[reactions.backwardrate;' '];
    reactionindex=reactionindex+1;
end

%% set general variables
variables.names=[variables.names;'elongrate_transcr'];
variables.vals=[variables.vals; 'hillfun(1,1,AGTP_half,AGTP)*hillfun(1,1,CUTP_half,CUTP) * transcr_speed / RNAP_width'];
    
variables.names=[variables.names;'elongrate_transl'];
variables.vals=[variables.vals; 'hillfun(1,1,AGTP_half,AGTP)*hillfun(1,1,AA_half,AA) * transl_speed / R_width'];

variables.names=[variables.names;'RNAP_o'];
variables.vals=[variables.vals; rnap_o_str ];

variables.names=[variables.names;'R_o'];
variables.vals=[variables.vals; r_o_str];

variables.names=[variables.names;'RNAP_f'];
variables.vals=[variables.vals; 'RNAP - RNAP_o'];

variables.names=[variables.names;'R_f'];
variables.vals=[variables.vals; 'R - R_o'];
    
%% set IC to zero except for RNAP and R_IC
states.IC=zeros(length(states.names),1);
states.IC(ismember(states.names,'RNAP'))=general.RNAP_IC;
states.IC(ismember(states.names,'R'))=general.R_IC;
states.IC(ismember(states.names,'AGTP'))=general.AGTP_IC;
states.IC(ismember(states.names,'CUTP'))=general.CUTP_IC;
states.IC(ismember(states.names,'AA'))=general.AA_IC;


