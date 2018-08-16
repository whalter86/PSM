function [modelstates,statecoding]=createPSM(general,genes,interactions,inputs,MN)
%statecoding: first column stands for geneID, second column for type of
%state!
% type = 1: dna-slot
% type = 2: mRNA
% type = 3: rna-slot
% type = 4: product
% type = 5: metabolite
%% define headers


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


%% handle gene data

if length(general.Name)<5 || ~strcmp(general.Name(end-5:end),'.txtbc')
    dotind=strfind(general.Name,'.');
    if isempty(dotind)
        dotind=length(general.Name)+1;
    end
    filename=[general.Name(1:dotind-1),'.txtbc'];
else
    filename=general.Name;
end




intertargets={};
interactionparamnames={};
interactionparamvals=[];
for i = 1:length(interactions)
    intertargets=[intertargets,interactions(i).Target];
    interactionparamnames=[interactionparamnames,strcat([interactions(i).Identifyer,'_'],interactions(i).ParamNames)];
    interactionparamvals=[interactionparamvals,interactions(i).ParamValues];
end


intargets={};
inparamnames={};
inparamvals=[];
for i = 1:length(inputs)
    intargets=[intargets,inputs(i).Target];
    inparamnames=[inparamnames,strcat([inputs(i).Identifyer,'_'],inputs(i).ParamNames)];
    inparamvals=[inparamvals,inputs(i).ParamValues];
end


% modelstates={'RNAP','R','si','a'};
modelstates={'RNAP','R'};
geneprodtext={'d/dt(RNAP) = 0 \n';'d/dt(R) = 0 \n' };
% modelstates={};
% generalparamnames={'RNAP_width','R_width','transcr_speed','decay_RNA',...
%     'transl_speed','decay_prot','dilution','inRate','metaRate',...
%     'energymulti','Kmet','Ksynth'};
% generalparamnames={'RNAP_width','R_width','transcr_speed','decay_RNA',...
%     'transl_speed','decay_prot','dilution'};
generalparamnames= fieldnames(general);
generalparamnames(ismember(generalparamnames,'Name'))=[]; %remove Name parameter

if ~ismember('dilution',generalparamnames)
    general.dilution=0;
    generalparamnames=[generalparamnames;'dilution'];
end

statecoding=zeros(length(modelstates),2);
reactionindex=1;    



% geneparams={'numgenes','genelength','initrate_transcr','initrate_transl','effectivity'};
geneparams=fieldnames(genes);
geneparams(ismember(geneparams,'product'))=[]; %remove product parameter
geneproducts={genes.product};


geneparamnames={};
geneparamvals=[];
genevariables={};
genedecay=cell(1,length(geneproducts));
boundprot={};
boundprotrate={};

rnap_o_str='RNAP_o = ';
r_o_str='R_o = ';
for i = 1:length(genes)
    %%% modelstates
    geneslots(i)=round(genes(i).genelength/general.RNAP_width)+1;
    RNAslots(i)=round(genes(i).genelength/general.R_width)+1;
    rnap_o_str=[rnap_o_str,'g',num2str(i),'_numgenes *('];
    for j =1:geneslots(i)
        modelstates=[modelstates,['x',num2str(i),'_',num2str(j)]];
        statecoding=[statecoding;i,1];
        rnap_o_str=[rnap_o_str,modelstates{end},'+'];
    end
    
    modelstates=[modelstates,['mRNA_g',num2str(i)]];
    statecoding=[statecoding;i,2];
    r_o_str=[r_o_str,modelstates{end},' *('];
    for j =1:RNAslots(i)
        modelstates=[modelstates,['y',num2str(i),'_',num2str(j)]];
        statecoding=[statecoding;i,3];
        r_o_str=[r_o_str,modelstates{end},'+'];
    end 
    rnap_o_str(end)=')';
    r_o_str(end)=')';
    if i ~=length(genes)
        rnap_o_str=[rnap_o_str,'+'];
        r_o_str=[r_o_str,'+'];
    end

    if ~ismember(genes(i).product,modelstates)
        modelstates=[modelstates,genes(i).product];
        statecoding=[statecoding;i,4];
    else
        tmplog=logical(ismember(modelstates,genes(i).product));
        statecoding(tmplog,:)=[i,4];
        geneprodtext(find(ismember(modelstates,genes(i).product)))=[];
    end

    %%%gene parameters
    for j = 1:length(geneparams)
        newparamstring=['g',num2str(i),'_',geneparams{j}];
        newparamval=genes(i).(geneparams{j});
        geneparamnames=[geneparamnames,newparamstring];
        geneparamvals=[geneparamvals,newparamval];
    end
    
    %%%gene variables
    relevantInteractions=interactions(ismember(intertargets,genes(i).product));
    relevantInputs=inputs(ismember(intargets,genes(i).product));
    newvar1=['g',num2str(i),'_initrate_transcr'];
    newvar2=['g',num2str(i),'_initrate_transl'];
    for j =1:length(relevantInteractions)
        sourceindex=find(ismember(geneproducts,relevantInteractions(j).Source));
        if strcmp(relevantInteractions(j).SourceType,'Protein')
            sstring=genes(sourceindex).product;
        elseif strcmp(relevantInteractions(j).SourceType,'mRNA')
            sstring=['mRNA_g',num2str(sourceindex)];
        else 
            error('wrong input');
        end
%         tmpstr=['hillfun(',num2str(relevantInteractions(j).V),',',num2str(relevantInteractions(j).h),',',num2str(relevantInteractions(j).K),',',sstring,')'];
        pnames='';
        for k=1:length(relevantInteractions(j).ParamNames)
            pnames=[pnames,',',relevantInteractions(j).Identifyer,'_',relevantInteractions(j).ParamNames{k}];
        end
        tmpstr=[relevantInteractions(j).Identifyer,'(',sstring,',time',pnames,')'];
        if strcmp(relevantInteractions(j).Mode,'tx')
            newvar1=[newvar1,' + ',tmpstr];
        elseif strcmp(relevantInteractions(j).Mode,'tl')
            newvar2=[newvar2,' + ',tmpstr];
        elseif strcmp(relevantInteractions(j).Mode,'direct') 
            if strcmp(relevantInteractions(j).SourceType,'mRNA')
                targetindex=find(ismember(geneproducts,relevantInteractions(j).Target));
                tstring=['mRNA_g',num2str(targetindex)];
                reactiontext=[reactiontext;[sstring,' + ',tstring, ' =>  :r',num2str(reactionindex),' \n']];
                reactiontext=[reactiontext;['    vf= ',tmpstr,' * ',sstring,' * ',tstring, ' \n']];
                reactiontext=[reactiontext;' \n'];
                reactionindex=reactionindex+1;
            elseif strcmp(relevantInteractions(j).SourceType,'Protein')
                newprot = [sstring,'_x_',genes(i).product];
                modelstates=[modelstates,newprot];
                statecoding=[statecoding;0,4];
                boundprot=[boundprot,newprot];
                tmpstr=[relevantInteractions(j).Identifyer,'(',sstring,',',genes(i).product,',',newprot,',time',pnames,')'];
                
                
                genedecay{i}=[genedecay{i},' -',tmpstr,' + g',num2str(sourceindex),'_decay_prot * ',newprot];
                genedecay{sourceindex}=[genedecay{sourceindex},' -',tmpstr,' + g',num2str(i),'_decay_prot * ',newprot];
                boundprotrate=[boundprotrate,tmpstr];
                boundprotdilrate=['g',num2str(sourceindex),'_decay_prot + g',num2str(i),'_decay_prot'];
            end
%             genedecay{i}=[genedecay{i},' -',relevantInteractions(j).Fun,'*',genes(i).product,' *',sstring];
            
        else 
            error('wrong input');
        end
    end
    for j =1:length(relevantInputs)
        pnames='';
        for k=1:length(relevantInputs(j).ParamNames)
            pnames=[pnames,',',strcat([relevantInputs(j).Identifyer,'_'],relevantInputs(j).ParamNames{k})];
        end
        tmpstr=[relevantInputs(j).Identifyer,'(time',pnames,')'];
        if strcmp(relevantInputs(j).Mode,'tx')
            newvar1=[newvar1,' + ',tmpstr];
        elseif strcmp(relevantInputs(j).Mode,'tl')
            newvar2=[newvar2,' + ',tmpstr];
        else 
            error('wrong input');
        end
    end
    newvar1=['g',num2str(i),'_init_transcr = max(',newvar1,',0)'];
    newvar2=['g',num2str(i),'_init_transl = max(',newvar2,',0)'];
    newvar3=['c',genes(i).product,' = ( ',genes(i).product ,'/(Vcell*NA))*1000'];
    genevariables=[genevariables,newvar1,newvar2,newvar3];
end

% variables={'elongrate_transcr = transcr_speed / RNAP_width',...
%     'elongrate_transl = transl_speed / R_width',...
%     rnap_o_str,...
%     r_o_str,...
%     'RNAP_f = RNAP - RNAP_o',...
%     'R_f = R - R_o',...
%     'energy_eff = a/(Ksynth+a)'};
variables={'elongrate_transcr = transcr_speed / RNAP_width',...
    'elongrate_transl = transl_speed / R_width',...
    rnap_o_str,...
    r_o_str,...
    'RNAP_f = RNAP - RNAP_o',...
    'R_f = R - R_o'};
variables=[variables,genevariables];

%handle inputs
for i = 1:length(inputs)
    pnames='';
    for k=1:length(inputs(i).ParamNames)
        pnames=[pnames,',',inputs(i).ParamNames{k}];
    end
    functiontext = [functiontext; [inputs(i).Identifyer,'(t',pnames,')=',inputs(i).Fun,' \n']];
end
for i = 1:length(interactions)
    pnames='';
    for k=1:length(interactions(i).ParamNames)
        pnames=[pnames,',',interactions(i).ParamNames{k}];
    end
    uwords=regexp(interactions(i).Fun,'\W?u\d*\W?','match');
    uwords=unique( cellfun(@(x) regexp(x,'u\d*','match'),uwords));
    utext='';
    for j =1:length(uwords) 
        utext=[utext,uwords{j},','];
    end
    functiontext = [functiontext; [interactions(i).Identifyer,'(',utext,'t',pnames,')=',interactions(i).Fun,' \n']];
end

%handle metabolic network
if ~isempty(MN)
    MNparamnames=MN.ParamNames;
    MNparamvals=MN.ParamValues;
    % check sizes
    vnum=size(MN.S,2); %number of reactions
    % if vnum ~= length(MN.kcat) || vnum ~= length(MN.kM) || vnum ~= length(MN.Enzymes)
    %     error('Size of metabolic parameters does not match!');
    % end

    if size(MN.Names,1)>size(MN.Names,2)
        newstates=MN.Names';
    else
        newstates=MN.Names;
    end
    modelstates=[modelstates,newstates];
    % for i = 1:vnum
    %     MNparamnames=[MNparamnames,['kcat_v',num2str(i)],['kM_v',num2str(i)]];
    %     MNparamvals=[MNparamvals,MN.kcat(i),MN.kM(i)];
    % end
else
    MNparamnames={};
    MNparamvals=[];
    vnum=0;
end

%% write data to file



for i = 1:length(modelstates)
%     if ismember(modelstates{i},{'RNAP','R','si','a'})
    if ismember(modelstates{i},{'RNAP','R'})
        statetext=[statetext;[modelstates{i},'(0)=100',' \n']];
    else
        statetext=[statetext;[modelstates{i},'(0)=0',' \n']];
    end
end
for i = 1:length(generalparamnames)
    paramtext=[paramtext;[generalparamnames{i},'=',num2str(general.(generalparamnames{i})),' \n']];
end
for i = 1:length(geneparamnames)
    paramtext=[paramtext;[geneparamnames{i},'=',num2str(geneparamvals(i)),' \n']];
end
for i = 1:length(interactionparamnames)
    paramtext=[paramtext;[interactionparamnames{i},'=',num2str(interactionparamvals(i)),' \n']];
end
for i = 1:length(inparamnames)
    paramtext=[paramtext;[inparamnames{i},'=',num2str(inparamvals(i)),' \n']];
end
for i = 1:length(MNparamnames)
    paramtext=[paramtext;[MNparamnames{i},'=',num2str(MNparamvals(i)),' \n']];
end

for i = 1:length(variables)
    vartext=[vartext;[variables{i},' \n']];
end

% reactiontext=[reactiontext;['  => si :r',num2str(reactionindex),' \n']];
% reactiontext=[reactiontext;'    vf= inRate * P_BG \n']; %other protein?
% reactiontext=[reactiontext;' \n'];
% reactionindex=reactionindex+1;

% reactiontext=[reactiontext;[' si => :r',num2str(reactionindex),' \n']];
% reactiontext=[reactiontext;'    vf= metaRate*si/(Kmet+si) * P_BG \n']; %other protein?
% reactiontext=[reactiontext;' \n'];
% adottext=['d/dt(a) = -dilution * a + ',num2str(general.energymulti),'*r',num2str(reactionindex)];
% reactionindex=reactionindex+1;

% reactiontext=[reactiontext;[' si =>  :r',num2str(reactionindex),' \n']];
% reactiontext=[reactiontext;'    vf= dilution  * si \n'];
% reactiontext=[reactiontext;' \n'];
% reactionindex=reactionindex+1;

for i = 1:length(genes)
    %dilution and decay
    reactiontext=[reactiontext;[' mRNA_g',num2str(i),' =>  :r',num2str(reactionindex),' \n']];
%     reactiontext=[reactiontext;['    vf= (dilution + decay_RNA) * mRNA_g',num2str(i), ' \n']];
    reactiontext=[reactiontext;['    vf= (dilution + g',num2str(i),'_decay_RNA) * mRNA_g',num2str(i), ' \n']];
    reactiontext=[reactiontext;' \n'];
    reactionindex=reactionindex+1;
    
    %transcription reactions
    reactiontext=[reactiontext;[' => x',num2str(i),'_1 :r',num2str(reactionindex),' \n']];
%     reactiontext=[reactiontext;['    vf=energy_eff * g',num2str(i),'_initrate_transcr * RNAP_f * (1-x',num2str(i),'_1) \n']];
    reactiontext=[reactiontext;['    vf= g',num2str(i),'_init_transcr * RNAP_f * (1-x',num2str(i),'_1) \n']];
    reactiontext=[reactiontext;' \n'];
    reactionindex=reactionindex+1;
    for j =1:geneslots(i)-1
        reactiontext=[reactiontext;['x',num2str(i),'_',num2str(j),' => x',num2str(i),'_',num2str(j+1),' :r',num2str(reactionindex),' \n']];
%         reactiontext=[reactiontext;['    vf=energy_eff * elongrate_transcr * x',num2str(i),'_',num2str(j),' * (1-x',num2str(i),'_',num2str(j+1),') \n']];
        reactiontext=[reactiontext;['    vf=elongrate_transcr * x',num2str(i),'_',num2str(j),' * (1-x',num2str(i),'_',num2str(j+1),') \n']];
        reactiontext=[reactiontext;' \n'];
        reactionindex=reactionindex+1;
    end
    reactiontext=[reactiontext;['x',num2str(i),'_',num2str(geneslots(i)),' => ',num2str(genes(i).numgenes),'*mRNA_g',num2str(i),' :r',num2str(reactionindex),' \n']];
%     reactiontext=[reactiontext;['    vf= energy_eff * elongrate_transcr * x',num2str(i),'_',num2str(geneslots(i)),' \n']];
    reactiontext=[reactiontext;['    vf= elongrate_transcr * x',num2str(i),'_',num2str(geneslots(i)),' \n']];
    reactiontext=[reactiontext;' \n'];
    reactionindex=reactionindex+1;
    
    % translation reactions
    reactiontext=[reactiontext;[' => y',num2str(i),'_1 :r',num2str(reactionindex),' \n']];
%     reactiontext=[reactiontext;['    vf=energy_eff * g',num2str(i),'_initrate_transl * R_f * (1-y',num2str(i),'_1) \n']];
    reactiontext=[reactiontext;['    vf=g',num2str(i),'_init_transl * R_f * (1-y',num2str(i),'_1) \n']];
    reactiontext=[reactiontext;' \n'];
%     adottext=[adottext,' - mRNA_g',num2str(i),'*r',num2str(reactionindex)];
    reactionindex=reactionindex+1;
    for j =1:RNAslots(i)-1
        reactiontext=[reactiontext;['y',num2str(i),'_',num2str(j),' => y',num2str(i),'_',num2str(j+1),' :r',num2str(reactionindex),' \n']];
%         reactiontext=[reactiontext;['    vf=energy_eff * elongrate_transl * y',num2str(i),'_',num2str(j),' * (1-y',num2str(i),'_',num2str(j+1),') \n']];
        reactiontext=[reactiontext;['    vf=elongrate_transl * y',num2str(i),'_',num2str(j),' * (1-y',num2str(i),'_',num2str(j+1),') \n']];
        reactiontext=[reactiontext;' \n'];
%         adottext=[adottext,'- mRNA_g',num2str(i),'*r',num2str(reactionindex)];
        reactionindex=reactionindex+1;
    end 
    reactiontext=[reactiontext;['y',num2str(i),'_',num2str(RNAslots(i)),' =>  :r',num2str(reactionindex),' \n']];
%     reactiontext=[reactiontext;['    vf=energy_eff * elongrate_transl * y',num2str(i),'_',num2str(RNAslots(i)),' \n']];
    reactiontext=[reactiontext;['    vf=elongrate_transl * y',num2str(i),'_',num2str(RNAslots(i)),' \n']];
    reactiontext=[reactiontext;' \n'];
%     adottext=[adottext,'- mRNA_g',num2str(i),'*r',num2str(reactionindex)];
%     geneprodtext=[geneprodtext;['d/dt(',genes(i).product,') = -(dilution + decay_prot) *',genes(i).product,'+ g',num2str(i),'_effectivity * mRNA_g',num2str(i),'*r',num2str(reactionindex),' \n']];
    geneprodtext=[geneprodtext;['d/dt(',genes(i).product,') =', genedecay{i} ,' -(dilution + g',num2str(i),'_decay_prot) *',genes(i).product,'+ mRNA_g',num2str(i),'*r',num2str(reactionindex),' \n']];
    reactionindex=reactionindex+1;
    
end

for i = 1:length(boundprot)
    geneprodtext=[geneprodtext;['d/dt(',boundprot{i},') =', boundprotrate{i} ,' -(dilution +',boundprotdilrate ,') *',boundprot{i},' \n']];
end


for i = 1: vnum
    substrateindices=MN.S(:,i)<0;
    productindices=MN.S(:,i)>0;
    substrates=MN.Names(substrateindices');
    subsSfakt=-MN.S(substrateindices,i);
    products=MN.Names(productindices');
    prodSfakt=MN.S(productindices,i);
    if isempty(substrates)
        lhs=' ';
    else
        lhs=[num2str(subsSfakt(1)),'*',substrates{1}];
        for j =2:length(substrates)
            lhs=[lhs,' + ',num2str(subsSfakt(j)),'*',substrates{j}];
        end
    end
    if isempty(products)
        rhs=' ';
    else
        rhs=[num2str(prodSfakt(1)),'*',products{1}];
        for j =2:length(products)
            rhs=[rhs,' + ',num2str(prodSfakt(j)),'*',products{j}];
        end
    end

    
   
    reactiontext=[reactiontext;[lhs,' => ',rhs,' :r',num2str(reactionindex),' \n']];
%     if isempty(substrates)
%         reactiontext=[reactiontext;['    vf= (kcat_v',num2str(i),') * ',MN.Enzymes{i}, ' \n']];
%     elseif length(substrates)==1
%         reactiontext=[reactiontext;['    vf= (kcat_v',num2str(i),') * ',MN.Enzymes{i},' * ',substrates{1},' / ( kM_v',num2str(i),' + ',substrates{1} ,' ) \n']];
%     else
%         error('Only single substrates are supported for now');
%     end
    reactiontext=[reactiontext;['    vf= ',MN.Fun{i},' \n']];

    reactiontext=[reactiontext;' \n'];
    reactionindex=reactionindex+1;
end
% adottext=[adottext,' \n'];

% alltext=[headertext;adottext;geneprodtext;statetext;paramtext;vartext;reactiontext;footertext];
alltext=[headertext;geneprodtext;statetext;paramtext;vartext;reactiontext;functiontext;footertext];

fid = fopen(filename,'w');
for row = 1:size(alltext,1)
    fprintf(fid,alltext{row,:});
end
fclose(fid);