clear all
close all

RNAPtot= 4600;
Rtot = 39400;

% general parameters
general.Name              = 'PSM';
general.RNAP_width        = 40; %nucleotides, BN-ID 107873
general.R_width           = 76; %nucleotides, BN-ID 100121, (checked)
general.transcr_speed     = 3300 ; %nucleotides/min, BN-ID 111871 (checked)
general.transl_speed      = 2970; %nucleotides/min, BN-ID 100059 (checked)
general.NA                = 6.0221e+23; %avogadro number
general.Vcell             = 6.7e-16; %cell volume
general.dilution          = 0; %growth rate
general.RNAP_IC           = RNAPtot; % total amount of RNAP at t=0
general.R_IC              = Rtot; % total amount of ribosomes at t=0

% gene parameters
gBG.ID                = 'gBG'; % gene identifier
gBG.numgenes          = 10; % genes, BN-ID 105751 %3000;  %genes, BN-ID 110942
gBG.genelength        = 1064; %nucleotides, BN-ID 105751
gBG.initrate_transcr  = 0.000985; %basal transcription rate
gBG.initrate_transl   = 0.0037; %basal translation rate
gBG.decay_RNA         = 1.518;
gBG.decay_prot        = .01858;
gBG.product           = 'P_BG';
% 
g1.ID                = 'g1'; % gene identifier
g1.numgenes          = 1; 
g1.genelength        = 1064;
g1.initrate_transcr  = 0; 
g1.initrate_transl   = 0.0037;
g1.decay_RNA         = 1.518;
g1.decay_prot        = .01858;
g1.product           = 'P1';
% 
g2.ID                = 'g2'; % gene identifier
g2.numgenes          = 1; 
g2.genelength        = 1064;
g2.initrate_transcr  = 0.000985; 
g2.initrate_transl   = 0.0037;
g2.decay_RNA         = 1.518;
g2.decay_prot        = .01858;
g2.product           = 'P2';


genes=[gBG,g1,g2];

% genetic network interactions
i1.Identifier = 'i1';
i1.Source= 'g1'; % source gene
i1.SourceType= 'Protein'; % Protein or mRNA
i1.Target= 'g1'; % target gene
i1.Mode= 'tx'; % 'tx' or 'tl' (for direct interactions between Proteins or mRNA, use metabolic network)
i1.ParamNames={'k1'};
i1.ParamValues=[.000001];
i1.Fun='-k1 * u1';  %function in parameters and u1, u2 ,..., un 

interactions=[i1];

% % Inputs
u1.Identifier= 'u1';
u1.Target= 'g1';
u1.Mode= 'tx'; % 'tx' or 'tl'
u1.ParamNames={'k1','k2'};
u1.ParamValues=[1e-4,.1];
u1.Fun='k1*cos(k2*t)';

inputs=[u1];

% attention! Units will not be adapted yet!
MN_test=convertSBMLMN('BIOMD0000000051.txt');

% Metabolic network
% example: A->B     B->C   , first reaction catalyzed by P1, second
% spontaneous
MN(1).S           = [-1 0; 1 -1; 0 1];  % Stoichometric matrix (# metabolites x # reactions)
MN(1).Names       = {'cA','cB', 'cC'};  % Name of metabolites 
MN(1).IC          = [1,0,0];
MN(1).ParamNames  = {'vP1_kcatA','vP1_KM','vconstC'};  % Parameter Names 
MN(1).ParamValues = [1e2,0.0017,1e3];  % Parameter values 
MN(1).Fun         = {'vP1_kcatA*cP1*cA / (vP1_KM+cA)','vconstC* cB '};

% direct Interactions between proteins
% example: P1 + P2 <-> P1_P2 with degradation reactions
MN(2).S           = [[-1;-1;1],[1;1;-1],[1;0;-1],[0;1;-1]];  % Stoichometric matrix (# metabolites x # reactions)
MN(2).Names       = {'P1','P2', 'P1_P2'};  % Name of metabolites 
MN(2).IC          = [0,0,0];
MN(2).ParamNames  = {'kfP1P2','kbP1P2'};  % Parameter Names 
MN(2).ParamValues = [1e3,1e2];  % Parameter values 
MN(2).Fun         = {'kfP1P2*P1*P2','kbP1P2*P1_P2','.01858 * P1_P2','.01858 * P1_P2'};


[modelstates,statecoding]=createPSM(general,genes,interactions,inputs,MN);
% [modelstates,statecoding]=createPSM(general,genes,interactions,inputs,[]);


%% model generation

model = IQMmodel([general.Name,'.txtbc']);
IQMmakeMEXmodel(model);
mexfh=eval(['@',general.Name]);

xnames={'RNAP','R'};
ICsoll=[RNAPtot,Rtot]; 
[~,~,XoI]=intersect(xnames,modelstates,'stable');
ICs = IQMinitialconditions(model);
ICs(XoI)=ICsoll;

paramnames=mexfh('parameters');
param=mexfh('parametervalues');
u1paramlog=ismember(paramnames,'u1_P1');

% %% drive into steady state
% pset=param;
% pset(u1paramlog)=0;
% tend=10000;
% simdata = mexfh(linspace(0,tend,20),ICs,pset);
% ICs=simdata.statevalues(end,:);
%% simulation



pset=param;
pset(u1paramlog)=1e-3;
tsim=120;
simdata = mexfh(linspace(0,tsim,200),ICs,pset);


%% plotting
[~,RNAPind]= ismember(modelstates,'RNAP');
RNAPlog=logical(RNAPind);
[~,Rind]= ismember(modelstates,'R');
Rlog=logical(Rind);



tsim=simdata.time;
RNAPsim=simdata.statevalues(:,RNAPlog);
Rsim=simdata.statevalues(:,Rlog);


[~,RNAPoind]= ismember(simdata.variables,'RNAP_o');
RNAPolog=logical(RNAPoind);
[~,Roind]= ismember(simdata.variables,'R_o');
Rolog=logical(Roind);

RNAPactsim=simdata.variablevalues(:,RNAPolog)./RNAPsim;
Ractsim=simdata.variablevalues(:,Rolog)./Rsim;

figure
subplot(2,1,1)
hold all
plot(tsim,RNAPactsim)
title('active RNAP')
subplot(2,1,2)
hold all
plot(tsim,Ractsim)
title('active ribosomes')

figure('Name','geneproducts')
numrow=length(genes);
numcol=3;

for i = 1:length(genes)
    [~,lam]= ismember(simdata.variables,[genes(i).ID,'_init_transcr']);
    lamlog=logical(lam);
    lamsim=simdata.variablevalues(:,lamlog);
    
    PG=statecoding(:,1)==i & statecoding(:,2)==4;
    mrnaPG=statecoding(:,1)==i & statecoding(:,2)==2;
    

    subplot(numrow,numcol,(i-1)*numcol+1)
    plot(tsim,lamsim)
    title('tx initiation rate')
    subplot(numrow,numcol,(i-1)*numcol+2)
    plot(tsim,simdata.statevalues(:,mrnaPG))
    title(modelstates{mrnaPG})
    subplot(numrow,numcol,(i-1)*numcol+3)
    plot(tsim,simdata.statevalues(:,PG))
    title(modelstates{PG})
    
    
end

%%% metabolites and protein complexes
figure('Name','metabolites and protein complexes')
statestoplot=modelstates(statecoding(:,2)==5);
numrow=length(statestoplot);
numcol=1;
for i = 1:length(statestoplot)
    [~,Mind]= ismember(modelstates,statestoplot{i});
    Mlog=logical(Mind);
    subplot(numrow,numcol,i)
    hold all
    plot(tsim,simdata.statevalues(:,Mlog))
    plot([tsim(1),tsim(end)],[1,1]*simdata.statevalues(1,Mlog),'r--')
    title(statestoplot{i})
end
% %%
% g2dnalog=statecoding(:,1)==2&statecoding(:,2)==1;
% g2dnanames=modelstates(g2dnalog);
% xnname=g2dnanames{end};
% g2mrnalog=statecoding(:,1)==2&statecoding(:,2)==3;
% g2mrnanames=modelstates(g2mrnalog);
% ynname=g2mrnanames{end};
% 
% ICset=ICs;
% ICset(XoI)=[5000,Rtot];
% pset=param;
% pset(u1paramlog)=1;
% tsim=120;
% simdata = mexfh(linspace(0,tsim,200),ICset,pset);
% xnlog=ismember(simdata.states,xnname);
% ynlog=ismember(simdata.states,ynname);
% xnmax=simdata.statevalues(end,xnlog);
% RNAPflog= ismember(simdata.variables,'RNAP_f');
% 
% 
% RNAPset=linspace(1,50,10);
% lamset=linspace(0,10,100);
% 
% 
% xnfun=@(x,K)xnmax*x./(x+K);
% xnfun2=@(x,tau)xnmax*tanh(tau*x);
% figure 
% ax=gca;
% hold all
% for k = 1:length(RNAPset)
%     ICset=ICs;
%     ICset(XoI)=[RNAPset(k),Rtot];
%     for i =1:length(lamset)
%         pset=param;
%         pset(u1paramlog)=lamset(i);
%         tsim=120;
%         simdata = mexfh(linspace(0,tsim,200),ICset,pset);
%         xn(i,k)=simdata.statevalues(end,xnlog);
%         yn(i,k)=simdata.statevalues(end,ynlog);
%         X(i,k)=sum(simdata.statevalues(end,g2dnalog));
%         RNAPf(i,k)=simdata.variablevalues(end,RNAPflog);
%     end
%     K=interp1(xn(:,k),lamset,1/8,'pchip');
%     tau(k)=atanh(1/2)/K;
%     
%     ci=ax.ColorOrderIndex;
%     plot(lamset,xn(:,k))
% %     ax.ColorOrderIndex=ax.ColorOrderIndex-1;
% %     plot(P1set,xnfun(P1set,K),'--')
%     ax.ColorOrderIndex=ci;
%     plot(lamset,xnfun2(lamset,tau(k)),'--')
% end
% 
% figure
% plot(RNAPset,tau)
% 
% figure('Name','X: sum over all x_i')
% surf(RNAPset,lamset,X)
% xlabel('RNAPtot')
% ylabel('lambda')
% 
% figure('Name','x_n')
% surf(RNAPset,lamset,xn)
% xlabel('RNAPtot')
% ylabel('lambda')
% 
% figure('Name','y_n')
% surf(RNAPset,lamset,yn)
% xlabel('RNAPtot')
% ylabel('lambda')
% 
% 
