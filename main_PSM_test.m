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
gBG.numgenes          = 10; % genes, BN-ID 105751 %3000;  %genes, BN-ID 110942
gBG.genelength        = 1064; %nucleotides, BN-ID 105751
gBG.initrate_transcr  = 0.000985; %basal transcription rate
gBG.initrate_transl   = 0.0037; %basal translation rate
gBG.decay_RNA         = 1.518;
gBG.decay_prot        = .01858;
gBG.product           = 'P_BG';
% 
g1.numgenes          = 1; 
g1.genelength        = 1064;
g1.initrate_transcr  = 0.000985; 
g1.initrate_transl   = 0.0037;
g1.decay_RNA         = 1.518;
g1.decay_prot        = .01858;
g1.product           = 'P1';



genes=[gBG,g1];

% genetic interactions
% i1.Source= 'P_G2';
% i1.SourceType= 'Protein'; % Protein or mRNA
% i1.Target= 'P_G3';
% i1.Mode= 'tx'; % 'tx' or 'tl'
% i1.V = -.015;
% i1.h = 2.1;
% i1.K = 40;

interactions=[];

% % Inputs
u1.Identifyer= 'u1';
u1.Target= 'P1';
u1.Mode= 'tx'; % 'tx' or 'tl'
u1.ParamNames={'P1','P2'};
u1.ParamValues=[1e-3,0];
u1.Fun='P1*cos(P2*t)';

inputs=[u1];

[modelstates,statecoding]=createPSM(general,genes);
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

%% drive into steady state
pset=param;
pset(u1paramlog)=0;
tend=10000;
simdata = mexfh(linspace(0,tend,20),ICs,pset);
ICs=simdata.statevalues(end,:);
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
    [~,lam]= ismember(simdata.variables,['g',num2str(i),'_init_transcr']);
    lamlog=logical(lam);
    lamsim=simdata.variablevalues(:,lamlog);
    
    PG=statecoding(:,1)==i & statecoding(:,2)==4;
    mrnaPG=statecoding(:,1)==i & statecoding(:,2)==2;
    

    subplot(numrow,numcol,(i-1)*numcol+1)
    plot(tsim,lamsim)
    title('tx initiation rate')
    subplot(numrow,numcol,(i-1)*numcol+2)
    plot(tsim,simdata.statevalues(:,mrnaPG))
    title('mRNA amount')
    subplot(numrow,numcol,(i-1)*numcol+3)
    plot(tsim,simdata.statevalues(:,PG))
    title('protein amount')
    
    
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
