function flag = writeToFile(general,geneproduction,states,parameters,variables,reactions,functions)
try
    %% filename
    if length(general.Name)<5 || ~strcmp(general.Name(end-5:end),'.txtbc')
        dotind=regexp(general.Name,'\.');
    %     dotind=strfind(general.Name,'.');
        if isempty(dotind)
            dotind=length(general.Name)+1;
        end
        filename=[general.Name(1:dotind-1),'.txtbc'];
    else
        filename=general.Name;
    end


    %% define standard blocks
    headertext={'********** MODEL NAME \n';...
    [general.Name,' \n'];...
    '********** MODEL NOTES \n';...
    '********** MODEL STATE INFORMATION \n'};
    geneprodtext={};
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

    %% append network specific lines
    for i = 1:length(geneproduction.names)
        geneprodtext=[geneprodtext;['d/dt(',geneproduction.names{i},')=',geneproduction.rhs{i},' \n']];
    end
    
    for i = 1:length(states.names)
        statetext=[statetext;[states.names{i},'(0)=',num2str(states.IC(i)),' \n']];
    end
    
    for i = 1:length(parameters.names)
        paramtext=[paramtext;[parameters.names{i},'=',num2str(parameters.vals(i)),' \n']];
    end
    
    for i = 1:length(variables.names)
        vartext=[vartext;[variables.names{i},'=',variables.vals{i},' \n']];
    end
    
    for i = 1:length(reactions.names)
        validatestring(reactions.operators{i},{'=>','<=>','<='});
        reactiontext= [reactiontext;[reactions.educts{i},' ', reactions.operators{i},' ' ,reactions.products{i},' :', reactions.names{i},' \n' ]];
        if strcmp(reactions.operators{i},'=>')
            reactiontext= [reactiontext; ['vf=',reactions.forwardrate{i}, ' \n']];
        elseif strcmp(reactions.operators{i},'<=>')
            reactiontext= [reactiontext; ['vf=',reactions.forwardrate{i}, ' \n']];
            reactiontext= [reactiontext; ['vr=',reactions.backwardrate{i}, ' \n']];
        elseif strcmp(reactions.operators{i},'<=')
            reactiontext= [reactiontext; ['vr=',reactions.backwardrate{i}, ' \n']];
        end
    end
    
    for i = 1:length(functions.names)
        functiontext=[functiontext;[functions.names{i},'(',functions.arguments{i},')=',functions.vals{i},' \n']];
    end
    

    %% write to file
    alltext=[headertext;geneprodtext;statetext;paramtext;vartext;reactiontext;functiontext;footertext];
    fid = fopen(filename,'w');
    for row = 1:size(alltext,1)
        fprintf(fid,alltext{row,:});
    end
    fclose(fid);
    flag = true;
catch
    flag = false;
end