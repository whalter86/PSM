function flag = writeToFile(header)
try
    headertext=header.headertext;
    geneprodtext=header.geneprodtext;
    statetext=header.statetext;
    paramtext=header.paramtext;
    vartext=header.vartext;
    reactiontext=header.reactiontext;
    functiontext=header.functiontext;
    footertext=header.footertext;
    
%     general
%     modelstates
%     generalparamnames
%     geneparamnames,geneparamvals
%     interactionparamnames,interactionparamvals,
%     paramtext
%     statetext
    
    for i = 1:length(modelstates)
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