function MN=convertSBMLMN(filename)

MNmodel = IQMmodel(filename);
MNstruct=IQMstruct(MNmodel);

MN(1).S=IQMstoichiometry(MNmodel);
MN(1).Names={MNstruct.states(:).name};
MN(1).IC=[MNstruct.states.initialCondition];
MN(1).ParamNames={MNstruct.parameters.name};
MN(1).ParamValues=[MNstruct.parameters.value];
MN(1).Fun={MNstruct.reactions(:).formula};
