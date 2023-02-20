function IntegrateAraCoreModel
%==========================================================================
%
% 'IntegrateAraCoreModel' function performs the integration of the Plant
% Lipid Module into the AraCore model published by Arnold and Nikoloski
% (2014) Plant Physiol. 165(3):1380-1391 (doi: 10.1104/pp.114.235358).
% The function also prints out a report with a list of unbalanced reactions
% named 'UnbalancedRxnsAraCore'.
%
% USAGE:
%       IntegrateAraCoreModel
%
%==========================================================================

% ###
% ##### 1). Verify the existence of AraCore model in the directory:
% ###

% Find file directory
findDir = what('ExampleModel');
getPath = fullfile(findDir.path);
pathModel = fullfile(getPath, 'AraCore.mat');

% Verify if AraCore model file is in directory:
ExistAraCoreModel = exist(pathModel, 'file');

if ExistAraCoreModel == 0
    
    fprintf('The adjusted AraCore model doesn´t exist!. \n')
    fprintf('Please wait until it is created in the following route ''../ExampleModel''!\n')

    %Run function to create adjusted AraCore model in folder:
    model = AdjustOriginalAraCore;

    oldPath = cd(getPath);
    writeCbModel(model,'format','XLS','fileName','AraCore.xls')
    writeCbModel(model,'format','MAT','fileName','AraCore.mat')
    cd(oldPath)
 
% If adjusted AraCore was already created it must be loaded into the workspace:    
else
    model = readCbModel(pathModel); 
end


% ###
% ##### 2). Search for unbalanced reactions and print report:
% ###

RxnBalancingList = VerifyMassChargeBalance(model);
UnbalancedRxnsAraCore = VerifyUnbalancedRxns(RxnBalancingList,model);

dirSave = what('IntegrateAraCoreModel');
dirSave = dirSave.path;
pathSave = fullfile(dirSave, 'OutputFiles', 'UnbalancedRxnsAraCore');
save(pathSave, "UnbalancedRxnsAraCore")
   

% ###
% ##### 3). Verify if the Plant Lipid Module reactions were already integrated into the AraCore model:
% ###

pathIntegModel = what('LipidModuleIntegration');
ExistIntModel = exist(fullfile(pathIntegModel.path, 'OutputFiles', 'OutputModelUnique.mat'), 'file');

switch ExistIntModel
    % Integrate the LipidModule into AraCore model:
    case 0
        fprintf('=====================================================================================================\n');
        fprintf('The ''LipidModuleIntegration'' function is going to be launched, so please follow the instructions!\n')
        fprintf('=====================================================================================================\n');
        LipidModuleIntegration
     
    % The LipidModule was already integrated!
    otherwise
        fprintf('=====================================================================================================\n');
        fprintf('The Plant Lipid Module was already integrated into the AraCore model!. \n')
        fprintf('Search the extended model in the following route ''../LipidModuleIntegration/OutputFiles''!\n')
        fprintf('=====================================================================================================\n');
    
        cd(oldPath)        
end

end