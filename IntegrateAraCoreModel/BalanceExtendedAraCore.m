function BalanceExtendedAraCore
%==========================================================================
%
% This function uses as input the AraCore model that was extended with the
% Plant Lipid Module. 
% The first step consists in removing the reactions with duplicated
% stoichiometry. Next, the mass- and charge-imbalanced reactions are
% identified and adjusted. 
% The output of the function is a model named 'OutputModelUniqueBalanced'
% that is a mass- and charge- balanced version of the model extended with
% the Plant Lipid Module.
% The function also prints out a report with a list of reactions that were
% not possible to balance in the case that they exist, named
% 'RxnBalancingListIntegrated'.
%
% USAGE:
%       BalanceExtendedAraCore
%
%==========================================================================

% Load .mat integrated model containing unique reactions:
pathIntegModel = what('LipidModuleIntegration');
pathInputModelFile1 = fullfile(pathIntegModel.path, 'OutputFiles', 'OutputModelUnique.mat');
OutputModelUnique = readCbModel(pathInputModelFile1);


% ###
% ##### 1). The usage of M-ACP[h] has to be harmonized in the model:
% ### 
% Find metabolite index:
MetIndx = strcmp(OutputModelUnique.mets, 'LM_M-ACP[h]');
GetSindx = OutputModelUnique.S(MetIndx, :);

% Find non-zero entries in S and add to M-ACP[h]:
GetCoefficients = OutputModelUnique.S(MetIndx, GetSindx ~= 0);
MetIdx = strcmp(OutputModelUnique.mets, 'M-ACP[h]'); % Find M-ACP[h] index
OutputModelUnique.S(MetIdx, GetSindx ~= 0) = GetCoefficients;

% Set all coefficients to zero and remove metabolite from model:
OutputModelUnique.S(MetIndx, :) = 0;
OutputModelUnique = removeMetabolites(OutputModelUnique, 'LM_M-ACP[h]');


% ###
% ##### 2). Eliminate duplicated light reactions:
% ###
RmvRxnsList = {'LM_PSI_h','LM_PSII_h','Ex_pho_loss'};
OutputModelUnique = removeRxns(OutputModelUnique, RmvRxnsList);


% ###
% ##### 3). Eliminate unused metabolite:
% ###
OutputModelUnique = removeMetabolites(OutputModelUnique, 'LM_pho_loss[h]');


% ###
% ##### 4). Double-check for duplicated reactions:
% ###
% FR –> checks F + R matrix, where S:=?F+R, which ignores reaction direction
[OutputModelUnique2, removedRxnInd, keptRxnInd] = checkDuplicateRxn(OutputModelUnique, 'FR', 1, 1);

%Adjust reversibility of reactions removed: 
if ~isempty(removedRxnInd)
    RxnsToAdjust = CheckRxnDirectionality(keptRxnInd, removedRxnInd, OutputModelUnique);
    for j = 1:size(RxnsToAdjust,1)
        KeptRxnID = RxnsToAdjust{j};
        KeptRxnIdx = strcmp(OutputModelUnique2.rxns, KeptRxnID);
        OutputModelUnique2.lb(KeptRxnIdx) = -1000;
    end
end


% ###
% ##### 5). Double-check for unbalanced reactions in the integrated model and print report:
% ###
RxnBalancingListIntegrated = VerifyMassChargeBalance(OutputModelUnique2);


% ###
% ##### 6). If there are still mass- and/or charge-unbalanced reactions act accordingly:
% ###
if ~isempty(RxnBalancingListIntegrated)
    OutputModelUnique2 = AdjustUnbalancedAtoms(OutputModelUnique2,RxnBalancingListIntegrated);
end

ModelBalanced = OutputModelUnique2;


% ###
% ##### 7). Save integrated model into the destination folder in .xls, .mat and .xml formats:
% ###
dirSave = what('IntegrateAraCoreModel');
dirSave = dirSave.path;
pathSave = fullfile(dirSave, 'OutputFiles');
oldPath = cd(pathSave);

writeCbModel(ModelBalanced,'format','XLS','fileName','OutputModelUniqueBalanced.xls')
writeCbModel(ModelBalanced,'format','MAT','fileName','OutputModelUniqueBalanced.mat')
%writeCbModel(ModelBalanced,'format','SBML','fileName','OutputModelUniqueBalanced.xml')
save('RxnBalancingListIntegrated', 'RxnBalancingListIntegrated')

cd(oldPath)

fprintf('============================================================================\n')
fprintf('The mass- and charge-balanced model is saved in the following route:\n')
fprintf('''../IntegrateAraCoreModel/OutputFiles''!\n')
fprintf('============================================================================\n')

end