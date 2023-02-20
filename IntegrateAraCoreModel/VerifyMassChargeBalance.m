function RxnBalancingList = VerifyMassChargeBalance(model)
%==========================================================================
% 'VerifyMassChargeBalance' function creates a list withh all reactions
% which are mass- and/or-charge-imbalanced, excluding reactions naturally
% imbalanced such as exchange, import, siphon, demand, biomass and SLIME
% reactions.
%
% USAGE:
%
%       RxnStoichiometry = VerifyUnbalancedRxns(UnbalancedRxnsList,model)
%
% INPUTS:       
%       model:              COBRA model structure containing the AraCore 
%                           model published by Arnold and Nikoloski (2014)
%
% OUTPUTS:      
%       RxnBalancingList:   List containing all the reaction which are
%                           mass- and/or charge-imbalanced reaction,
%                           excluding reactions naturally imbalanced such
%                           as exchange, import, siphon, demand, biomass
%                           and SLIME reactions.
%==========================================================================

[~, imBalancedMass, imBalancedCharge, imBalancedRxnBool, ~, ~, ~] = checkMassChargeBalance(model, 1);

RxnBalancingList = [string(model.rxns), string(imBalancedMass), string(imBalancedCharge), string(imBalancedRxnBool)];

% Eliminate mass- and charge-balanced reactions:
IdxBalanced = strcmp(RxnBalancingList(:,4), 'false');
RxnBalancingList = RxnBalancingList(IdxBalanced == 0,:);

% Eliminate reactions naturally imbalanced: exchange, import, siphon, 
% demand, biomass and SLIME reactions 
IdxNatImb = contains(RxnBalancingList(:,1), {'EX_','Ex_','Im_','Si_H','biomass','SLIMEr','DM','Bio_AA','Bio_CLim','Bio_NLim','Bio_opt'});
RxnBalancingList = RxnBalancingList(IdxNatImb == 0,:);

% Eliminate reactions where metabolites without known molecular formula are
% participating:
IdxNaN = strcmp(RxnBalancingList(:,2), 'NaN');
RxnBalancingList = RxnBalancingList(IdxNaN == 0,:);

end