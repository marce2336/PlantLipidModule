function BI = FindMatchingComps(A, B)
% ==============================================================================
% Cell strings: Find positions of strings of A in B
%
% [AI, BI] = FindMatchingComps(A, B)
%
%  Example: [~, matched] = FindMatchingComps(compNames, compListTemplate{1,2});
%
% INPUT:
%   A:             Cell strings with list of abbreviations for LM compartment to match
%   B:             Cell strings with list of compartments included in the Template Model
%   caseSensitive: If this string starts with 'i', the upper/lower case of the
%                  strings is ignored. Optional, default: 'sensitive'.
% OUTPUT:
%   AI:            Indices of common strings in A as [1 x N] vector.
%                  Each occurence of repeated strings in A is considered.
%                  AI is sorted from low to high indices.
%   BI:            Indices of common strings in B as [1 x N] vector. If B is not unique,
%                  the first occurrence is used.
%                  such that A{AI} == B{BI}.
%
% NOTES:
% - If A or B are multi-dimensional cell arrays, AI and BI are the linear
%   indices (see IND2SUB).
% - If the name of the compartment is composed (e.g. contains more than one
%   word), then each of the entries will be verified in a loop.
%
% ==============================================================================

% ###
% ##### 1). Find out if compartment names are composed:
% ###

c = cellfun(@(x) strsplit(x, ' '), A', 'UniformOutput', false);
isComposed = cellfun(@(x) size(x,2), c, 'UniformOutput', false);
isComposed = str2double(string(isComposed)) > 1;

% ###
% ##### 2). Depending on structure of name, proceed to find coincidences:
% ###

% Note: when synonyms of selected compartment are composed, make sure that the 
%       conincidences found in the Template Model contain an equivalente name.
    
% Find number of elements in lists of compartments:
nB = numel(B);
nA = numel(A);


% Collect the index of the first occurrence in B for every A:
M = zeros(1, nB);
copyA = A;
count = 1;

for iA = 1:nA

 switch sum(isComposed) > 0 %---------------------------------------------------
     case 0 % Compartment is identified with a single word
         findSynonyms = find(contains(B, A{iA}, 'IgnoreCase',true));
         %idx = findSynonyms(strlength(B(findSynonyms)) <= max(sizeA)+2);

     case 1 % Compartment name contains more than two words
         comp_i = split(A{iA}, ' ');

         switch numel(comp_i) %-------------------------------------------------
             case 1
                 findSynonyms = find(contains(B, comp_i{1}, 'IgnoreCase',true));

                 %idx = findSynonyms(strlength(B(findSynonyms)) <= max(sizeA));

             case 2
                 findSynonyms = find(contains(B, comp_i{1}, 'IgnoreCase',true).*...
                     contains(B, comp_i{2}, 'IgnoreCase',true));

                 %idx = findSynonyms(strlength(B(findSynonyms)) <= max(sizeA));
                 
             case 3
                 findSynonyms = find(contains(B, comp_i{1}, 'IgnoreCase',true).*...
                     contains(B, comp_i{2}, 'IgnoreCase',true).*...
                     contains(B, comp_i{3}, 'IgnoreCase',true));
                 %idx = findSynonyms(strlength(B(findSynonyms)) <= max(sizeA));
         end

 end

 if sum(findSynonyms) ~= 0 %sum(idx) ~= 0
    M(count:(numel(findSynonyms)+count)-1) = findSynonyms; %idx;
    count = count+numel(findSynonyms);
 else
     copyA{iA} = 'NA';
 end
end


% ###
% ##### 3). Obtain indexes of matching compartment:
% ###

matchingB = unique(M(M~=0));  % If any occurrence was found in A, get index in B
matchedNames = B(matchingB);

% ###
% ##### 4). Double check the identity of matching compartment:
% ###

copyA(strcmp(copyA, 'NA')) = ''; % Eliminate compartments not matched
sizeA = strlength(copyA);

doubleCheck = strlength(matchedNames) <= max(sizeA)+2;
BI = matchingB(doubleCheck);

end
