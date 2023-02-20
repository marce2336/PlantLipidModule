function [indexA, indexB] = CompareArrays(arrayA, arrayB, flagCase)
%% ========================================================================
% 'CompareArrays' compares the elements of array A and B, and finds the
% index of common elements. 
%
% USAGE: 
%  [indexA, indexB] = CompareArrays(arrayA, arrayB, flagCase)
%
% INPUTS:
%       arrayA, arrayB: Cell strings of any size.
%
%       flagCase:       (i) - Ignore the upper/lower case of strings during
%                             the comparison
%                       (s) - Default. The comparison is 'case sensitive'
%         
% OUTPUT:
%       indexA:         Index of common elements in arrayA
%       indexB:         Index of common elements in arrayB
%
% =========================================================================

if ~exist('flagCase', 'var')
    flagCase = 's';
end

% ###
% ##### 1). For each element in array A, search for the first occurrence in array B:
% ###
nA = numel(arrayA);
nB = numel(arrayB);

if nargin > 2 && strncmpi(flagCase, 'i', 1)  % Ignore the upper/lower case
   if nA < nB
      % Retrieve the index of the first common element in array B for every element of array A
      commonElements = zeros(1, nA);
      for iA = 1:nA
         Ind = find(strcmpi(arrayA{iA}, arrayB));
         if ~isempty(Ind)
            commonElements(iA) = Ind(1);
         end
      end
      
   else  % When array-B <= array-A, B, iterate over B
      commonElements = zeros(1, nA);
      for iB = nB:-1:1
         commonElements(strcmpi(arrayB{iB}, arrayA)) = iB;
      end
   end
   
else  % Consider the upper/lower case
   if nA <= nB
      % Retrieve the index of the first common element in array B for every element of array A
      commonElements = zeros(1, nA);
      for iA = 1:nA
         Ind = find(strcmp(arrayA{iA}, arrayB));
         if ~isempty(Ind)
            commonElements(iA) = Ind(1);
         end
      end
      
   else  % nB <= nA, B is smaller, so better loop over B:
      % Mark every A which equal the current B
      commonElements = zeros(1, nA);
      for iB = nB:-1:1
         commonElements(strcmp(arrayB{iB}, arrayA)) = iB;  % Same as M(find(strcmp()))
      end
   end
end


% ###
% ##### 2). Get indexes of common elements:
% ###
indexA = find(commonElements); % Common elements in array A found
indexB = commonElements(indexA); % Index of common elements of array-A in array-B

end
