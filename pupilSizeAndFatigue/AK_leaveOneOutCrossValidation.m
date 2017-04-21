function [ accuracy ] = AK_leaveOneOutCrossValidation( predictors, classification )
%AK_leaveOneOutCrossValidation cross-validates the accuracy of linear
%discriminant analysis for a set of predictors and set of classifications
%   INPUT:
%       predictors: a matrix of predictor variables n observations by n predictor variables in size;
%       classification: a vector of classifications for each observation (should be a cell array of strings)
%   OUTPUT:
%       accuracy: this cross validation cycles through observations one at
%       a time, leaving one out of the algorithm training and trying to
%       classify it based on the others; this accuracy statistic represents
%       the proportion of observations the algorithm correctly classifies
%       based on all the other observations


% check inputs
if nargin < 2 || ~iscell(classification) || ~all(cellfun(@ischar,classification))
    error(['AK_leaveOneOutCrossValidation requires two arguments as input: '... 
        'predictors: a matrix of predictor variables n observations by n predictor variables in size; '...
        'classification: a vector of classifications for each observation (should be a cell array of strings)']);
end

% check accuracy of leave-one-out validation
for iO = 1:length(classification) % cycle through observations
    % leave out the current observation
    test = predictors(iO,:);
    ansClass = classification(iO);
    training = predictors; training(iO,:) = [];
    trainingClass = classification; trainingClass(iO) = [];
    % run linear discriminant analysis
    try
        guessClass = classify(test,training,trainingClass,'linear');
    catch
        guessClass = 'LDA failed';
    end
    % grade performance for this observation
    correct(iO) = strcmp(ansClass,guessClass);
end

% calculate accuracy
accuracy = sum(correct)/length(correct);

end
