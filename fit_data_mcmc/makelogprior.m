function [logprior] = makelogprior(lb,ub,nonlcon)
% This function will return logprior, a logical function which will be true
% if evaluated for a feasible walker.

% Inputs:
%   lb - lower bound vector for all parameters
%   ub - upper bound vector for all parameters
%   nonlcon - nonlinear constraints for the model of interest which takes a
%   walker as an input

% Output:
%   logprior - function handle which can be evaluated for a set of walkers.

% Check bounds are the same length
Nparams = length(lb);
if Nparams ~= length(ub)
    error('Upper and lower bounds are not the same length')
end

logprior = @(x) logical(true);

for i = 1:Nparams
    additional_constraint = @(x) x(i) <= ub(i) & x(i) >= lb(i);
    logprior = @(x) logprior(x) & additional_constraint(x);
end

% Add nonlinear constraints
if ~isempty(nonlcon)
    logprior = @(x) logprior(x) && nonlcon(x);
end

end

