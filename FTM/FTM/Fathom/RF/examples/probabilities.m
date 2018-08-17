D = 1:5
S = (max(D) - D) + min(D)
L = S/sum(S)
PP = 1 - (D/max(D))

S        = (max(D(:)) - D) + min(D(:)); % convert distance to similarity
lik(i,:) = S/sum(S);                    % normalized similarities
