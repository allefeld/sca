function [ci, qsRanked, cs, cp] = scaM(R, q)

% synchronization cluster analysis via Markov coarse-graining
%
% [ci, qsRanked, cs, cp] = scaM(R)
% [ci, qsRanked, cs, cp] = scaM(R, q)
%
% R:            matrix of bivariate synchronization strengths
%               (values from [0, 1], symmetric)
% q:            number of clusters, optional
%       q is normally determined from R. Via this parameter,
%       the automatic selection can be overridden.
%
% ci:           cluster indices, for each element the cluster it belongs to
% qsRanked:     the possible q-values in ranking order
%
% optional additional output:
%       simply defined cluster strengths and participation indices
% cs:           cluster strengths of identified clusters
% cp:           corresponding cluster participation index vectors
%
% If no output parameters are specified, a summary of the analysis results
% is displayed.
%
% Code for: C. Allefeld and S. Bialonski, "Detecting synchronization clusters
% in multivariate time series via coarse-graining of Markov chains", submitted.
% arXiv:0707.2479
%
% Copyright (C) 2007 Carsten Allefeld
%
% This program is free software: you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation, either version 3 of the License, or (at your
% option) any later version. This program is distributed in the hope that
% it will be useful, but without any warranty; without even the implied
% warranty of merchantability or fitness for a particular purpose. See the
% GNU General Public License <http://www.gnu.org/licenses/> for more details.


    % check input
    if (min(R(:)) < 0) || (max(R(:)) > 1)
        error('scaM: synchronization strengths have to be from [0, 1]!')
    end
    if norm(R - R') > norm(R) * eps
        error('scaM: synchronization matrix has to be symmetric!')
    end


    % transition matrix of a reversible Markov process
    P = R * diag(1 ./ sum(R));


    % eigenvalue decomposition, left eigenvectors
    [A, l] = eig(P');
    l = diag(l);
%     [dummy, ind] = sort(abs(l), 'descend');
    [dummy, ind] = sort(-abs(l));
    l = l(ind);
    A = A(:, ind);

    % stationary measure
    p0 = sum(R)' / sum(sum(R));

    % normalize left eigenvectors
    A = A * diag(1 ./ sqrt(sum(diag(p0) * A .^ 2)));


    % modulus of eigenvalues
    al = abs(l);

    % regularize small eigenvalues close to numerical precision
    al(al < eps) = eps;

    % compute timescales and separation factors
    zeta = 0.01;
    tau = [Inf; log(zeta) ./ log(al(2 : end))];
    F = tau(1 : end - 1) ./ tau(2 : end);

    % timescale-based ranking of q-values
%     [dummy, qsRanked] = sort(F, 'descend');
    [dummy, qsRanked] = sort(-F);

    % if q is not specified explicitely, select it based on the ranking
    if nargin < 2
        q = qsRanked(2);
    end

    % determine timescale
    tauq = tau(q + 1);


    % eigenvector space for selected timescale
    o = A * diag(al .^ tauq);
    o = o(:, 2 : q);


    % k-means clustering, initialized using approximate simplex vertices
    ci = kMeans(o, o(simplexGuess(o), :));


    % additionally: simply defined cluster strengths and participation indices
    cp = zeros(size(R, 1), q);
    for i = 1 : q
        ind = (ci == i);
        cp(ind, i) = mean(R(ind, ind), 2);
    end
    cs = sum(cp);

    % remap such that cluster indices correspond to cluster strengths
%     [cs, ind] = sort(cs, 'descend');
    [cs, ind] = sort(-cs); cs = -cs;
    cp = cp(:, ind);
    [dummy, map] = sort(ind);
    ci = map(ci);


    % diagnostic output
    if nargout == 0
        fprintf('scaM\n\n')
        fprintf(' ranking of possible numbers of clusters:\n');
        disp([['  '; '  '; '  '], num2str([qsRanked'; tau(qsRanked + 1)'; F(qsRanked)'], ' %7.3g')])
        fprintf('\n');
        fprintf(' looking for %d clusters at timescale %.3g\n\n', q, tauq)
        fprintf(' clustering:\n')
        fprintf('  %d', ci);
        fprintf('\n')
    end
