function ci = kMeans(x, cc, maxIter)

% simple implementation of the k-means algorithm
%
% ci = kMeans(x, cc, maxIter = 100)
%
% x:        data points (n x d)
% cc:       initial cluster centers (k x d)
% ci:       cluster indices (n)
% maxIter:  maximum number of iterations
%
% Separates a data set x of n points in a d-dimensional Euclidean space into
% k clusters, trying to minimize the sum of within-cluster variances. The
% algorithm starts from the k initial cluster centers given in cc.
%
% J. B. MacQueen, "Some Methods for Classification and Analysis of
% Multivariate Observations", Proceedings of 5th Berkeley Symposium on
% Mathematical Statistics and Probability (1967), pp. 281-297.
% See also <http://en.wikipedia.org/wiki/K-means>
%
% helper function for scaM.m
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


    if nargin < 3, maxIter = 100; end

    n = size(x, 1);
    k = size(cc, 1);

    d = zeros(n, k);
    for iter = 1 : maxIter
        % compute distances of all points to all cluster centers
        for i = 1 : k
            d(:, i) = sum((x - ones(n, 1) * cc(i, :)) .^ 2, 2);
        end

        % assign points to cluster with minimum distance
        [dummy, ci] = min(d');

        % if a cluster is empty, the data point farthest
        % from its cluster center is assigned to it
        dc = diag(d(1 : n, ci));
        reassigned = true;
        while reassigned
            reassigned = false;
            for i = 1 : k
                if ~any(ci == i)
                    [dummy, outlier] = max(dc);
%                     fprintf('empty cluster %d, assigned %d\n', i, outlier)
                    ci(outlier) = i;
                    dc(outlier) = 0;
                    reassigned = true;
                end
            end
        end

        % cluster centers are moved to mean of assigned points
        for i = 1 : k
            cc(i, :) = mean(x(ci == i, :));
        end

        % cluster assignments change no longer: stop
        if exist('ciOld')
            if all(ciOld == ci), return, end
        end
        ciOld = ci;
    end
    error('kMeans: no solution found within %d iterations', maxIter)
