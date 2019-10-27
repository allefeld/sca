function [cs, cp, Rt, e, v] = scaE(R)

% generalized Synchronization Cluster Analysis via eigenvalue decomposition
%
% [cs, cp, Rt, e, v] = scaE(R)
%
% R:    synchronization matrix
%
% cs:   cluster strengths of identified clusters,
%       in ascending order (row vector)
% cp:   corresponding cluster participation index vectors
%       one column per cluster
%
% additional data:
% Rt:   trimmed synchronization matrix
%  results of the original eigenvalue decomposition
% e:    eigenvalues larger than 1, in ascending order (row vector)
% v:    corresponding normalized eigenvectors, one column per eigenvalue
%
% Code for C. Allefeld, M. Müller, and J. Kurths. Eigenvalue decomposition as a
% generalized synchronization cluster analysis. International Journal of 
% Bifurcation and Chaos, 17, 2007 (in press).
%
% Copyright (C) 2005 Carsten Allefeld
%
% This program is free software: you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation, either version 3 of the License, or (at your
% option) any later version. This program is distributed in the hope that
% it will be useful, but without any warranty; without even the implied
% warranty of merchantability or fitness for a particular purpose. See the
% GNU General Public License <http://www.gnu.org/licenses/> for more details.


    N = size(R, 1);
    
    % eigenvalue decomposition of original matrix
    [v, d] = eig(R);
    e = diag(d)';
    [e, ind] = sort(e);
    v = v(:, ind);

    % select components with eigenvalues larger than 1
    ind = find(e > 1);
    e = e(ind);
    v = v(:, ind);
    if size(ind, 2) == 0,
%         fprintf('No relevant eigenvectors found!\n')
        cs = [];
        cp = [];
        Rt = R;
        return
    end

    % identify clusters
    % and remove inter-cluster synchronizations from the matrix
%     global method
%     switch method
%         case 'sign'     % pattern of eigenvector signs
%             id = (sign(v) == 1) * 2 .^ (0 : size(v, 2) - 1)';
%             uid = unique(id)';
%         case 'maxpart'  % maximum participation in cluster
            [dummy, id] = max(ones(N, 1) * e .* v .^ 2, [], 2);
            uid = 1 : size(v, 2);
%             % ! in this and the next case, uid is not really the list of
%             % ! unique ids, because it may contain ids that do not exist
%         case 'maxdist'  % maximum contribution to eigenvector
%             [dummy, id] = max(v .^ 2, [], 2);
%             uid = 1 : size(v, 2);
%         case 'none'     % no selection
%             id = ones(N, 1);
%             uid = 1;
%     end
    
    NC = size(uid, 2);      % number of clusters
    Rt = R;
    for i = 1 : NC
        Rt(id == uid(i), id ~= uid(i)) = 0;
        Rt(id ~= uid(i), id == uid(i)) = 0;
    end

    % repeat eigenvalue decomposition for the trimmed matrix
    [vt, d] = eig(Rt);
    et = diag(d)';
    [et, ind] = sort(et);
    vt = vt(:, ind);

    % compute cluster participation indices
    % and select "proper" clusters
    cp = ones(N, 1) * et .* vt .^ 2;
    ind = [];
    for i = 1 : NC
        rcs = sum(cp(id == uid(i), :), 1);
        [maxval, j] = max(rcs);
        if maxval > 1
            ind = [ind, j];
        end
    end
    cp = cp(:, ind);
    cs = sum(cp);
    
