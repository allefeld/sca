function l = simplexGuess(X)

% determine nodes of a simplex approximating the data
%
% l = simplexGuess(X)
%
% X:    data points (points x dimensions)
% l:    indices of points comprising the simplex
%
% Following: P. Deuflhard and M. Weber, "Robust Perron cluster analysis in
% conformation dynamics", Linear Algebra and its Applications 398 (2005),
% pp. 161-184. doi:10.1016/j.laa.2004.10.026
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


    [N, dim] = size(X);

    % point farthest from origin
    [dummy, l] = max(sum(X .^ 2, 2));
    % translate it to the origin
    X0 = X - ones(N, 1) * X(l, :);

    for i = 1 : dim
        % point farthest from origin
        [dummy, ln] = max(sum(X0 .^ 2, 2));
        l = [l, ln];
        % remove corresponding subspace
        v = X0(ln, :);
        v = v / sqrt(sum(v .^ 2));
        X0 = X0 - X0 * v' * v;
    end
