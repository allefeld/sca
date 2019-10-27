function [rhoEst, r, MSE, E] = sca(R)

% Synchronization Cluster Analysis (SCA)
%
% [rhoEst, r, MSE, E] = sca(R)
%
% R:        matrix of empirical synchronization strengths between oscillators (square)
% rhoEst:   estimate of the synchronization strength between oscillator and cluster
% r:        cluster strength estimate
% MSE:      mean square of normalized errors
% E:        matrix of normalized errors in the reconstruction of R
%
% Code for C. Allefeld and J. Kurths. An approach to multivariate phase synchronization
% analysis and its application to event-related potentials. International Journal of 
% Bifurcation and Chaos, 14(2):417-426, 2004.
%
% Copyright (C) 2003 Carsten Allefeld
%
% This program is free software: you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation, either version 3 of the License, or (at your
% option) any later version. This program is distributed in the hope that
% it will be useful, but without any warranty; without even the implied
% warranty of merchantability or fitness for a particular purpose. See the
% GNU General Public License <http://www.gnu.org/licenses/> for more details.
    
    % number of oscillators
    N = size(R, 1);
    
    % initial estimate
    rhoEst = max(R - diag(diag(R)));
    
    % loop control initialization
    MSEOld = Inf;
    
    for it = 1 : 20
        % variance-dependent weighting factor
        factor = 1 ./ vw(rhoEst' * rhoEst);
        factor = factor - diag(diag(factor));
        
        % mean square of normalized errors
        MSE = mean(mean(factor .* (R - rhoEst' * rhoEst) .^ 2));
        
        % loop control
        if MSE > MSEOld, rhoEst = rhoEstOld; MSE = MSEOld; break, end
        MSEOld = MSE;
        rhoEstOld = rhoEst;

        % new estimate
        rhoEstNext = (rhoEst * (factor .* R)) ./ (rhoEst .^ 2 * factor);
        % smooth out oscillating iteration behavior
        rhoEst = (rhoEst + rhoEstNext) / 2;
    end

    % normalized errors, MSE = mean(mean(E .^ 2))
    factor = 1 ./ vw(rhoEst' * rhoEst);
    factor = factor - diag(diag(factor));
    E = sqrt(factor) .* (R - rhoEst' * rhoEst);
    
    % cluster strength
    A = inv_rhoM(rhoEst);
    r = sum(A .* rhoEst) / sum(abs(A));

        
function ss = vw(rho)
% estimation variance of R as a function of rho

    ss = 1/2 * (1 - rho .^ 2) .^ 2;


function rho = rhoM(kappa)

% population mean resultant length of the von Mises distribution M(mu, kappa)
%
% rho = rhoM(kappa)

    rho = besseli(1, kappa) ./ besseli(0, kappa);


function kappa = inv_rhoM(rho)

% inverse of the rhoM function
%
% kappa = inv_rhoM(rho)
%
% kappa such that rho == rhoM(kappa)

    mlock
    persistent pp rhoMax
    
    if isempty(pp)
        fprintf('inv_rhoM: generating piecewise polynomial\n')
        kappas = [0 :0.01: 10, 10.03 :0.03: 40, 40.1 :0.1: 140];
        rhos = rhoM(kappas);
        pp = spline(rhos, kappas);
        rhoMax = max(rhos);
    end
    
    kappa = ppval(pp, abs(rho)) .* sign(rho);
    
    ind = find(abs(rho) > rhoMax);
    if ~isempty(ind)
        warning(['rho exceeds maximum (' num2str(rhoMax) ')!'])
        kappa(ind) = NaN;
    end
