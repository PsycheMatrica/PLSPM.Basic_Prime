function [L, B, Covf] = Dijktra_correction(z0, W0, B0, W, L, Gamma, correct_mode)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dijktra_correction() - MATLAB function to correct bias in the Partial   %
%   Least Sqaures Path Modeling (PLSPM) estimators when constructs are    % 
%   latent variables.                                                     %
% Author: Heungsun Hwang                                                  %
% Contributor: Gyeongcheol Cho                                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments:                                                        %
% W = PLS estimates                                                       %
% L = PLS loadings                                                        %
% Gamma = PLS component scores                                            %
% correct_mode = 1 by P (false for component, true for factor)            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P = size(W0,2);
R = corrcoef(z0);                 % sample correlation matrix
rho = ones(P,1);                % rho_A per block
for p = 1:P
    ind_wp = W0(:,p);
    N_wp=sum(ind_wp,1);    
    Rii = R(ind_wp,ind_wp);
    wp = W(ind_wp,p);
    rho_j = 1;    
    if N_wp > 1 
        if correct_mode(p) == 1
           cj2 = (wp'*(Rii - diag(diag(Rii)))*wp)*pinv(wp'*(wp*wp' - diag(diag(wp*wp')))*wp);
           rho_j = cj2*(wp'*wp).^2;
        end
    end
    if correct_mode(p) == 1
       rho(p) = rho_j;
       L(p,ind_wp) = (sqrt(rho(p))*wp)/(wp'*wp);
    end
end
Covf = eye(P);                  % corrected factor correlations
for q = 1:P
    for p = 1:P
        if q ~= p
           Covf(q,p) = Gamma(:,q)'*Gamma(:,p)/sqrt(rho(q)*rho(p));
        end
    end
end
B = double(B0);                         % path coefficients based on corrected latent correlations
for p = 1:P
    bindexj = find(B0(:,p));
    if ~isempty(bindexj)
       xx = Covf(bindexj,bindexj);
       xy = Covf(p,bindexj)';
       B(bindexj,p) = pinv(xx)*xy;
    end
end
