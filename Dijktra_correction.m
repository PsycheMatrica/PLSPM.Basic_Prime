function [W,C, B, Cov_CV, rho] = Dijktra_correction(Z, W0, B0, W, C, Gamma, correct_mode,Cov_CV_s1,Rho_A_s1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dijktra_correction() - MATLAB function to correct bias in the Partial   %
%   Least Sqaures Path Modeling (PLSPM) estimators                        %
% Author: Gyeongcheol Cho & Heungsun Hwang                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments:                                                        %
%   Z            - Data on standardized indicators                        %
%   W0           - Logical matrix for weights                             %
%   B0           - Logical matrix for path coefficients                   %
%   W            - Matrix of weight estimates                             %
%   C            - Matrix of loading estimates                            %
%   Gamma        - PLSPM component scores                                 %
%   correct_mode - 1 x P2 vector of correction strategies)                %
%   Cov_CV_s1    - Consistent covariance matrix of 1st-order constructs   %    
%   Rho_A_s1     - 1 × P1 vector of Dijktra's construct reliabilities     % 
% Output arguments:                                                       %
%   W      - Matrix of corrected weight estimates                         %
%   C      - Matrix of corrected loading estimates                        %
%   B      - Matrix of corrected path coefficient estimates               %
%   Cov_CV - Construct covariance matrix                                  %
%   rho    - 1 × P vector of construct reliabilities                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P = size(W0,2);
S = Z'*Z; %corrcoef(z0);              
rho = ones(1,P);               
for p = 1:P
    ind_wp = W0(:,p);
    N_wp=sum(ind_wp,1);    
    S_pp = S(ind_wp,ind_wp);
    wp = W(ind_wp,p);
    if N_wp > 1 
        if correct_mode(1,p) == 1 % 1st-order factor with observed indicators
                                % or 2nd-order factor with 1st-order components  
           cj2 = (wp'*(S_pp - diag(diag(S_pp)))*wp)*pinv(wp'*(wp*wp' - diag(diag(wp*wp')))*wp);
           rho(1,p) = cj2*(wp'*wp).^2;

           cp=(sqrt(rho(1,p))*wp)/(wp'*wp);
           C(p,ind_wp) = cp; 
        elseif correct_mode(1,p) == 2 % 2nd-order component with 1st-order factors  
           rho_a=Rho_A_s1(1,ind_wp);
           S_pp(eye(N_wp)==1)=rho_a;
           rho(1,p) = wp'*S_pp*wp;   

           q_p=C(p,ind_wp)./sqrt(rho_a);
           Sp=Cov_CV_s1(ind_wp,ind_wp);
           wp=Sp\(q_p');
           wp_std=(1/sqrt(wp'*Sp*wp))*wp;
           W(ind_wp,p)=wp_std;
           C(p,ind_wp)=wp_std'*Sp;
        elseif correct_mode(1,p) == 3 % 2nd-order factor with 1st-order factors  
           rho_a=Rho_A_s1(1,ind_wp);
           S_pp(eye(N_wp)==1)=rho_a;
           rho(1,p) = wp'*S_pp*wp;  

           cp=(sqrt(rho(1,p))*wp)/(wp'*wp);
           C(p,ind_wp) = cp; 
        end
    else
        if correct_mode(1,p) == 2 || correct_mode(1,p) == 3 
           rho(1,p) = Rho_A_s1(1,ind_wp);
        end
    end
end
Cov_CV = eye(P);                  
for q = 1:P
    for p = 1:P
        if q ~= p
           Cov_CV(q,p) = Gamma(:,q)'*Gamma(:,p)/sqrt(rho(1,q)*rho(1,p));
        end
    end
end
B = double(B0);                        
for p = 1:P
    bindexj = find(B0(:,p));
    if ~isempty(bindexj)
       xx = Cov_CV(bindexj,bindexj);
       xy = Cov_CV(p,bindexj)';
       B(bindexj,p) = pinv(xx)*xy;
    end
end
