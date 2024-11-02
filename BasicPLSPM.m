function [INI,TABLE,ETC] = BasicPLSPM(z0, W0, B0, modetype,scheme,ind_sign,N_Boot,Max_iter,Min_limit,Flag_Parallel)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BasicPLSPM() - MATLAB function to perform a basic version of Partial    %
%               Least Sqaures Path Modeling  (PLSPM).                     %
% Author: Heungsun Hwang & Gyeongcheol Cho                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments:                                                        %
%   Data = an N by J matrix of scores for N individuals on J indicators   %
%   W = a J by P matrix of weight parameters                              %
%   B = a P by P matrix of path coefficients                              %
%   modetype = a vector of length P indicating the mode of each latent    %
%              variable (1 = mode A, 2 = mode B)                          %
%   scheme = an integer indicating the scheme for updating inner weights  % 
%              (1 = centroid, 2 = factorial, 3 = path weighting)          %
%   ind_sign = a vector of length P representing sign-fixing indicator of %
%              each latent variable                                       %
%   N_Boot = Integer representing the number of bootstrap samples for     %
%            calculating standard errors (SE) and 95% confidence          %
%            intervals (CI)                                               %
%   Max_iter = Maximum number of iterations for the Alternating Least     % 
%              Squares (ALS) algorithm                                    %
%   Min_limit = Tolerance level for ALS algorithm                         %
%   Flag_Parallel = Logical value to determine whether to use parallel    %
%                   computing for bootstrapping                           %
% Output arguments:                                                       %
%   INI: Structure array containing goodness-of-fit values, R-squared     % 
%        values, and matrices parameter estimates                         %
%     .Converge = Logical value indicating whether the ALS algorithm      %
%                 converges within the maximum number of iterations       %
%     .iter = Number of iterations for the ALS algorithm                  %
%     .W: a J by P matrix of weight estimates                             %
%     .C: a P by J matrix of loading estimates                            %
%     .B: a P by P matrix of path coefficient estimates                   %
%     .CVscore: an N by P matrix of component scores                      % 
%  TABLE: Structure array containing tables of parameter estimates, their %
%         SEs, 95% CIs,and other statistics                               %
%     .W: Table for weight estimates                                      %
%     .C: Table for loading estimates                                     %
%     .B: Table for path coefficients estimates                           %
%  ETC: Structure array including bootstrapped parameter estmates         %
%     .W_Boot: Matrix of bootstrapped weight estimates                    %
%     .C_Boot: Matrix of bootstrapped loading estimates                   %
%     .B_Boot: Matrix of bootstrapped path coefficient estimates          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Note:                                                                   %
% 1. This function employs the cost function unifying Mode A and Mode B   %
%    of the PLS algorithm as follows:                                     %
%    "f = alpha_j*SS(Z - FW') + beta_j*SS(F - ZW), where F = Gamma*B"     %
%      alpha_j + beta_j = 1                                               %
%      alpha_j = 1: Mode A                                                %
%       beta_j = 1: Mode B                                                %
%      Z = N by J matrix of indicators                                    %
%      F = N by P matrix of sign-weighted composites of latents           %
%      B = P by P scheme matrix for connected latents                     %
%      A = P by J matrix of weights (~ loadings)                          %
%      W = J by nvlv matrix of weights                                    %
%      Gamma = N by P matrix of latents                                   %
% 2. The detailed description of the algorithm can be found in the        %
%    following paper:                                                     %
%      Hwang, H., Takane, Y. & Tenenhaus, A. An Alternative Estimation    %
%      Procedure for Partial Least Squares Path Modeling. Behaviormetrika %
%      42, 63–78 (2015). https://doi.org/10.2333/bhmk.42.63               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

P = size(W0,2);                    % num of LVs
[N, J] = size(z0);              % num of cases
%% Initialization 
W0=W0~=0; Nw=sum(sum(W0,1),2);
C0=W0';   Nc=sum(sum(C0,1),2);
B0=B0~=0; Nb=sum(sum(B0,1),2);
ind_Bdep=sum(B0,1)>0; Py = sum(ind_Bdep,2);


[est_W,est_C,est_B,it,Converge, est_Gamma] = ALS_BasicPLSPM(z0,W0,B0,modetype,scheme,ind_sign,Max_iter,Min_limit,N,J,P);
INI.iter = it;
INI.Converge=Converge;
INI.W = est_W;
INI.C = est_C;
INI.B = est_B;
INI.CVscore = est_Gamma*sqrt(N);

if N_Boot<100
   TABLE.W=[est_W(W0),NaN(Nw,5)];
   TABLE.C=[est_C(C0),NaN(Nc,5)];
   TABLE.B=[est_B(B0),NaN(Nb,5)];
   ETC.W_Boot=[];
   ETC.C_Boot=[];
   ETC.B_Boot=[];  
else
   W_Boot=zeros(Nw,N_Boot);
   C_Boot=zeros(Nc,N_Boot);
   B_Boot=zeros(Nb,N_Boot);
   if Flag_Parallel
       parfor b=1:N_Boot
           Z_ib = GC_Boot(z0);
           [W_b,C_b,B_b,~,~]=ALS_BasicPLSPM(Z_ib,W0,B0,modetype,scheme,ind_sign,Max_iter,Min_limit,N,J,P);           
           W_Boot(:,b)=W_b(W0);
           C_Boot(:,b)=C_b(C0);
           B_Boot(:,b)=B_b(B0);
       end
   else
       for b=1:N_Boot
           if rem(b,100)==1; fprintf("Bootstrapping %d\n", b); end
           Z_ib = GC_Boot(z0);
           [W_b,C_b,B_b,~,~]=ALS_BasicPLSPM(Z_ib,W0,B0,modetype,scheme,ind_sign,Max_iter,Min_limit,N,J,P);
           W_Boot(:,b)=W_b(W0);
           C_Boot(:,b)=C_b(C0);
           B_Boot(:,b)=B_b(B0);
       end
   end
% (5) Calculation of statistics
   alpha=.05;
   CI=[alpha/2,alpha,1-alpha,1-(alpha/2)];
   loc_CI=round(CI*(N_Boot-1))+1; % .025 .05 .95 .975
% basic statistics for parameter
   TABLE.W=para_stat(est_W(W0),W_Boot,loc_CI);
   TABLE.C=para_stat(est_C(C0),C_Boot,loc_CI); 
   if Py>0; TABLE.B=para_stat(est_B(B0),B_Boot,loc_CI); end
   ETC.W_Boot=W_Boot;
   ETC.C_Boot=C_Boot;
   ETC.B_Boot=B_Boot;    
end
end
function Table=para_stat(est_mt,boot_mt,CI_mp)
   boot_mt=sort(boot_mt,2);
   SE=std(boot_mt,0,2);
   Table=[est_mt,SE,boot_mt(:,CI_mp(1,1)),boot_mt(:,CI_mp(1,4))]; 
end
function in_sample=GC_Boot(Data)
   N=size(Data,1); 
   index=ceil(N*rand(N,1));
   in_sample=Data(index,:); 
end
