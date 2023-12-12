% grid value calculations of cost with alpha and lambda space
% and plot showing the transition of grad of cost to zero
ATOM = "Mg"; lambdaFile = "1_9"; STRUCTURE ="HCP_5-100";

[g,~,magG,KSKS_U,OFOF_U,KSKS,OFOF,nGoU,nGoP,V_UC]...
    = Init2(ATOM,lambdaFile,STRUCTURE);



% varying lambda, alpha, costFunction
lambda_W = linspace(-2,2,100);
alpha_W = linspace(0.1,15,100);


% Initialize the array to store cost and gradients
signNABLA = zeros(length(lambda_W),length(alpha_W));
NABLA = zeros(length(lambda_W), length(alpha_W));
gradLAMBDA = zeros(length(lambda_W), length(alpha_W));
gradALPHA = zeros(length(lambda_W), length(alpha_W));


% Nested loop to calculate cost and gradients
for j = 1:length(alpha_W)
    for i = 1:length(lambda_W)
        [NABLA(i,j),~,~,signNABLA(i,j)]= GRADnabla(lambda_W(i), alpha_W(j));
    end
end

sign_deltacost = signNABLA; % this has +ve as well -ve element
deltacost = NABLA;
%gradalpha = gradALPHA;
%gradlambda = gradLAMBDA;
[Alpha_W,Lambda_W] = meshgrid(alpha_W,lambda_W);


%%%%%%%%%%%%%%%%%%%%%% log10 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure;
% contour(Lambda_W,Alpha_W,deltacost,20); colormap(jet);
% xlabel('\lambda');
% ylabel('\alpha');
% info_text = '$\Delta$';
% annotation('textbox', [0.35, 0.7, 0.10, 0.05],...
%     'String', info_text, 'Interpreter', 'latex', 'BackgroundColor', 'w');
% colorbar;
% figure;
% contour(Lambda_W,Alpha_W,gradalpha,20); colormap(jet);
% xlabel('\lambda');
% ylabel('\alpha');
% info_text = '$\partial \Delta/\partial \alpha$';
% annotation('textbox', [0.35, 0.7, 0.10, 0.05],...
%     'String', info_text, 'Interpreter', 'latex', 'BackgroundColor', 'w');
% colorbar;
% figure;
% contour(Lambda_W,Alpha_W,gradlambda,20); colormap(jet);
% xlabel('\lambda');
% ylabel('\alpha');
% info_text = '$\partial \nabla/\partial \lambda$';
% annotation('textbox', [0.35, 0.7, 0.10, 0.05],...
%     'String', info_text, 'Interpreter', 'latex', 'BackgroundColor', 'w');
% colorbar;
%%%%%%%%%%%%%%%%%%%%%%%%% with log10 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% variable I defined may not be useful- see the write up.
% maglog10sign_deltacost = abs(log10(sign_deltacost));
% phaslog10sign_deltacost = atan(imag(log10...
%     (sign_deltacost))./real(log10(sign_deltacost)));
% 
% log10deltacost = log10(NABLA);
% 
% maglog10gradalpha = abs(log10(gradalpha)); 
% phaslog10gradalpha = atan(imag(log10(gradalpha)./real(log10(gradalpha))));
% 
% maglog10gradlambda = abs(log10(gradlambda));
% phaslog10gradlamdbda = atan(imag(log10...
%     (gradlambda)./real(log10(gradlambda))));


% contour cost
% figure;
% zlev = linspace(-0.0001,0.0001,11);
% contour(Lambda_W,Alpha_W, deltacost,'LevelList',zlev,...
%     'LabelSpacing',100,'ShowText','on');
% colormap(jet);
% xlabel('\lambda');
% ylabel('\alpha');
% info_text = '$\Delta_{sign}$';
% annotation('textbox', [0.35, 0.7, 0.10, 0.05],...
%     'String', info_text, 'Interpreter', 'latex', 'BackgroundColor', 'w');
% colorbar;
% title([ATOM, '-', STRUCTURE]);

%contour Cost_{sign}
% figure;
% contour(Lambda_W,Alpha_W,sign_deltacost,20); colormap(jet);
% xlabel('\lambda');
% ylabel('\alpha');
% info_text = '$\Delta_{sign}$';
% annotation('textbox', [0.35, 0.7, 0.10, 0.05],...
%     'String', info_text, 'Interpreter', 'latex', 'BackgroundColor', 'w');
% colorbar;

%contour log10gradAlpha, log10gradLambda
% figure;
% contour(Lambda_W,Alpha_W,log10gradalpha,20); colormap(jet);
% xlabel('\lambda');
% ylabel('\alpha');
% info_text = '$\partial \Delta/\partial \alpha$';
% annotation('textbox', [0.35, 0.7, 0.10, 0.05],...
%     'String', info_text, 'Interpreter', 'latex', 'BackgroundColor', 'w');
% colorbar;
% figure;
% contour(Lambda_W,Alpha_W,log10gradlambda,20); colormap(jet);
% xlabel('\lambda');
% ylabel('\alpha');
% info_text = '$\partial \nabla/\partial \lambda$';
% annotation('textbox', [0.35, 0.7, 0.10, 0.05],...
%     'String', info_text, 'Interpreter', 'latex', 'BackgroundColor', 'w');
% colorbar;

%surf plot cost Cost_{sign}
figure;
surf(Lambda_W,Alpha_W,deltacost);
xlabel('\lambda');
ylabel('\alpha');
zlabel('\nabla')
% 
% figure;
% surf(Lambda_W,Alpha_W,mod_deltacost);
% xlabel('\lambda');
% ylabel('\alpha');
% zlabel('modified \nabla')

% surf plot log Cost  grad_alpha grad_lambda
% figure;
% surf(Lambda_W,Alpha_W,log10gradalpha);
% xlabel('\lambda');
% ylabel('\alpha');
% zlabel('\partial \nabla/\partial \alpha')
% 
% figure;
% surf(Lambda_W,Alpha_W,log10gradlambda);
% xlabel('\lambda');
% ylabel('\alpha');
% zlabel('\partial \nabla/\partial \lambda')


clearvars -except Alpha_W Lambda_W sign_deltacost deltacost STRUCTURE 


% DIRECTORY = 'contourPlot/Al/VariableAl/';
% fileName = fullfile(DIRECTORY,fileStructure);
% save(fileName,'Alpha_W','Lambda_W','sign_deltacost','deltacost', ...
%     'gradalpha','log10deltacost','maglog10sign_deltacost',...
%     "phaslog10sign_deltacost", 'phaslog10gradalpha','maglog10gradalpha',...
%     'phaslog10gradlamdbda','maglog10gradlambda')

function [NABLA_Val,GRAD_NABLA_LAMBDA,GRAD_NABLA_ALPHA,signNABLA] ...
    = GRADnabla(lambda,alpha)
% GRADnabla - This is a function 
% input - alpha and lambda 
% output - NABLA_Val >> Cost function with NL correction
%       - GRAD_NABLA_LAMBDA >> Gradient of Cost Function with NL wrt lambda
%       - GRAD_NABLA_ALPHA >> Gradient of Cost Function with NL wrt alpha
%       - signNABLA >> this is calculated 
%       a. whether signNABLA^2 = NABLA_Val >> and it was
% This is the test script for doing the average of all the perturbation
% and cost function

% global delKSKS delOFOF delOFOFNL;
% Retriving the value of variable from the base workspace and 
% assign as a local variable
g = evalin("base","g");
magG = evalin("base","magG");
OFOF_U = evalin("base","OFOF_U");
OFOF_U = real(OFOF_U);
KSKS_U = evalin("base","KSKS_U");
KSKS = evalin("base","KSKS");
OFOF = evalin("base","OFOF");
OFOF = real(OFOF);
nGoU = evalin("base","nGoU");
nGoP = evalin("base","nGoP");
V_UC = evalin("base","V_UC");
% note magG_cases is defined global as the reciprocal lattice
% vector remains static in frozen phonon case.

m = size(g,1);
diff = g(m)-g(1);
len = length(g)-1;
del_g = diff/len;

% preallocating variable size
del_val_NL = zeros(length(g),1);
del_val = zeros(length(g),1);
nab_val_NL = zeros(length(g),1);
nab_val = zeros(length(g),1);
del_alpha = zeros(length(g),1);
del_lambda = zeros(length(g),1);
TNL = zeros(length(g),1);
delKSKS = zeros(length(g),1);
delOFOF = zeros(length(g),1);
delOFOFNL = zeros(length(g),1);
gradTNL = zeros(length(g),2);


for k = 1:length(g)
    
    nG = nGoP{k};
    Vuc = V_UC{k};
    [dim1,dim2,dim3] = size(nG);
    pos1 = ceil(dim1/2);
    pos2 = ceil(dim2/2);
    pos3 = ceil(dim3/2);
    
    %% neff is the average density i.e. n_G=0
    %% Kf is the fermi wave-vector
    % density
    neff = nG(pos1,pos2,pos3);
    
    
    
    % optimization parameter alpha defined
    % scaling fermi wave vector kf
    kf = (3*pi^2*neff)^(1./3.);  %>>> scalar
    %kfeta = alpha*(3*pi^2*neff)^(1./3.);
    
    
    precision = 1e-10;
    eta = magG/(2*kf) + precision; %>>> is 3D array
    %eta = magG/(2*kfeta) + precision;
     
    
    %% Lindhard Response for homogeneous electron gas
    % dimension >> 1/([E]*[L^3])
    % dimension >> density/[E]
    %Sus_L = -1/2 * kf/pi^2 * ...
    %    (1 + (1-eta.^2)./(2*eta) .* log(abs((1+eta)./(1-eta))));

    Sus_L_prime = -1/2 * kf/pi^2 * ...
        (1 + (alpha^2-eta.^2)./(2*eta*alpha) .* ...
        log(abs((alpha+eta)./(alpha-eta))));
    
   
    %% Susceptibility for VW limit
    %Sus_W = -kf/pi^2 * 1./(1 + 3* lambda*eta.^2);

    Sus_W_prime = -kf/pi^2 * alpha^2./(alpha^2 + 3* lambda*eta.^2);
    
    
    %% Non local susceptibility
    inv_Sus_NL = -(1./Sus_L_prime - 1./Sus_W_prime);
    
    %% Kernel 
    K_NL =   inv_Sus_NL;
  
    expr_lind1 = 1./(2*eta).*log(abs((alpha+eta)./(alpha-eta)));
    expr_lind2 = alpha./(alpha^2-eta.^2);
    expr_lind_factor = 0.5*kf/pi^2.*(1-eta.^2/alpha^2);
    
    expr_vW1 = (2/alpha - 2*alpha./(alpha^2+3*lambda*eta.^2));
    
    EXPR_LIND_alpha = ...
       1./Sus_L_prime.*expr_lind_factor.*(expr_lind1-expr_lind2);  %eqn15

    EXPR_vW_alpha = -1./Sus_W_prime.*expr_vW1;   % eqn 14
    parK_parAlpha = EXPR_vW_alpha - EXPR_LIND_alpha; %EQN 16

    expr_vW_Lambda_1 = 1./Sus_W_prime;
    expr_vW_Lambda_2 = 3*eta.^2;
    expr_vW_Lamdba_3 = alpha^2+3*lambda*eta.^2;
    expr_vW_Lambda = -expr_vW_Lambda_1.*expr_vW_Lambda_2./expr_vW_Lamdba_3;

    parK_parLambda = expr_vW_Lambda;  % EQN 20
    
    delta_nG = (nG - nGoU);
 
    TNLelement = (delta_nG.^2).* K_NL*Vuc;  %eq-2
    parTNLparAlpha = (delta_nG.^2).* parK_parAlpha * Vuc;  % eq7-only-summand
    parTNLparLambda = (delta_nG.^2).* parK_parLambda * Vuc; %eq-18-only-summand
    
    %% Non-local Energy
    TNL(k) = sum(TNLelement,"all"); %eq -2
   
    %% Gradient Non-local Energy wrt lambda and alpha
    
    gradTNL(k,1) = sum(parTNLparLambda,"all"); %eq-7
    gradTNL(k,2) = sum(parTNLparAlpha,"all"); % eq-18
    
    %% calculating E - E_{min};
    delKSKS(k) = KSKS(k) - KSKS_U;
    delOFOF(k) = OFOF(k) - OFOF_U;
    OFOFNL = OFOF(k) + TNL(k);
    OFOFNL_U = OFOF_U + 0;
    delOFOFNL(k) = OFOFNL - OFOFNL_U;
    
    %% deviation nabla;
   
    del_val_NL(k) = delOFOFNL(k)-delKSKS(k);
    del_val(k) = delOFOF(k) - delKSKS(k);

    nab_val(k) = (delOFOF(k) - delKSKS(k))^2;
    nab_val_NL(k) = (delOFOFNL(k)-delKSKS(k))^2;
    
    %% gradeint calculation of cost function for alpha lambda space
    del_alpha(k) = 2*(delOFOFNL(k)-delKSKS(k))*gradTNL(k,2);
    %eq-6 integrand
    del_lambda(k) = 2*(delOFOFNL(k)-delKSKS(k))*gradTNL(k,1);
    % similar to eq-6 integrand
end

NABLA_Val = trapz(del_g,nab_val_NL);

%NABLA_Val_wo = trapz(del_g,nab_val);
%RatioNABLA = NABLA_Val/NABLA_Val_wo;


GRAD_NABLA_LAMBDA = trapz(del_g,del_lambda);
GRAD_NABLA_ALPHA = trapz(del_g,del_alpha);

% sign dependent nablas calculations
%DELTA = trapz(del_g,del_val);
DELTA_Val = trapz(del_g,del_val_NL); 

                                    
sign_DELTA = sign(DELTA_Val);
signNABLA = sign_DELTA*sqrt(NABLA_Val);

end