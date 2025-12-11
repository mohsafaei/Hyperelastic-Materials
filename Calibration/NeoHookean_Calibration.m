close all;
clear; clc;
% -------------------------------------------------------------
% Mooney–Rivlin parameter identification from test data
%% -------------------------------------------------------------
% This script identifies the two Mooney–Rivlin constants (c1, c2) for an
% incompressible hyperelastic material by minimizing a least–squares
% objective function based on experimental stress–stretch data.
%
% Tests used:
%   - Uniaxial tension
%   - Pure shear
% (Equi–biaxial test is included as a function and can be activated.)
%
% The optimization is performed with fmincon, and the calibrated model is
% compared against the experimental data in combined plots.
%
%
% Author: Mohammad Ali Safaei
% Google Scholar: https://scholar.google.com/citations?user=jD_-4JcAAAAJ&hl
% email: ma.safaei@ut.ac.ir
%
%-------------------------------------------------------------
%% 
resp = input('Include equi-biaxial test? (1 = yes, 0 = no): ');

useEquiBiaxial = (resp == 1);   %boolean
% -------------------------------------------------------------
% Experimental data: uniaxial tension (MPa)
% -------------------------------------------------------------
uniaxial_data = [ ...
    1.00   0.00;
    1.05   0.12;
    1.10   0.21;
    1.15   0.28;
    1.20   0.325;
    1.30   0.37;
    1.40   0.39;
    1.50   0.40;
    1.55   0.41;
    1.60   0.43;
    1.70   0.46;
    1.75   0.48;
    1.80   0.50;
    1.85   0.53;
    1.90   0.56;
    2.00   0.63;
    2.10   0.70;
    2.15   0.74;
    2.20   0.79;
    2.25   0.84;
    2.30   0.89;
    2.35   0.95;
    2.40   1.02;
    2.45   1.10;
    2.50   1.18;
    2.55   1.27;
    2.60   1.37;
    2.65   1.47;
    2.70   1.59;
    2.75   1.71;
    2.85   2.00;
    2.90   2.15];

A1  = uniaxial_data(:,1);
SiggSE_exp  = uniaxial_data(:,2);

%% -------------------------------------------------------------
%% Experimental data: pure shear
%% -------------------------------------------------------------
shear_data = [ ...
    1.000   0.000;
    1.020   0.040;
    1.030   0.070;
    1.040   0.090;
    1.050   0.110;
    1.070   0.160;
    1.080   0.180;
    1.090   0.210;
    1.100   0.230;
    1.110   0.250;
    1.125   0.270;
    1.140   0.290;
    1.155   0.320;
    1.170   0.340;
    1.185   0.360;
    1.200   0.370;
    1.215   0.390;
    1.225   0.410;
    1.240   0.420;
    1.260   0.450;
    1.275   0.460;
    1.300   0.480;
    1.320   0.500;
    1.340   0.515;
    1.360   0.530;
    1.375   0.540;
    1.400   0.550;
    1.420   0.560;
    1.440   0.575;
    1.460   0.590;
    1.480   0.600;
    1.500   0.610];

C1 = shear_data(:,1);
SiggP_exp = shear_data(:,2);

%{
%%A  is  EIGENVALUE stretches in simple Extention
filename1 = 'simpleExtention.xlsx';
A = xlsread(filename1,'D2:E114');
A1 = A(:,1);
A2 = sqrt(1./A1);
SiggS_exp = A(:,2)*10e6;   %Stress experimental data for Simple Extention in Pa unit
%%x is the constant of neo-hookean model
%% 
%%B  is  EIGENVALUE stretches in Equi-Biaxial
filename2 = 'EquiBiaxial.xlsx';
B = xlsread(filename2,'A2:B87');
B1 = B(:,1);
B3 = 1./(B1.^2);
SiggE_exp = B(:,2)*10e6;

%%M_SiggE_exp = mean(SiggE_exp);

 %% ------
%%C  is  EIGENVALUE stretches in Pure Shear----------------
filename3 = 'PureShear.xlsx';
C = xlsread(filename3,'A2:B80');
C1 = C(:,1);
C2 = 1./(C1);
SiggP_exp = C(:,2)*10e6;
%SiggP_exp = (SiggP_exp)*10e6
%}

%% T IS THE OBJECTIVE FUNCTION = OPTIMAZATION Function is defined as RSS(Residual Sum of Squares)
if useEquiBiaxial
    % Uniaxial + pure shear + equi-biaxial
    ObjF = @(x) gge(B1,SiggE_exp,x)+ ggs(A1,SiggSE_exp,x)+ ggp(C1,SiggP_exp,x);
else
    % Uniaxial + pure shear only
    ObjF = @(x) ggs(A1,SiggSE_exp,x)+ ggp(C1,SiggP_exp,x);
end

x0=0;
lb =0.001;
ub =10.0;
[C10,ycte] = fmincon(ObjF,x0,-1,0,[],[],lb,ub);

nnn = 'The constant array of the strain energy function of neo-hookean will be:';
disp(nnn);
disp('xcte');
disp(C10);
bbb = 'The value of the RSS objective function at the end will be(non-important value):';
disp(bbb);
disp('ycte:');
disp(ycte);
%%
figure('Position',[100 100 1200 400]);   % adjust to taste [web:66]

% Create a tight 1x3 layout:
t = tiledlayout(1,3);
t.TileSpacing = 'tight';                 % or 'compact' [web:61][web:62]
t.Padding     = 'tight';                 % or 'compact' [web:61][web:65]
% 

Simple_fun = 2*C10*((A1)-(1./((A1).^2)));   % First Piola Equation
ax1 = nexttile;
Pfit(2) = plot(A1,Simple_fun,'r','LineWidth',3); hold on; grid on;
PData(2) = plot(A1,SiggSE_exp,'Marker','o','MarkerFaceColor', 'b','Color','b');
xlabel('stretch');
ylabel('stress  (in Pa)');
title('Uniaxial');
legend([Pfit(2),PData(2)],'Approximation function','Uniaxial Data','Location','southeast');
ylim([0 1.2*max(SiggSE_exp)]);

if useEquiBiaxial

    ax2 = nexttile;
    Equi_fun= 2*C10*((B1)-(1./(B1.^5)));
    Pfit(2) = plot(B1,Equi_fun,'r','LineWidth',3); hold on; grid on;
    PData(2) = plot(B1,SiggE_exp,'Marker','o','MarkerFaceColor', 'b','Color','b');
    xlabel('stretch');
    ylabel('stress  (in Pa)');
    title('Equi-biaxial');
    legend([Pfit(2),PData(2)],'Approximation function','Equi-biaxial Data','Location','southeast');
    ylim([0 1.2*max(SiggE_exp)]);

else
    nexttile; axis off;
end
    
ax3 = nexttile;

Pure_fun = 2*C10*((C1)-(1./(C1)));
Pfit(3) = plot(C1,Pure_fun,'r','LineWidth',3); hold on; grid on;
PData(3) = plot(C1,SiggP_exp,'Marker','o','MarkerFaceColor', 'b','Color','b');
xlabel('stretch');
ylabel('stress  (in Pa)');
title('Pure shear');
legend([Pfit(3),PData(3)],'Approximation function','Pure Shear Data','Location','southeast');
ylim([0 1.2*max(SiggP_exp)]);

%
%End of the script

% -------------------------------------------------------------
%   Helper functions for Neo-Hookean parameter identification
% -------------------------------------------------------------
% Each function below evaluates the sum of squared errors (SSE) between
% experimental nominal stresses and the Mooney–Rivlin model prediction
% for a specific loading mode. These SSE values are used by the main
% objective function TT in the fmincon optimization.
% -------------------------------------------------------------


%%FOR NEO-HOOKEAN METHOD

function [SiggE_for] = gge(B1,SiggE_exp,x)
FFT=zeros(86,1);
TTE1=zeros(86,1);
Sig_ex = SiggE_exp;
M_Sig_ex = mean(Sig_ex);
numfor=size(B1,1)

for i= 1:numfor
    FFT(i) =(1/(B1(i)))*(2*x*((B1(i)^2)-(1/B1(i)^(4))));
    TTE1(i) =   (Sig_ex(i) - FFT(i))^2;
end 

SiggE_for = sum(TTE1(i));   %%True Stress formulation for Equi-biaxiall in Pa unit


end

 %%FOR NEO-HOOKEAN METHOD

function [SiggS_for] = ggs(A1,SiggS_exp,x)
FFT=zeros(113,1);
TTS1=zeros(113,1);
Sig_ex=SiggS_exp;
M_Sig_ex = mean(Sig_ex);
numfor = size(A1,1)

for i= 1:numfor
    FFT(i) =  (1/(A1(i)))*(2*x*((A1(i)^2)-((1/A1(i)))));
    TTS1(i) = (Sig_ex(i) - FFT(i))^2;
end
SiggS_for = sum(TTS1(i));   %%True Stress formulation for Simple Extention in Pa unit
end


function [SiggP_for] = ggp(C1,SiggP_exp,x)
FFT=zeros(79,1);
TTE1=zeros(79,1);
Sig_ex=SiggP_exp;
M_Sig_ex = mean(Sig_ex);
numfor = size(C1,1)
for i= 1:numfor
    FFT(i) =  (1/(C1(i)))*(2*x*((C1(i)^2)-1));
    TTE1(i) =  (Sig_ex(i) - FFT(i))^2;
end
SiggP_for = sum(TTE1(i));  %%True Stress formulation for Pure Shear in Pa unit
end
