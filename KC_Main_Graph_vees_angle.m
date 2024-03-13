%% Description
%
% This script calculates the LIMITING COEFFICIENT OF FRICTION of a generic
% Kinematic Coupling.
% The method is explained in detail in:
% F. Patti and J. M. Vogels, "Self-alignment of kinematic couplings: Effects of deformations,"
% Precision Engineering, vol. 60, no. September, pp. 348–354, 2019
% doi: 10.1016/j.precisioneng.2019.08.013.
% The naming convention is as much as possible similar to the one used in
% the paper.
% The output of the script is the graph reported in Figure 6 in the paper.
%
% Inputs:
% the geometry of the coupling is defined by the coordinates of the 6
% points of contact in the matrix "r" and the 6 normals in the matrix "n".
% The load is defined with a number of forces with the matrix "FTot":
% colomn "i" represent the force "i", applied in the point of coordinates
% stored in the column "i" of the matrix RForce.

%%  Overall settings
set(0,'defaulttextinterpreter','latex')
set(0,'defaultaxesticklabelinterpreter','latex')
clear
clc
%% Global variables
global n r FTot RForce DeltaS A NormalReactionsNested mu mu_def;
%% Input
% Standard kc: 3 balls on V-grooves at 120 degress;
% the 3 ball centers are on the xy plane; the 3 lines of the vees direcions intersect at the origin.
% this model is parametric so we can study the effect of geometry variations on the LCF
DBall=20e-3;    %[m] Ball diameter 
Arm=250e-3;     %[m] distance of the ball center from the origin
Beta=45*pi/180; %[rad] half angle of the V 

BetaVect=5:2:65;              %range of the plot
BetaVectRad=BetaVect*pi/180;

% external force expressed as a matrix of all contributors. 
% This matrix has as many columns as the applied forces.
FTot= [ 0     %[N]
        0     
        -100];

% Points of application of all external forces.
% Column 'i' is the point of application of force 'i' in matrix FTOT
RForce=1e-3*[0	%[m]
             0	
             0];

%% Geometry definition
n=zeros(3,6);       % normals in the 6 points of contact. each column is a unit vector,
                    % pointing from the vee surface to the ball center
r=zeros(3,6);       % coordinates of the 6 points of contact
Radius=zeros(3,6);  % column "i" is the vector from ball center to the contact point "i"

% Calculation of the Limiting Coefficient of Friction
MuVect=zeros(length(BetaVect),1);
Mu_defVect=zeros(length(BetaVect),1);

for idx1=1:length(BetaVect)
    Radius(:,1)=RotateY(-(pi/2-BetaVectRad(1,idx1)))*[0; 0; -DBall/2];
    Radius(:,2)=RotateY(pi/2-BetaVectRad(1,idx1))*[0; 0; -DBall/2];
    r(:,1)=[0;Arm;0]+Radius(:,1);    % contact point 1 obtained by vector sum
    r(:,2)=[0;Arm;0]+Radius(:,2);    % contact point 2 obtained by vector sum
    for idx2=3:6                 % contact points by rotating the previous ones 
        Radius(:,idx2)=RotateZ(120*pi/180)*Radius(:,idx2-2);
        r(:,idx2)=RotateZ(120*pi/180)*r(:,idx2-2);
    end
    for idx3=1:6 % the radius vectors are converted into contact normals
        n(:,idx3)=-Radius(:,idx3)/norm(Radius(:,idx3));
    end
    KC_LimitingCoefficientOfFriction
    MuVect(idx1,1)=mu;
    Mu_defVect(idx1,1)=mu_def;
end

figure;
plot(BetaVect,MuVect);
hold on;
plot(BetaVect,Mu_defVect);
plot(BetaVect,BetaVectRad);
ylim([0, 0.4]);
title('Effect of vee angle on LCF');  % Title and labeling to be assigned after the plot
xlabel('$\beta$');
ylabel('LCF');
legend('Without deformations', 'With deformations', 'Self-locking condition');
grid minor;
box on;

%% Functions

function f=RotateX(phi) %rotation matrix around X
f= [1   0         0
    0   cos(phi)  -sin(phi)
    0   sin(phi)  cos(phi)];
end

function f=RotateY(psi) %rotation matrix around Y
f= [cos(psi)    0   sin(psi)
    0           1   0
    -sin(psi)   0   cos(psi)];
end

function f=RotateZ(theta)   %rotation matrix around Z
f= [cos(theta)  -sin(theta) 0
    sin(theta)  cos(theta)  0 
    0           0           1];
end
