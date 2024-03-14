%% Description
%
% This script calculates the LIMITING COEFFICIENT OF FRICTION of the typical 
% Kinematic Coupling: 3 balls on vees distributed at 120 degrees around an
% axis. The 3 ball centers are on the XY plane, equally distant from the Z axis;
% the 3 lines of the vees directions intersect at the origin.
% The overall geometry is defined by assigning a ball radius, a
% distance between the axis and the balls centers and the angle of the vee.
%
% The method is explained in detail in:
% F. Patti and J. M. Vogels, "Self-alignment of kinematic couplings: Effects of deformations,"
% Precision Engineering, vol. 60, no. September, pp. 348–354, 2019
% doi: 10.1016/j.precisioneng.2019.08.013.
% 
% The naming convention is as much as possible similar to the one used in
% the paper.
%
% Inputs:
% Ball radius: DBall
% Distance from z axis: Arm
% Half angle of the vees: Beta
% 
% The load is defined with a number of forces with the matrix "FTot":
% colomn "i" represent the force "i", applied in the point of coordinates
% stored in the column "i" of the matrix RForce.
%
% Outputs:
% LCF calculated with the classical approach;
% LCF calculated with the approach proposed in the paper (worst case scenario 
% due to deformations); in both cases the corresponding removed 
% constraint is provided;
% Reaction forces when the coupling is in nominal position;
% Reaction forces in the equilibrium configuration that defines the LCF,
% where a certain constraint is removed and friction forces are present.
% This is done in both cases: classical approach and worst case scenario.
%
% Author: Francesco Patti - VDL ETG Eindhoven (NL)%
%
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
DBall=6e-3;    %[m] Ball diameter 
Arm=18.5e-3;     %[m] distance of the ball center from the origin
Beta=45*pi/180; %[rad] half angle of the V 

% external force expressed as a matrix of all contributors. 
% This matrix has as many columns as the applied forces.
FTot= [ 0   0     %[N]
        0   0  
        -330  -6];

% Points of application of all external forces.
% Column 'i' is the point of application of force 'i' in matrix FTOT
RForce=1e-3*[0  0	%[m]
             0	-16.519   
             0  15.348];

%% Geometry definition
n=zeros(3,6);       % normals in the 6 points of contact. each column is a unit vector,
                    % pointing from the vee surface to the ball center
r=zeros(3,6);       % coordinates of the 6 points of contact
Radius=zeros(3,6);  % column "i" is the vector from ball center to the contact point "i"

% The first ball is placed along the y axis, the others are obtained by
% rotating the first one by 120 and 240 degrees.
Radius(:,1)=RotateY(-(pi/2-Beta))*[0; 0; -DBall/2];% definition of radius vector of point 1
Radius(:,2)=RotateY(pi/2-Beta)*[0; 0; -DBall/2];   % definition of radius vector of point 2 
r(:,1)=[0;Arm;0]+Radius(:,1);    % contact point 1 obtained by vector sum of Radius and Arm
r(:,2)=[0;Arm;0]+Radius(:,2);    % contact point 2 obtained by vector sum of Radius and Arm
for idx=3:6                 % contact points by rotating the previous ones 
    Radius(:,idx)=RotateZ(120*pi/180)*Radius(:,idx-2);
    r(:,idx)=RotateZ(120*pi/180)*r(:,idx-2);
end

for idx=1:6 % the radius vectors are converted into contact normals
    n(:,idx)=-Radius(:,idx)/norm(Radius(:,idx));
end

%% Fricton calculations
KC_LimitingCoefficientOfFriction % this solves the equations
% and returns the liming coefficient of friction mu and mu_def
% it also provides the constraint removed, the normal reactions in both
% cases and the normal reactions for the nested position.

%% Displaying the results
Case=["LCF without deformations";"LCF with deformations"];
FrictionCoefficient=[mu ;mu_def];
ConstraintRemoved=[k;k_def];
FrictionTable=table(Case,FrictionCoefficient,ConstraintRemoved);
disp(FrictionTable)
disp(table(NormalReactions,NormalReactions_def))
disp(table(NormalReactionsNested))

%% Exporting to Excel
%FileName='LCF.xlsx';
%writetable(T,FileName);
%winopen(FileName);
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
