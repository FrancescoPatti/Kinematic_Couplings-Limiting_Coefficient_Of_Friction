%% Description
%
% This script calculates the LIMITING COEFFICIENT OF FRICTION of a generic 
% Kinematic Coupling. 
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
% Since any layout can be studied, the geometry is defined by simply assigning the
% coordinates of the 6 points of contact with the matrix "r", where each
% column "i" is the coordinates of point "i". Moreover, the matrix "n"
% contains the normals of the contacts (pointing towards the moving body).
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
%%  Overall settings
set(0,'defaulttextinterpreter','latex')
set(0,'defaultaxesticklabelinterpreter','latex')
clear
clc
%% Global variables
global n r FTot RForce DeltaS A NormalReactionsNested mu mu_def;
%% Input
% Standard kc with vees of 90deg
% normals in the 6 points of contact
n=[ -0.707	0.707	0.354	-0.354	0.354	-0.354
    0	    0	    -0.612	0.612	0.612	-0.612
    0.707	0.707	0.707	0.707	0.707	0.707];

% coordinates of the 6 points of contact
r=[ 7.071	-7.071	-220.042	-212.971	212.971	    220.042
    250	    250	    -118.876	-131.124	-131.124	-118.876
    -7.071	-7.071	-7.071	    -7.071	    -7.071	    -7.071 ];

% external force expressed as a matrix of all contributors. 
% This matrix has as many columns as the applied forces.
FTot= [ 0     
        0     
        -100];

% Points of application of all external forces.
% Column 'i' is the point of application of force 'i' in matrix FTOT
RForce=[0	
        0	
        0];

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