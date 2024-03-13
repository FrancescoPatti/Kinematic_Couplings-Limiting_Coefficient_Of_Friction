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

% external force expressed as a matrix of all contributors. 
% This matrix has as many columns as the applied forces.
FTot= [ 0     %[N]
        0     
        -1];

% Points of application of all external forces.
% Column 'i' is the point of application of force 'i' in matrix FTOT
RForce=1e-3*[0	%[m]
             0	
             100];

%% Geometry definition
n=zeros(3,6);       % normals in the 6 points of contact. each column is a unit vector,
                    % pointing from the vee surface to the ball center
r=zeros(3,6);       % coordinates of the 6 points of contact
Radius=zeros(3,6);  % column "i" is the vector from ball center to the contact point "i"

Radius(:,1)=RotateY(-(pi/2-Beta))*[0; 0; -DBall/2];
Radius(:,2)=RotateY(pi/2-Beta)*[0; 0; -DBall/2];
r(:,1)=[0;Arm;0]+Radius(:,1);    % contact point 1 obtained by vector sum
r(:,2)=[0;Arm;0]+Radius(:,2);    % contact point 2 obtained by vector sum
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
