%
% Author: Francesco Patti - VDL ETG Eindhoven (NL)%
%
 Input:
% This scripts receives the following quantities in input:
% n:        normals of the 6 points of contact, pointing outward of the
% grooves. Column k is the normal of constraint k.
%
% r:        coordinates of the 6 points of contact. Column k is the
% coordinates of constraint k.
%
% FTot:     external force expressed as a matrix of all contributors. This
% matrix has as many columns as the applied forces,
% 
% RForce:   points of application of all external forces. Column 'i' is the
% point of application of force 'i' in matrix FTOT.

%% Global variables
global n r FTot RForce DeltaS A NormalReactionsNested mu mu_def;

%% Sum of and moment of external load

%FE is column vector, sum of all forces
FE= [0;0;0];
for idx= 1:size(FTot,2)
    FE=FE+FTot(1:3,idx);
end

%Sum of moments w.r.t. origin
SumMom= [0;0;0];
for idx= 1:size(FTot,2)
    SumMom=SumMom+cross(RForce(1:3,idx),FTot(1:3,idx));
end

%% Determining the 6 rigid motions

% This motions are found by solving a linear system of equations in the
% form AMot*X=B and the solution is X=AMot\B

P=zeros(3,6);   % this matrix contains the cross product r_i X n_i, needed to calculate the 6 rigid motions
for idx= 1:6
    P(:,idx)=cross(r(:,idx),n(:,idx));
end

% Building the coefficients matrix AMot
AMot=zeros(6,6);
for i=1:6
    for j=1:3
        AMot(i,j)=n(j,i);
    end
end
for i=1:6
    for j=1:3
        AMot(i,j+3)=P(j,i);
    end
end

%DeltaS{i,k} is the displacement of point of contact "i" when the constraint "k" is removed
DeltaS{6,6}=zeros(3,1); %creating a cell array (array of vectors)
% Building the matrix BMot
% Solving the motion problem
% Obtaining all 6 displacements in points r_i, for each removed constrained k
for k=1:6
    BMot=zeros(6,1);
    BMot(k)=-1;
    X=AMot\BMot;               %Solution of the system of equations
    DeltaS0=X(1:3);         %Displacement of the origin
    DeltaTheta=X(4:6);      %Rotation
    for i=1:6
    DeltaS{i,k}=DeltaS0+cross(DeltaTheta,r(:,i)); %DeltaS{i,k} is the displacement of point of contact "i" when the constraint "k" is removed.
    end
end

%% Definition of the eigenvalue problem

% Building matrix AT

AT=zeros(6,7);
for i=1:3
    for j=1:6
        AT(i,j)=n(i,j);
    end
end

for i=4:6
    VectorKron=[kronecker(i,4);kronecker(i,5);kronecker(i,6)];
    for j=1:6
        AT(i,j)=dot(cross(r(:,j),n(:,j)),VectorKron);
    end
end

AT(1:3,7)=FE(1:3,1);

for i=4:6
    VectorKron=[kronecker(i,4);kronecker(i,5);kronecker(i,6)];
    AT(i,7)=dot(SumMom,VectorKron);
end

% Building Ak: AT with column k removed
A{7}=zeros(6,6);%creating a cell array (array of vectors)
matr=zeros(6,6);
A{1}=AT(:,2:7);  %A1
for k=2:6
    matr(:,1:k-1)=AT(:,1:k-1);  %Ak
    matr(:,k:6)=AT(:,k+1:7);    %Ak
    A{k}=matr;
end
A{7}=AT(1:6,1:6); 
% A7 is not needed for the eigenvalue calcualtion, but it's the matrix of the system
% of equations used to calculate the equilibrium for all 6 constraints active.
% It is used at the end of the calculations.
Omega=zeros(6,1);
Omega1=2*pi*[1;1;1;1;1;1];
%% Finding Equilibrium reactions in nested position
% Without friction, we have 6 unknown normal reactions and 6 equations: 3
% for the forces sum, 3 for the moments sum.
% the problem is expressed in the form AEqX=BEq and solved as AEq\BEq
AEq=A{7};
BEq=zeros(6,1);
BEq(1:3,1)=-FE(1:3,1);
for i=4:6
    VectorKron=[kronecker(i,4);kronecker(i,5);kronecker(i,6)];
    BEq(i,1)=-dot(SumMom,VectorKron);
end
NormalReactionsNested=AEq\BEq;
%% Finding the smallest LCF for all constraints removed

% Friction against incipient motion
% This uses the function FindLCF_K to compare all the 6 friction values and
% provide the smallest of them
mu=inf;
EigenvalueIsFound=0;
for i=1:6
    [EigV, Friction]= FindLCF_K(Omega,i);
    if Friction < mu
    mu=Friction;
    NormalReactions=EigV(1:5)/EigV(6);
    k=i;
    EigenvalueIsFound=1;
    end
end

% Friction directions worst case scenario
% This uses the function FindLCF_K to compare all the 6 friction values and
% provide the smallest of them. In each case the friction directions are
% obtained with the function Omega_min
mu_def=inf;
EigenvalueIsFound=0;
for i=1:6
    [EigV, Friction]= FindLCF_K(Omega_min(i),i);
    if Friction < mu
    mu_def=Friction;
    NormalReactions_def=EigV(1:5)/EigV(6);
    k_def=i;
    EigenvalueIsFound=1;
    end
end

%% Functions definitions

function f= kronecker(m,n) % Kronecker Delta defined
if m==n 
    f=1;
else
    f=0;
end
end

function f= B(Om,k)
%Building the matrix B of the eingenvalue problem
% B is a matrix 6x6, function of the removed constraint k and function of
% omega. Ut is therefore convenient to create a function B(Omega,k).
global r
matr=zeros(6,7); 
f=zeros(6,6);
for i=1:3
	for j=1:6    
		VectorKron=[kronecker(i,1);kronecker(i,2);kronecker(i,3)];
		matr(i,j)= -dot(Xi(j,k)*cos(Om(j))+Eta(j,k)*sin(Om(j)),VectorKron);
	end
end
for i=4:6
	VectorKron=[kronecker(i,4);kronecker(i,5);kronecker(i,6)];
	for j=1:6
		matr(i,j)=-dot(cross(r(:,j),Xi(j,k)*cos(Om(j))+Eta(j,k)*sin(Om(j))),VectorKron);
	end
end
for i=1:6
	matr(i,7)=0;
end
if k==1
    f=matr(:,2:7);  %remove column 1
else
    f(:,1:k-1)=matr(:,1:k-1);%remove column k
    f(:,k:6)=matr(:,k+1:7);
end
end

function f= Xi(i,k)
% Xi is the unit vector of the displacement in point i when constraint k is
% removed
global DeltaS;
f=zeros(3,1);
f=1/norm(DeltaS{i,k})*DeltaS{i,k};
end

function f= Eta(i,k)
% Eta is the unit vector of the direction normal to both displacement and n
% in point i when constraint k is removed.
global DeltaS n;
f=zeros(3,1);
f=cross(1/norm(DeltaS{i,k})*DeltaS{i,k},n(:,i));
end

function [Eigenvector, Friction] = FindLCF_K(Om,k)
% This function solves the generalised problem of eigenvalues, returning the
% LFC for constraint k removed, with friction directions Om.
% It returns the associated eigenvector as well.
global A
EigVectors=zeros(6,6);
mu=zeros(6,6);
[EigVectors, mu]= eig(A{k},-B(Om,k));
% the eigenvalues represent the coefficients of friction for which the equilibrium is verified;
% there is only one physically meaninful value: the smallest real positive. The next routine is selecting this value. 
MinRealEig=inf; % Initial value
EigenvalueIsFound=0; % switch to signal if there is no meaningful value found.
for i=1:6
    if imag(mu(i,i))==0 && real(mu(i,i))>0 && real(mu(i,i)) < MinRealEig
    MinRealEig=mu(i,i);
    Eigenvector=EigVectors(:,i);
    EigenvalueIsFound=1;
    end
end
if EigenvalueIsFound==1
Friction=MinRealEig;
else Friction=NaN;
    Eigenvector=NaN;
end
end

function f = LCF_K(Om,k)
% This function is identical to FindLCF_K, but it returns only the friction
% coefficient. It is needed for finding Omega_min: the function to minimize
% can have only 1 output.
global A
EigVectors=zeros(6,6);
mu=zeros(6,6);
[~, mu]= eig(A{k},-B(Om,k));
% the eigenvalues represent the coefficients of friction for which the equilibrium is verified;
% there is only one physically meaninful value: the smallest real positive. The next routine is selecting this value. 
MinRealEig=inf; % Initial value
EigenvalueIsFound=0; % switch to signal if there is no meaningful value found.
for i=1:6
    if imag(mu(i,i))==0 && real(mu(i,i))>0 && real(mu(i,i)) < MinRealEig
    MinRealEig=mu(i,i);
    EigenvalueIsFound=1;
    end
end
if EigenvalueIsFound==1
f=MinRealEig;
else f=NaN;
end
end

% The following functions are needed for finding the minimum LCF with
% fminsearch: this allow to minimize only for the directions and not the
% removed constraint.
function f=LCF1(Om)
f=LCF_K(Om,1);
end
function f=LCF2(Om)
f=LCF_K(Om,2);
end
function f=LCF3(Om)
f=LCF_K(Om,3);
end
function f=LCF4(Om)
f=LCF_K(Om,4);
end
function f=LCF5(Om)
f=LCF_K(Om,5);
end
function f=LCF6(Om)
f=LCF_K(Om,6);
end

function f=Omega_min(k)
options = optimset('MaxFunEval', 5000,'TolX',1e-4); % setting max iteration number 
if k==1
[minimized_result, ~] = fminsearch(@LCF1, [0;0;0;0;0;0], options);
end
if k==2
[minimized_result, ~] = fminsearch(@LCF2, [0;0;0;0;0;0], options);
end
if k==3
[minimized_result, ~] = fminsearch(@LCF3, [0;0;0;0;0;0], options);
end
if k==4
[minimized_result, ~] = fminsearch(@LCF4, [0;0;0;0;0;0], options);
end
if k==5
[minimized_result, ~] = fminsearch(@LCF5, [0;0;0;0;0;0], options);
end
if k==6
[minimized_result, ~] = fminsearch(@LCF6, [0;0;0;0;0;0], options);
end
f= minimized_result;
end