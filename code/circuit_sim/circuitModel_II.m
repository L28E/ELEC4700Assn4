function [G,C] = circuitModel(Cn)
%circuitModel Builds the matrices for modelling the circuit

    R1 = 1;
    R2 = 2;
    R3 = 4.702629;    
    R4 = 0.1;
    Ro = 1000;
    C1 = 0.25;       
    L = 0.2;
    a = 100;

    %% The G-matrix represents the voltage relationships bewteen the linear components of the circuit;
    % It gets filled with the coeffcients of the set of equations 
    
    % The C-matrix represents the 1st order voltage relationships between the
    % circuit components (I think)

    % For each equation, let's write out the coeffcients in the following order: V1, V2, V3, V4, VO, IL
    G=zeros(6);
    C=zeros(6);

    % Node 1:
    G(1,1)=1;

    % Node 2:
    G(2,1)=1/R1;
    G(2,2)=-(1/R1+1/R2);
    G(2,6)=-1;
    C(2,1)=C1;
    C(2,2)=-C1;

    % Node 3:
    % I'll include the capacitor, and then later Ill put the noisy current guy in the forcing vector 
    G(3,3)=-1/R3;
    G(3,6)=1;
    C(3,3)=-Cn;    

    % Node 4:
    G(4,3)=-a*1/R3;
    G(4,4)=1;
    C(4,3)=-a*Cn;

    % Node 5:
    G(5,4)=1/R4;
    G(5,5)=-(1/R4+1/Ro);

    % The inductor:
    G(6,2)=1;
    G(6,3)=-1;
    C(6,6)=-L;

end

