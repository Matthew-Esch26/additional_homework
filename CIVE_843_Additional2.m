
% Coordinates 
x1 = 0;  y1 = 0;
x2 = 10; y2 = 20;
x3 = 35; y3 = 20;

% Properties
E = 29000; % ksi 
A1 = 15;   % in^2
A2 = 15;   % in^2
I1 = 1000; % in^4
I2 = 1000; % in^4
w = 2/12;  % kip-in 

% Supports (1 = Fixed, 0 = Pinned)
Support1 = 1;
Support2 = 0;

% Member 1 stiffness
[K1, k1, T1, L1, a] = MemberStiffness(x1, y1, x2, y2, E, A1, I1);

% Member 2 stiffness
[K2, k2, T2, L2, b] = MemberStiffness(x2, y2, x3, y3, E, A2, I2);

Q_f1 = [
    (w*sin(a)*L1)/2;
    (w*cos(a)*L1)/2;
    (w*cos(a)*L1^2)/12;
    (w*sin(a)*L1)/2;
    (w*cos(a)*L1)/2;
    (-w*cos(a)*L1^2)/12
];

Q_f2 = [
    0;
    (7*w*L2)/20 + (3*1.5*w*L2)/20;
    (3*w*L2^2)/60 + (2*1.5*w*L2^2)/60;
    0;
    ((w + 1.5*w)/2)*L2 - ((7*w*L2)/20 + (3*1.5*w*L2)/20);
    (L2/6)*(w*(-2*L2) - 1.5*w*L2) + ((7*w*L2)/20 + (3*1.5*w*L2)/20)*L2 - ((3*w*L2^2)/60 + (2*1.5*w*L2^2)/60)
];

F_f1 = T1' * Q_f1; 
F_f2 = T2' * Q_f2;

if Support1 == 1 && Support2 == 1
    s = [
        K1(4,4) + K2(1,1),   K1(4,5) + K2(1,2),   K1(4,6) + K2(1,3);
        K1(5,4) + K2(2,1),   K1(5,5) + K2(2,2),   K1(5,6) + K2(2,3);
        K1(6,4) + K2(3,1),   K1(6,5) + K2(3,2),   K1(6,6) + K2(3,3)
    ];
    P = [0; 0; 0];
    P_f = [
        F_f1(4,1) + F_f2(1,1);
        F_f1(5,1) + F_f2(2,1);
        F_f1(6,1) + F_f2(3,1)
    ];
    d = s \ (P - P_f);
    
    v1_global = [0;0;0; d(1); d(2); d(3)];
    v2_global = [d(1); d(2); d(3); 0; 0; 0];

elseif Support1 == 1 && Support2 == 0
    s = [
        K1(4,4) + K2(1,1),   K1(5,4) + K2(2,1),   K1(6,4) + K2(3,1),   K2(6,1);
        K1(4,5) + K2(1,2),   K1(5,5) + K2(2,2),   K1(6,5) + K2(3,2),   K2(6,2);
        K1(4,6) + K2(1,3),   K1(5,6) + K2(2,3),   K1(6,6) + K2(3,3),   K2(6,3);
        K2(1,6),             K2(2,6),             K2(3,6),             K2(6,6)
    ];
    P = [0; 0; 0; 0];
    P_f = [
        F_f1(4,1) + F_f2(1,1);
        F_f1(5,1) + F_f2(2,1);
        F_f1(6,1) + F_f2(3,1);
        F_f2(6,1)
    ];
    d = s \ (P - P_f);
    
    v1_global = [0; 0; 0; d(1); d(2); d(3)]; 
    v2_global = [d(1); d(2); d(3); 0; 0; d(4)];
    
else 
    s = [
        K1(3,3),            K1(3,4),                 K1(3,5),                 K1(3,6);
        K1(4,3),            K1(4,4) + K2(1,1),       K1(4,5) + K2(1,2),       K1(4,6) + K2(1,3);
        K1(5,3),            K1(5,4) + K2(2,1),       K1(5,5) + K2(2,2),       K1(5,6) + K2(2,3);
        K1(6,3),            K1(6,4) + K2(3,1),       K1(6,5) + K2(3,2),       K1(6,6) + K2(3,3)
    ];    
    P_f = [
        F_f1(4,1) + F_f2(1,1);
        F_f1(5,1)+F_f2(2,1);
        F_f1(6,1)+F_f2(3,1);
        F_f2(6,1);
    ];

    P = [0;0;0;0]; 
    d = s^-1 * (P-P_f);

        v1_global = [0;0;0; d(1); d(2); d(3)];
        v2_global = [d(1); d(2); d(3); 0; 0; 0];
end

% F Matrix for Member 1
[u1, Q1, F1] = getMemberForces(T1, k1, v1_global, Q_f1);

% F Matrix for Member 2
[u2, Q2, F2] = getMemberForces(T2, k2, v2_global, Q_f2);


if Support1 == 1 && Support2 == 1
    R = [F1(1,1);
    F1(2,1);
    F1(3,1);
    F2(4,1);
    F2(5,1);
    F2(6,1)];  

elseif Support1 == 1 && Support2 == 0
    R = [F1(1,1);
    F1(2,1);
    F1(3,1);
    F2(4,1);
    F2(5,1);];
else 
     R = [F1(1,1);
    F1(2,1);
    F2(4,1);
    F2(5,1);
    F2(6,1)];
end

disp('displacements')
disp(d)
disp('Reactions');
disp(R);


% FUNCTIONS
function [u, Q, F] = getMemberForces(T, k, v_global, Q_f)
    u = T * v_global;          % Local displacements
    Q = k * u + Q_f;           % Local member forces
    F = T' * Q;                % Global member forces
end

function [K, k, T, L, theta] = MemberStiffness(x_b, y_b, x_e, y_e, E, A, I)
    
    % Geometry
    dx = x_e - x_b;
    dy = y_e - y_b;
    L = sqrt(dx^2 + dy^2) * 12; 
    
    theta = atan2(dy, dx); 
    c = cos(theta);
    s = sin(theta);
    
    % Local Stiffness Matrix (k)
    k = [
        E*A/L,   0,           0,           -E*A/L,  0,           0;
        0,       12*E*I/L^3,  6*E*I/L^2,    0,      -12*E*I/L^3, 6*E*I/L^2;
        0,       6*E*I/L^2,   4*E*I/L,      0,      -6*E*I/L^2,  2*E*I/L;
        -E*A/L,  0,           0,            E*A/L,   0,           0;
        0,      -12*E*I/L^3, -6*E*I/L^2,    0,       12*E*I/L^3, -6*E*I/L^2;
        0,       6*E*I/L^2,   2*E*I/L,      0,      -6*E*I/L^2,  4*E*I/L
    ];
    % Transformation Matrix (T)
    T = [
        c,  s,  0,  0,  0,  0;
       -s,  c,  0,  0,  0,  0;
        0,  0,  1,  0,  0,  0;
        0,  0,  0,  c,  s,  0;
        0,  0,  0, -s,  c,  0;
        0,  0,  0,  0,  0,  1
    ];
    % Global Stiffness Matrix (K)
    K = T' * k * T;
end