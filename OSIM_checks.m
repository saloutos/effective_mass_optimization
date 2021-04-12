% Andrew SaLoutos
% 3/30/2021

% script to check effective mass operations

% reference 1: "Safe Robotic Grasping: Minimum Impact-Force Grasp Selection" by Mavrakis et al, 2017
% reference 2: "Inertial Properties in Robotic Manipulation: An Object-Level Framework" by Khatib, 1995

syms m Ixx Ixy Ixz Iyy Iyz Izz real % inertial params
syms rx ry rz real % vector from GP to CoM
syms x y z real % end-effector positions
syms ex ey ez real % euler angles
syms vx vy vz wx wy wz real % linear and angular velocities

% operational space mass matrix at CoM
Ic = [Ixx, Ixy, Ixz; Ixy, Iyy, Iyz; Ixz, Iyz, Izz];
LLoc = [m*eye(3), zeros(3,3); zeros(3,3), Ic];

% translate to grasp point, GP
rhat = [0, -rz, ry; rz, 0, -rx; -ry, rx, 0];
rhat_T = -rhat;
T = [eye(3), rhat; zeros(3,3), eye(3)];
T_T = [eye(3), zeros(3,3); rhat_T, eye(3)];
LLoGP = T_T*LLoc*T; % equation 6 in Mavrakis, can check against equation 59 in Khatib

% change coordinates
% note: inverted E matrices from ETH notes, since they are defined as w = E(x)*xdot
% from Khatib, E matrices are defined as xdot = E(x)*w, which seems more intuitive
Ep = eye(3);
Ep_inv = eye(3);


ER = [1, sin(ex)*sin(ey)/cos(ey), -cos(ex)*sin(ey)/cos(ey); 0, cos(ex), sin(ex); 0, -sin(ex)/cos(ey), cos(ex)/cos(ey)];
ER_T = [1, 0, 0; sin(ex)*sin(ey)/cos(ey), cos(ex), -sin(ex)/cos(ey); -cos(ex)*sin(ey)/cos(ey), sin(ex), cos(ex)/cos(ey)];
ER_inv = [1, 0, sin(ey); 0, cos(ex), -cos(ey)*sin(ex); 0, sin(ex), cos(ex)*cos(ey)];
ER_inv_T = [1, 0, 0; 0, cos(ex), sin(ex); sin(ey), -cos(ey)*sin(ex), cos(ex)*cos(ey)];

E = [Ep, zeros(3,3); zeros(3,3), ER];
E_T = [Ep, zeros(3,3); zeros(3,3), ER_T];
E_inv = [Ep, zeros(3,3); zeros(3,3), ER_inv];
E_inv_T = [Ep, zeros(3,3); zeros(3,3), ER_inv_T];

LLo_Khatib1 = simplify(E_inv_T*LLoGP*E_inv) % equation 61 in Khatib...I'd like to keep things in this form since this is the standard formulation

LLo_Khatib2 = simplify(E_T*LLoGP*E); % this would be correct using the original ETH notes, where E matrices are inverted

LLo_Mavrakis = simplify(E_T*LLoGP*E_inv_T); % equation 7 in Mavrakis, weird combination of both? not sure

% results: 

% upper 3x3 matrix is m*I3 for all of the LLo calculations, so the
% translational effective mass should be exactly the same and equal to the
% object's mass...so the variation in effective mass along the trajectories
% reported in [1] must be due to variation in the arm pose?

% other blocks are different for all of the LLo calculations, I would trust
% the original Khatib derivation, after making sure that E(x) is defined in the
% same way
