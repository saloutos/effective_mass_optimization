clear
name = 'arm';

% Define variables for time, generalized coordinates + derivatives, controls, and parameters 
syms t th1 th2 th3 dth1 dth2 dth3 ddth1 ddth2 ddth3 real % coords
syms m1 m2 m3 m_motor I1 I2 I3 I_motor real % inertial params for links and motors
syms l_O_m1 l_A_m2 l_B_m3 g real % CoM values for links
syms l_OA l_AB l_BC real % link lengths, motors are located at points O,A,B respectively
syms tau1 tau2 tau3 Fx Fy real % controls and external forces
syms Ir N real % motor parameters

% Group them
q   = [th1  ; th2 ;  th3];     % generalized coordinates
dq  = [dth1 ; dth2; dth3];     % first time derivatives
ddq = [ddth1;ddth2;ddth3];     % second time derivatives
u   = [tau1 ; tau2; tau3];     % controls
F   = [Fx ; Fy];               % forces

% parameters
p   = [m1 m2 m3 m_motor I1 I2 I3 I_motor Ir N l_O_m1 l_A_m2 l_B_m3 l_OA l_AB l_BC g]';        

% Generate Vectors and Derivatives
ihat = [1; 0; 0];
jhat = [0; 1; 0];
khat = cross(ihat,jhat);

xhat = [1; 0; 0];
yhat = [0; 1; 0];

e1hat =  cos(th1)*ihat + sin(th1)*jhat;
e2hat =  cos(th1+th2)*ihat + sin(th1+th2)*jhat;
e3hat = cos(th1+th2+th3)*ihat + sin(th1+th2+th3)*jhat;

ddt = @(r) jacobian(r,[q;dq])*[dq;ddq]; % a handy anonymous function for taking time derivatives

rA = l_OA*e1hat;
rB = rA + l_AB*e2hat;
rC = rB + l_BC*e3hat;

% r_motor1 = rO; % located at origin
r_motor2 = rA;
r_motor3 = rB;

r_m1 = l_O_m1 * e1hat;
r_m2 = rA + l_A_m2 * e2hat;
r_m3 = rB + l_B_m3 * e3hat;

drA = ddt(rA);
drB = ddt(rB);
drC = ddt(rC);

dr_motor2 = ddt(r_motor2);
dr_motor3 = ddt(r_motor3);

dr_m1 = ddt(r_m1);
dr_m2 = ddt(r_m2);
dr_m3 = ddt(r_m3);

% Calculate Kinetic Energy, Potential Energy, and Generalized Forces
F2Q = @(F,r) simplify(jacobian(r,q)'*(F));    % force contributions to generalized forces
M2Q = @(M,w) simplify(jacobian(w,dq)'*(M));   % moment contributions to generalized forces

omega1 = dth1;
omega2 = dth1 + dth2;
omega3 = dth1 + dth2 + dth3;

% links
T1 = (1/2)*m1*dot(dr_m1,dr_m1) + (1/2)*I1*omega1^2;
T2 = (1/2)*m2*dot(dr_m2,dr_m2) + (1/2)*I2*omega2^2;
T3 = (1/2)*m3*dot(dr_m3,dr_m3) + (1/2)*I3*omega3^2;
% motor bodies
T_m2 = (1/2)*m_motor*dot(dr_motor2,dr_motor2) + (1/2)*I_motor*omega1^2;
T_m3 = (1/2)*m_motor*dot(dr_motor3,dr_motor3) + (1/2)*I_motor*omega2^2;
% rotors
T_mr1 = (1/2)*Ir*(N*th1)^2;
T_mr2 = (1/2)*Ir*(dth1 + N*dth2)^2;
T_mr3 = (1/2)*Ir*(dth1 + dth2 + N*dth3)^2;

% links
Vg1 = m1*g*dot(r_m1, khat);
Vg2 = m2*g*dot(r_m2, khat);
Vg3 = m3*g*dot(r_m3, khat);
% motors
Vg_m2 = m_motor*g*dot(r_motor2, khat);
Vg_m3 = m_motor*g*dot(r_motor3, khat);

T = simplify(T1 + T2 + T3 + T_m2 + T_m3 + T_mr1 + T_mr2 + T_mr3);
Vg = Vg1 + Vg2 + Vg3 + Vg_m2 + Vg_m3;

% generalized forces (necessary?)
Q_tau1 = M2Q(tau1*khat,omega1*khat);
Q_tau2 = M2Q(tau2*khat,omega2*khat); 
Q_tau3 = M2Q(tau3*khat,omega3*khat);

Q_tau2R = M2Q(-tau2*khat,omega1*khat);
Q_tau3R = M2Q(-tau3*khat,omega2*khat);

Q_tau = Q_tau1+(Q_tau2+Q_tau2R)+(Q_tau3+Q_tau3R);

Q_F = F2Q([Fx;Fy;0],rC);

Q = Q_tau;

% Assemble the array of cartesian coordinates of the key points
keypoints = [rA(1:2) rB(1:2) rC(1:2)];

%% All the work is done!  Just turn the crank...
% Derive Energy Function and Equations of Motion
E = T+Vg;
L = T-Vg;
eom = ddt(jacobian(L,dq).') - jacobian(L,q).' - Q;

% Rearrange Equations of Motion
A = jacobian(eom,ddq);
b = A*ddq - eom;

% Equations of motion are
% eom = A *ddq + (coriolis term) + (gravitational term) - Q = 0
Mass_Joint_Sp = A;
Grav_Joint_Sp = simplify(jacobian(Vg, q)');
Corr_Joint_Sp = simplify( eom + Q - Grav_Joint_Sp - A*ddq);

% Compute endpoint jacobian
phi = th1 + th2 + th3;
dphi = dth1 + dth2 + dth3;
qC = [rC(1); rC(2); phi]; % define [x,y,phi] of end-effector as world coordinates
dqC = [drC(1); drC(2); dphi];
J = jacobian(qC,q);

% Partition endpoint Jacobian
Jv = J(1:2,:);
Jw = J(3,:);

% Compute ddt( J )
dJ= reshape( ddt(J(:)) , size(J) );

% Compute operational space eom as well
% From equations 51 and 52 in Khatib paper
% Mass_Joint_Sp_inv = inv(A); % takes a while, but should be okay
% J_pinv = pinv(J);
% J_inv = inv(J);
% Mass_Op_Sp = J_inv' * Mass_Joint_Sp * J_inv;
% Mass_Op_Sp_inv = J * Mass_Joint_Sp_inv * J';
% Mass_Trans_Op_Sp_inv = Jv * Mass_Joint_Sp_inv * Jv';
% 
% 
% 
% Mass_Trans_Op_Sp = inv(Mass_Trans_Op_Sp_inv);
% 
% Mass_Rot_Op_Sp_inv = Jw * Mass_Joint_Sp_inv * Jw';
% 
% Mass_Rot_Op_Sp = inv(Mass_Rot_Op_Sp_inv);
% 
% 
% J_bar = Mass_Joint_Sp_inv * J' * Mass_Op_Sp;
% Corr_Op_Sp = J_bar' * Corr_Joint_Sp - Mass_Op_Sp * dJ * dq;
% Grav_Op_Sp = J_bar' * Grav_Joint_Sp;


%% Generate necessary functions to simulate the arm
z  = [q ; dq]; % state variables

% % directory = '../AutoDerived/'; % Write functions to a separate folder because we don't usually have to see them
% matlabFunction(A,'file',['A_' name],'vars',{z p});
% matlabFunction(b,'file',['b_' name],'vars',{z u p});
% matlabFunction(E,'file',['energy_' name],'vars',{z p});
% matlabFunction(qC,'file',['position_tip'],'vars',{z p});
% matlabFunction(dqC,'file',['velocity_tip'],'vars',{z p});
% matlabFunction(J ,'file',['jacobian_tip'],'vars',{z p});
% matlabFunction(Jv ,'file',['jacobian_v_tip'],'vars',{z p});
% matlabFunction(Jw ,'file',['jacobian_w_tip'],'vars',{z p});
% matlabFunction(dJ ,'file',['jacobian_dot_tip'],'vars',{z p});

% matlabFunction(Grav_Joint_Sp,'file',['Grav_arm_joint'],'vars',{z p});
% matlabFunction(Corr_Joint_Sp,'file',['Corr_arm_joint'],'vars',{z p});
% matlabFunction(keypoints,'file',['keypoints_' name],'vars',{z p});
 
% matlabFunction(Mass_Op_Sp,'file',['Mass_arm_op'],'vars',{z p});
% matlabFunction(Mass_Op_Sp_inv,'file',['LL_arm_op_inv'],'vars',{z p});
% matlabFunction(Mass_Trans_Op_Sp_inv,'file',['LLv_arm_op_inv'],'vars',{z p});
% matlabFunction(Mass_Trans_Op_Sp,'file',['LLv_arm_op'],'vars',{z p});
% matlabFunction(Mass_Rot_Op_Sp_inv,'file',['LLw_arm_op_inv'],'vars',{z p});
% matlabFunction(Mass_Rot_Op_Sp,'file',[LLw_arm_op'],'vars',{z p});
% matlabFunction(Grav_Op_Sp,'file',['Grav_arm_op'],'vars',{z p});
% matlabFunction(Corr_Op_Sp,'file',['Corr_arm_op'],'vars',{z p});

%% Derive cost and constraint functions, gradients and hessians for optimization
% syms e1 e2 real
% 
% e = [e1;e2];
% cv = e' * Mass_Trans_Op_Sp_inv * e;
% grad_cv = jacobian(cv,q); % equivalent to gradient(c,q)

% matlabFunction(cv,'file',['cost_cv'],'vars',{z e p});
% matlabFunction(grad_cv,'file',['grad_cv'],'vars',{z e p});

