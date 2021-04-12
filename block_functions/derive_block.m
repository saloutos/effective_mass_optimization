clear
name = 'cyl_block';

% Define variables for time, generalized coordinates + derivatives, controls, and parameters 
syms t x y th dx dy dth ddx ddy ddth real % coords
syms m I real % inertial params for single block
syms rx ry real % vector from applied force to center of mass (in world frame?)
syms rb real % radius of cylindrical block
syms Fx Fy tau real % forces at center of mass?
syms mu_s g real % coefficient of friction?

% Group them
q   = [    x;    y;   th];      % generalized coordinates
dq  = [   dx;   dy;  dth];      % first time derivatives
ddq = [  ddx;  ddy; ddth];      % second time derivatives
u   = [   Fx;   Fy;  tau];      % controls are really forces
F   = [];                       % forces

% parameters
p  = [m I rb mu_s g]'; % what parameters are necessary here?        

% Generate Vectors and Derivatives
ihat = [1; 0; 0];
jhat = [0; 1; 0];
khat = cross(ihat,jhat);

xhat = [1; 0; 0];
yhat = [0; 1; 0];

e1hat =  cos(th)*ihat + sin(th)*jhat; % x-axis in object frame
e2hat =  -sin(th)*ihat + cos(th)*jhat; % y-axis in object frame

ddt = @(r) jacobian(r,[q;dq])*[dq;ddq]; % a handy anonymous function for taking time derivatives

r_com = [x; y]; % vector to center of mass
dr_com = ddt(r_com);

% r_GP = r_com - [rx; ry]; % vector to point of applied force

% Calculate Kinetic Energy, Potential Energy, and Generalized Forces
F2Q = @(F,r) simplify(jacobian(r,q)'*(F));    % force contributions to generalized forces
M2Q = @(M,w) simplify(jacobian(w,dq)'*(M));   % moment contributions to generalized forces

% kinetic energy
T = (1/2)*m*dot(dr_com,dr_com) + (1/2)*I*dth^2;
% no potential energy
Vg = 0;

% generalized forces (necessary?)
% TODO: add friction here?
% TODO: add force at non-com point here?
Q = u; % already defined above

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
Mass = A;
Grav = simplify(jacobian(Vg, q)');
Corr = simplify( eom + Q - Grav - A*ddq);

% Compute other quantities?

% generalized forces from gripping point? friction?
% impact law?


%% Generate necessary functions to simulate the arm
z  = [q ; dq]; % state variables

% % directory = '../AutoDerived/'; % Write functions to a separate folder because we don't usually have to see them
matlabFunction(A,'file',['A_' name],'vars',{z p});
matlabFunction(b,'file',['b_' name],'vars',{z u p});
matlabFunction(E,'file',['energy_' name],'vars',{z p});
matlabFunction(q,'file',['position_' name],'vars',{z p});
matlabFunction(dq,'file',['velocity_' name],'vars',{z p});
% matlabFunction(J ,'file',['jacobian_tip'],'vars',{z p});
% matlabFunction(Jv ,'file',['jacobian_v_tip'],'vars',{z p});
% matlabFunction(Jw ,'file',['jacobian_w_tip'],'vars',{z p});
% matlabFunction(dJ ,'file',['jacobian_dot_tip'],'vars',{z p});

matlabFunction(Grav,'file',['Grav_' name],'vars',{z p});
matlabFunction(Corr,'file',['Corr_' name],'vars',{z p});
% matlabFunction(fric,'file',['fric_' name],'vars',{z u p});
 

