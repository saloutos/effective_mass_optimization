% plotting collision scenarios with exit vel metrics

% Andrew SaLoutos
% 5/5/2021

%% add functions

addpath(genpath('arm_floating_functions'))

%% test parameters
m1 = 1;                 m2 = 1;
m3 = 0.5;               m_motor = 0.5;
I1 = 0.02;              I2 = 0.02;
I3 = 0.01;              I_motor = 0.000625; % assuming radius r=0.05m
Ir = 6.25e-6;           N = 6;
l_O_m1 = 0.25;          l_A_m2 = 0.25;
l_B_m3 = 0.25;          l_OA = 0.5;
l_AB = 0.5;             l_BC = 0.5;
g = 9.81; % do I want gravity to start?

mb = 2; Ib = 0.5*mb*0.1^2;

% arrange in vector
p   = [m1 m2 m3 m_motor I1 I2 I3 I_motor Ir N l_O_m1 l_A_m2 l_B_m3 l_OA l_AB l_BC g mb Ib]';

% define object inertia
mo = 1.0;
LLo = mo*eye(2);

%% test configuration
xb = 0; yb = 0; thb = 0;

a = 1.1;

th1 = a; %0; 
th2 = -2*a; %0.1; 
th3 = a; %0;

q = [xb; yb; thb; th1; th2; th3];

% along a specified finger velocity direction
vf = [ 0.5; 0.0];
evf = vf/norm(vf);

%% get post-collision object velocities
[vo, vol, vfp, vfpl] = object_vels(q,p,LLo,vf);

%% eval exit vel metrics


% evaluate H, J
H = H_arm_floating(q,p);
J = J_arm_floating(q,p);
Hbb = H(1:3,1:3);
Jb = J(:,1:3);

% evaluate OSIMs
LL_inv = J/H*J';
LLf = inv(LL_inv);
LLfl_inv = Jb/Hbb*Jb';
LLfl = inv(LLfl_inv);

% evaluate partial OSIM ratios
LLpf = (LLf+LLo)\LLf;
LLpfl = (LLfl+LLo)\LLfl;



[phi_vf, Pvf, phi_vo, Pvo, phi_vol, Pvol] = eval_exit_vel_metrics(q,p,LLo);


evo = vo/norm(vo);
evol = vol/norm(vol);

eta_vf = evf'*Pvf*evf;
eta_vo = evo'*Pvo*evo;
eta_vol1 = evf'*Pvol*evf;
eta_vol2 = evo'*Pvol*evo;
eta_vol3 = evol'*Pvol*evol;

eta_vol3b = 1 - ((vf'*LLpfl'*LLpf*vf)/(vf'*LLpfl'*LLpfl*vf));

%% evaluate along a circle of unit directions for finger velocity?

[V,D] = eig(LL_inv);
% using parametric equations of ellipse
th_ellipse = atan2(V(2,1),V(1,1)); % angle between first eigenvector and positive x axis
gamma_ell = 0.2; %/(max([D(1,1),D(2,2)])); % ellipse scaling
l_x = gamma_ell*sqrt(D(1,1)); % TODO: check if sqrt should be here? 
l_y = gamma_ell*sqrt(D(2,2)); 
jj = linspace(0, 2*pi, 100);
% make this belted ellipsoid?
x_ell = (l_x*cos(jj))*cos(th_ellipse) - (l_y*sin(jj))*sin(th_ellipse);
y_ell = (l_x*cos(jj))*sin(th_ellipse) + (l_y*sin(jj))*cos(th_ellipse);

x_belt = (l_x*cos(jj));
y_belt = (l_y*sin(jj));
r_belt = x_belt.^2 + y_belt.^2;

x_belt2 = (r_belt.*cos(jj))*cos(th_ellipse) - (r_belt.*sin(jj))*sin(th_ellipse);
y_belt2 = (r_belt.*cos(jj))*sin(th_ellipse) + (r_belt.*sin(jj))*cos(th_ellipse);

% circle for reference? makes sense with Pvol only
x_circ = gamma_ell * ( cos(jj)*cos(th_ellipse) - sin(jj)*sin(th_ellipse) );
y_circ = gamma_ell * ( cos(jj)*sin(th_ellipse) + sin(jj)*cos(th_ellipse) );


%% plot results?
kp = keypoints_arm_floating(q,p);

figure(1); clf; hold on;
% object
obj_plt = plot(kp(1,end),kp(2,end),'ok','MarkerSize',15,'MarkerFaceColor','k');
% arm keypoints
arm_plt1 = plot(kp(1,:),kp(2,:),'Color',[0.5, 0.5, 0.5],'LineWidth',2);
arm_plt2 = plot(kp(1,:),kp(2,:),'o','MarkerSize',5,'Color',[0.5, 0.5, 0.5],'MarkerFaceColor',[0.5, 0.5, 0.5]);
% Pvol ellipse and unit circle
% plot([kp(1,end)+x_ell],[kp(2,end)+y_ell],'LineWidth',1.25,'Color',[0.85, 0.33, 0.]);
% plot([kp(1,end)+x_circ],[kp(2,end)+y_circ],'--','Color',[0.85, 0.33, 0.]);
meff_plt = plot([kp(1,end)+x_belt2],[kp(2,end)+y_belt2],'LineWidth',1.25,'Color',[0.85, 0.33, 0.]);
% vf
vf0 = kp(:,end)-vf;
vf_plt = quiver(vf0(1),vf0(2),vf(1),vf(2),'r','LineWidth',1.5);
% vo, vol
vo_plt = quiver(kp(1,end),kp(2,end),vo(1),vo(2),'Color','b','LineWidth',1.5);
% vol_plt = quiver(kp(1,end),kp(2,end),vol(1),vol(2),'Color','g','LineWidth',1.5);
% vfp, vfpl
vfp_plt = quiver(kp(1,end),kp(2,end),vfp(1),vfp(2),'Color','g','LineWidth',1.5);
% vfpl_plt = quiver(kp(1,end),kp(2,end),vfpl(1),vfpl(2),'Color','m','LineWidth',1.5);

axis equal; 
xlim([-0.5,2.0]); %xlim([-0.5,2.5]); 
ylim([-1.0,1.0]); %ylim([-1.5,1.5]); 
legend([meff_plt, vf_plt, vo_plt, vfp_plt], {'m_{eff}','vf-','vo+','vf+'});
% legend([vf_plt, vo_plt, vfp_plt, vol_plt, vfpl_plt], {'vf-','vo+','vf+','vol+','vfl+'});

% text for metrics
ax_lengths = [sqrt(D(1,1)),sqrt(D(2,2))];
eta_max = max(ax_lengths);
eta_min = min(ax_lengths);
str1 = sprintf('\\eta_{max} : %.3f', eta_max); 
str2 = sprintf('\\eta_{min} : %.3f', eta_min);
str3 = sprintf('\\phi_{vol} : %.3f', phi_vol);
str4 = sprintf('\\eta_{vol} : %.3f', eta_vol3);
% text(-0.3,-1.1,str1);
% text(-0.3,-1.3,str2);
% text(0.5,-1.1,str3);
% text(0.5,-1.3,str4);

% title
xlabel('X'); ylabel('Y'); 
str = sprintf('Collision at q = [%.1f, %.1f, %.1f] and vf = [%.1f, %.1f]', q(4),q(5),q(6),vf(1),vf(2));
title(str);




