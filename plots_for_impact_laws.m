% plots for task space impact laws section of RQE presentation

% Andrew SaLoutos
% 5/11/2021

%% add functions

addpath(genpath('arm_functions'))

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

% arrange in vector
p   = [m1 m2 m3 m_motor I1 I2 I3 I_motor Ir N l_O_m1 l_A_m2 l_B_m3 l_OA l_AB l_BC g]';

% define object inertia
mo = 1.0;
LLo = mo*eye(2);

%% test configuration
a = 1.1;

th1 = a; %0; 
th2 = -2*a; %0.1; 
th3 = a; %0;

q = [th1; th2; th3];

% along a specified finger velocity direction
vf = [ 0.5; 0.5];
evf = vf/norm(vf);

%% get post-collision object velocities
[vo, vfp] = object_vels_fixed_base(q,p,LLo,vf);

%% eval exit vel metrics


% evaluate H, J
H = A_arm(q,p);
J = jacobian_v_tip(q,p);

% evaluate OSIMs
LL_inv = J/H*J';
LLf = inv(LL_inv);
% evaluate partial OSIM ratio
LLpf = (LLf+LLo)\LLf;


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
kp = keypoints_arm(q,p);
kp = [zeros(2,1), kp];

figure(1); clf; hold on;
% object
obj_plt = plot(kp(1,end),kp(2,end),'ok','MarkerSize',15,'MarkerFaceColor','k');
% arm keypoints
arm_plt1 = plot(kp(1,:),kp(2,:),'Color',[0.5, 0.5, 0.5],'LineWidth',2);
arm_plt2 = plot(kp(1,:),kp(2,:),'o','MarkerSize',5,'Color',[0.5, 0.5, 0.5],'MarkerFaceColor',[0.5, 0.5, 0.5]);
% belted ellipsoid
meff_plt = plot([kp(1,end)+x_belt2],[kp(2,end)+y_belt2],'LineWidth',1.25,'Color',[0.85, 0.33, 0.]);
% vf
vf0 = kp(:,end)-vf;
vf_plt = quiver(vf0(1),vf0(2),vf(1),vf(2),1,'r','LineWidth',1.5);
% vo
vo_plt = quiver(kp(1,end),kp(2,end),vo(1),vo(2),1,'Color','b','LineWidth',1.5);
% vfp
vfp_plt = quiver(kp(1,end),kp(2,end),vfp(1),vfp(2),1,'Color','g','LineWidth',1.5);

axis equal; 
xlim([-0.1,1.6]); %xlim([-0.5,2.5]); 
ylim([-0.6,0.6]); %ylim([-1.5,1.5]); 
% legend([meff_plt, vf_plt, vo_plt, vfp_plt], {'m_{eff}','vf-','vo+','vf+'});

% title
xlabel('X'); ylabel('Y'); 
str = sprintf('Collision at q = [%.1f, %.1f, %.1f] and v_f^- = [%.1f, %.1f]', q(1),q(2),q(3),vf(1),vf(2));
title(str);

%% plot LLfp vs scaled LLf and LLo

LLf_scale = logspace(-2,2,101);
LLo_scale = logspace(-2,2,101);

[X,Y] = meshgrid(LLf_scale,LLo_scale);
Z1 = zeros(size(X)); % det of LLpf
Z2 = zeros(size(X)); % magnitude of object exit velocity
Z3 = zeros(size(X)); % magnitude of finger exit velocity
for ii=1:length( LLf_scale )
    llf_temp = LLf_scale(ii);
    for jj=1:length( LLo_scale )
        llo_temp = LLo_scale(jj);    
        LLpf_temp = (llf_temp*LLf + llo_temp*LLo)\(llf_temp*LLf);
        vo_temp = 2*LLpf_temp*vf;
        vfp_temp = (2*LLpf_temp-eye(2))*vf;
        Z1(jj,ii) = det( LLpf_temp );
        Z2(jj,ii) = norm(vo_temp)/norm(vf);
        Z3(jj,ii) = norm(vfp_temp)/norm(vf);     
    end
end

figure; hold on;
surf(X,Y,Z1);
xlabel('\Lambda_f scale'); ylabel('\Lambda_o scale'); zlabel('|\Lambda^p_f|');
set(gca,'xscale','log'); set(gca,'yscale','log');
% colorbar;
% figure;
% surf(X,Y,Z2);
% xlabel('\Lambda_f scale'); ylabel('\Lambda_o scale'); zlabel('Normalized Object Exit Velocity');
% set(gca,'xscale','log'); set(gca,'yscale','log');
% figure;
% surf(X,Y,Z3);
% xlabel('\Lambda_f scale'); ylabel('\Lambda_o scale'); zlabel('Normalized Finger Exit Velocity');
% set(gca,'xscale','log'); set(gca,'yscale','log');

%% slice for LLf

yplt1 = LLo_scale(51)*ones(size(LLo_scale));
yplt2 = LLo_scale(76)*ones(size(LLo_scale));
yplt3 = LLo_scale(26)*ones(size(LLo_scale));

zplt1 = Z1(51,:);
zplt2 = Z1(76,:);
zplt3 = Z1(26,:);

% figure; hold on;
plt1 = plot3(LLf_scale, yplt1, zplt1,'r','LineWidth',2);
plt2 = plot3(LLf_scale, yplt2, zplt2,'r','LineWidth',2);
plt3 = plot3(LLf_scale, yplt3, zplt3,'r','LineWidth',2);

% set(gca,'xscale','log');
% xlabel('\Lambda_f scale'); ylabel('|\Lambda^p_f|');
lg1 = sprintf('\\Lambda_o scale %.2f',LLo_scale(51));
lg2 = sprintf('\\Lambda_o scale %.2f',LLo_scale(71));
lg3 = sprintf('\\Lambda_o scale %.2f',LLo_scale(31));
% legend([plt2,plt1,plt3],{lg2,lg1,lg3});

view(-37.5,30);