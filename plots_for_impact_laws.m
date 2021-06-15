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
% [vop, vfp] = object_vels_fixed_base(q,p,LLo,vf);
[vop, vfp] = object_vels_fixed_base_meff(q,p,mo,vf)

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
gamma_ell = 0.25; %/(max([D(1,1),D(2,2)])); % ellipse scaling
l_x = gamma_ell*1/sqrt(D(1,1)); % TODO: check if sqrt should be here? 
l_y = gamma_ell*1/sqrt(D(2,2)); 
jj = linspace(0, 2*pi, 100);
u = [cos(jj); sin(jj)];
q_ell = zeros(2,100);
for ii=1:100
    q_ell(:,ii) = gamma_ell*u(:,ii)/(u(:,ii)'*LL_inv*u(:,ii));
end

% make this belted ellipsoid?
x_ell = (l_x*cos(jj))*cos(th_ellipse) - (l_y*sin(jj))*sin(th_ellipse);
y_ell = (l_x*cos(jj))*sin(th_ellipse) + (l_y*sin(jj))*cos(th_ellipse);

x_belt = (l_x*cos(jj));
y_belt = (l_y*sin(jj));
r_belt = x_belt.^2 + y_belt.^2;

x_belt2 = (r_belt.*cos(jj))*cos(th_ellipse) - (r_belt.*sin(jj))*sin(th_ellipse);
y_belt2 = (r_belt.*cos(jj))*sin(th_ellipse) + (r_belt.*sin(jj))*cos(th_ellipse);


%% plot results?
kp = keypoints_arm(q,p);
kp = [zeros(2,1), kp];

figure(1); clf; hold on;
% object
obj_plt = plot(kp(1,end),kp(2,end),'ok','MarkerSize',15,'MarkerFaceColor','k');
% arm keypoints
arm_plt0 = plot(0,0,'sk','MarkerSize',20,'MarkerFaceColor','k');
arm_plt1A = plot(kp(1,:),kp(2,:),'Color',[0.5, 0.5, 0.5],'LineWidth',2);
arm_plt1B = plot(kp(1,:),kp(2,:),'o','MarkerSize',5,'Color',[0.5, 0.5, 0.5],'MarkerFaceColor',[0.5, 0.5, 0.5]);
% ellipsoids
% meff_plt1 = plot([kp(1,end)+x_ell],[kp(2,end)+y_ell],'--','LineWidth',1.25,'Color',[0.8500, 0.3250, 0.0980]);
meff_plt2 = plot([kp(1,end)+q_ell(1,:)],[kp(2,end)+q_ell(2,:)],'LineWidth',1.25,'Color',[0.8500, 0.3250, 0.0980]);

% vf
vf0 = kp(:,end)-vf;
vf_plt = quiver(vf0(1),vf0(2),vf(1),vf(2),1,'r','LineWidth',1.5);
% vo
vo_plt = quiver(kp(1,end),kp(2,end),vop(1),vop(2),1,'Color','b','LineWidth',1.5);
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










% %% plot LLfp vs scaled LLf and LLo
% 
% LLf_scale = logspace(-2,2,101);
% LLo_scale = logspace(-2,2,101);
% 
% [X,Y] = meshgrid(LLf_scale,LLo_scale);
% Z1 = zeros(size(X)); % det of LLpf
% Z2 = zeros(size(X)); % magnitude of object exit velocity
% Z3 = zeros(size(X)); % magnitude of finger exit velocity
% for ii=1:length( LLf_scale )
%     llf_temp = LLf_scale(ii);
%     for jj=1:length( LLo_scale )
%         llo_temp = LLo_scale(jj);    
%         LLpf_temp = (llf_temp*LLf + llo_temp*LLo)\(llf_temp*LLf);
%         vo_temp = 2*LLpf_temp*vf;
%         vfp_temp = (2*LLpf_temp-eye(2))*vf;
%         Z1(jj,ii) = det( LLpf_temp );
%         Z2(jj,ii) = norm(vo_temp)/norm(vf);
%         Z3(jj,ii) = norm(vfp_temp)/norm(vf);     
%     end
% end
% 
% figure; hold on;
% surf(X,Y,Z1);
% xlabel('\Lambda_f scale'); ylabel('\Lambda_o scale'); zlabel('|\Lambda^p_f|');
% set(gca,'xscale','log'); set(gca,'yscale','log');
% % colorbar;
% % figure;
% % surf(X,Y,Z2);
% % xlabel('\Lambda_f scale'); ylabel('\Lambda_o scale'); zlabel('Normalized Object Exit Velocity');
% % set(gca,'xscale','log'); set(gca,'yscale','log');
% % figure;
% % surf(X,Y,Z3);
% % xlabel('\Lambda_f scale'); ylabel('\Lambda_o scale'); zlabel('Normalized Finger Exit Velocity');
% % set(gca,'xscale','log'); set(gca,'yscale','log');
% 
% %% slice for LLf
% 
% yplt1 = LLo_scale(51)*ones(size(LLo_scale));
% yplt2 = LLo_scale(76)*ones(size(LLo_scale));
% yplt3 = LLo_scale(26)*ones(size(LLo_scale));
% 
% zplt1 = Z1(51,:);
% zplt2 = Z1(76,:);
% zplt3 = Z1(26,:);
% 
% % figure; hold on;
% plt1 = plot3(LLf_scale, yplt1, zplt1,'r','LineWidth',2);
% plt2 = plot3(LLf_scale, yplt2, zplt2,'r','LineWidth',2);
% plt3 = plot3(LLf_scale, yplt3, zplt3,'r','LineWidth',2);
% 
% % set(gca,'xscale','log');
% % xlabel('\Lambda_f scale'); ylabel('|\Lambda^p_f|');
% lg1 = sprintf('\\Lambda_o scale %.2f',LLo_scale(51));
% lg2 = sprintf('\\Lambda_o scale %.2f',LLo_scale(71));
% lg3 = sprintf('\\Lambda_o scale %.2f',LLo_scale(31));
% % legend([plt2,plt1,plt3],{lg2,lg1,lg3});
% 
% view(-37.5,30);



%% helper function

function q_ik = inverse_kinematics_init(p,ee_i,q1_des)
    % pull out necessary parameters
    l1 = p(14);
    l2 = p(15);
    l3 = p(16);

    % pull out desired end-effector position
    % note: will have to do interpolation eventually here
    x_ee = ee_i(1);
    y_ee = ee_i(2);

    % set q1 initially
    % calculate q2 and q3
    q1_ik = q1_des;
    x_A = l1*cos(q1_ik);
    y_A = l1*sin(q1_ik);

    % alternatively, could set th3 = q1+q2+q3 so that the third link is
    % perpendicular to the desired velocity (to start)
    % could also warm start with a first solve to find the configuration
    % with the third link close to parallel to the desired velocity

    x_ee = x_ee-x_A; % distances taken from second joint
    y_ee = y_ee-y_A;

    q3_ik = acos( (x_ee.^2 + y_ee.^2 - l2^2 - l3^2) / (2 * l2 * l3) ); % this can be +/-, leave as just + for now
    q2_ik = atan2(y_ee,x_ee) - atan2( (l3*sin(q3_ik)), (l2+(l3*cos(q3_ik))) ) - q1_ik;

    % outputs
    q_ik = [q1_ik;q2_ik;q3_ik];

end


