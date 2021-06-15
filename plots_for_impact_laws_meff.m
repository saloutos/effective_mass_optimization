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
ee_des = [1,0];
q_ik1 = inverse_kinematics_init(p,ee_des,1.1);
q_ik2 = inverse_kinematics_init(p,ee_des,-0.3);

% along a specified finger velocity direction
vf = [ 0.5; 0.5];
evf = vf/norm(vf);

%% get post-collision object velocities
[vo1, vfp1] = object_vels_fixed_base(q_ik1,p,LLo,vf);
[vo2, vfp2] = object_vels_fixed_base(q_ik2,p,LLo,vf);

%% eval exit vel metrics

% evaluate H, J
H1 = A_arm(q_ik1,p);
J1 = jacobian_v_tip(q_ik1,p);

% evaluate OSIMs
LL_inv1 = J1/H1*J1';
LLf1 = inv(LL_inv1);
% evaluate partial OSIM ratio
LLpf1 = (LLf1+LLo)\LLf1;

% evaluate H, J
H2 = A_arm(q_ik2,p);
J2 = jacobian_v_tip(q_ik2,p);

% evaluate OSIMs
LL_inv2 = J2/H2*J2';
LLf2 = inv(LL_inv2);
% evaluate partial OSIM ratio
LLpf2 = (LLf2+LLo)\LLf2;

% check structure
LLf1
[U,S,V] = svd(LLf1)
LLf2
[U,S,V] = svd(LLf2)




%% evaluate along a circle of unit directions for finger velocity?

[V1,D1] = eig(LL_inv1);
% using parametric equations of ellipse
th_ellipse = atan2(V1(2,1),V1(1,1)); % angle between first eigenvector and positive x axis
gamma_ell = 0.4; %/(max([D(1,1),D(2,2)])); % ellipse scaling
l_x = gamma_ell*1/sqrt(D1(1,1)); % TODO: check if sqrt should be here? 
l_y = gamma_ell*1/sqrt(D1(2,2)); 

jj = linspace(0, 2*pi, 100);
u = [cos(jj); sin(jj)];
q_ell1 = zeros(2,100);
for ii=1:100
    q_ell1(:,ii) = gamma_ell*u(:,ii)/(u(:,ii)'*LL_inv1*u(:,ii));
end
q_ell1B = zeros(2,100);
for ii=1:100
    q_ell1B(:,ii) = gamma_ell*u(:,ii)*(u(:,ii)'*LLf1*u(:,ii));
end
% make this belted ellipsoid?
x_ell1 = (l_x*cos(jj))*cos(th_ellipse) - (l_y*sin(jj))*sin(th_ellipse);
y_ell1 = (l_x*cos(jj))*sin(th_ellipse) + (l_y*sin(jj))*cos(th_ellipse);

[V2,D2] = eig(LL_inv2);
% using parametric equations of ellipse
th_ellipse = atan2(V2(2,1),V2(1,1)); % angle between first eigenvector and positive x axis
gamma_ell = 0.4; %/(max([D(1,1),D(2,2)])); % ellipse scaling
l_x = gamma_ell*1/sqrt(D2(1,1)); % TODO: check if sqrt should be here? 
l_y = gamma_ell*1/sqrt(D2(2,2)); 
jj = linspace(0, 2*pi, 100);
u = [cos(jj); sin(jj)];
q_ell2 = zeros(2,100);
for ii=1:100
    q_ell2(:,ii) = gamma_ell*u(:,ii)/(u(:,ii)'*LL_inv2*u(:,ii));
end
q_ell2B = zeros(2,100);
for ii=1:100
    
    q_ell2B(:,ii) = gamma_ell*u(:,ii)*(u(:,ii)'*LLf2*u(:,ii));
end
% make this belted ellipsoid?
x_ell2 = (l_x*cos(jj))*cos(th_ellipse) - (l_y*sin(jj))*sin(th_ellipse);
y_ell2 = (l_x*cos(jj))*sin(th_ellipse) + (l_y*sin(jj))*cos(th_ellipse);


% x_belt = (l_x*cos(jj));
% y_belt = (l_y*sin(jj));
% r_belt = x_belt.^2 + y_belt.^2;
% 
% x_belt2 = (r_belt.*cos(jj))*cos(th_ellipse) - (r_belt.*sin(jj))*sin(th_ellipse);
% y_belt2 = (r_belt.*cos(jj))*sin(th_ellipse) + (r_belt.*sin(jj))*cos(th_ellipse);
% 
% % circle for reference? makes sense with Pvol only
% x_circ = gamma_ell * ( cos(jj)*cos(th_ellipse) - sin(jj)*sin(th_ellipse) );
% y_circ = gamma_ell * ( cos(jj)*sin(th_ellipse) + sin(jj)*cos(th_ellipse) );


%% plot results?
kp1 = keypoints_arm(q_ik1,p);
kp1 = [zeros(2,1), kp1];

kp2 = keypoints_arm(q_ik2,p);
kp2 = [zeros(2,1), kp2];

figure(1); clf; hold on;
% object
obj_plt = plot(kp1(1,end),kp1(2,end),'ok','MarkerSize',15,'MarkerFaceColor','k');
% arm keypoints
arm_plt0 = plot(0,0,'sk','MarkerSize',20,'MarkerFaceColor','k');
arm_plt1A = plot(kp1(1,:),kp1(2,:),'Color',[0, 0, 0.55],'LineWidth',2);
arm_plt1B = plot(kp1(1,:),kp1(2,:),'o','MarkerSize',5,'Color',[0, 0, 0.55],'MarkerFaceColor',[0, 0, 0.55]);
% arm keypoints
arm_plt2A = plot(kp2(1,:),kp2(2,:),'Color',[0.55, 0, 0],'LineWidth',2);
arm_plt2B = plot(kp2(1,:),kp2(2,:),'o','MarkerSize',5,'Color',[0.55, 0, 0],'MarkerFaceColor',[0.55, 0, 0]);
% ellipsoids
meff_plt1 = plot([kp1(1,end)+q_ell1(1,:)],[kp1(2,end)+q_ell1(2,:)],'--b','LineWidth',1.25);
meff_plt2 = plot([kp2(1,end)+q_ell2(1,:)],[kp2(2,end)+q_ell2(2,:)],'--r','LineWidth',1.25);


meff_plt1B = plot([kp1(1,end)+x_ell1],[kp1(2,end)+y_ell1],'-.','Color',[0.68, 0.85, 0.90],'LineWidth',1.25);
meff_plt2B = plot([kp2(1,end)+x_ell2],[kp2(2,end)+y_ell2],'-.','Color',[1.0, 0.8, 0.8],'LineWidth',1.25);


meff_plt1C = plot([kp1(1,end)+q_ell1B(1,:)],[kp1(2,end)+q_ell1B(2,:)],'-*','Color',[0.68, 0.85, 0.90],'LineWidth',1.25);
meff_plt2C = plot([kp2(1,end)+q_ell2B(1,:)],[kp2(2,end)+q_ell2B(2,:)],'-*','Color',[1.0, 0.8, 0.8],'LineWidth',1.25);

% eigenvectors
% meff_vec11 = quiver(kp1(1,end),kp1(2,end),-gamma_ell*V1(1,1)*(1/sqrt(D1(1,1))),-gamma_ell*V1(2,1)*(1/sqrt(D1(1,1))),'b','AutoScale','off','LineWidth',1.5);
% meff_vec11B = quiver(kp1(1,end),kp1(2,end),gamma_ell*V1(1,1)*(1/sqrt(D1(1,1))),gamma_ell*V1(2,1)*(1/sqrt(D1(1,1))),'b','AutoScale','off','LineWidth',1.5);
% 
% meff_vec21 = quiver(kp1(1,end),kp1(2,end),-gamma_ell*V2(1,1)*(1/sqrt(D2(1,1))),-gamma_ell*V2(2,1)*(1/sqrt(D2(1,1))),'r','AutoScale','off','LineWidth',1.5,'MaxHeadSize',0.6);
% meff_vec21B = quiver(kp1(1,end),kp1(2,end),gamma_ell*V2(1,1)*(1/sqrt(D2(1,1))),gamma_ell*V2(2,1)*(1/sqrt(D2(1,1))),'r','AutoScale','off','LineWidth',1.5,'MaxHeadSize',0.6);



% % vf
% vf0 = kp1(:,end)-vf;
% vf_plt = quiver(vf0(1),vf0(2),vf(1),vf(2),1,'r','LineWidth',1.5);
% % vo
% vo_plt = quiver(kp1(1,end),kp1(2,end),vo(1),vo(2),1,'Color','b','LineWidth',1.5);
% % vfp
% vfp_plt = quiver(kp1(1,end),kp1(2,end),vfp(1),vfp(2),1,'Color','g','LineWidth',1.5);

axis equal; 
xlim([-0.1,1.8]); %xlim([-0.5,2.5]); 
ylim([-0.7,0.7]); %ylim([-1.5,1.5]); 
% legend([meff_plt, vf_plt, vo_plt, vfp_plt], {'m_{eff}','vf-','vo+','vf+'});
legend([meff_plt1, meff_plt1B, meff_plt2, meff_plt2B], {'m_{eff} #1','\Lambda_v #1', 'm_{eff} #2', '\Lambda_v #2'});

% title
xlabel('X'); ylabel('Y'); 
str = sprintf('Inertial Ellipsoids at q1 = [%.1f, %.1f, %.1f] and q2 = [%.1f, %.1f, %.1f]', q_ik1(1),q_ik1(2),q_ik1(3),q_ik2(1),q_ik2(2),q_ik2(3));
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


