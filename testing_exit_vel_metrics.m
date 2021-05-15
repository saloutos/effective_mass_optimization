% testing floating base dynamics

% Andrew SaLoutos
% 5/4/2021

%% add functions

addpath(genpath('arm_floating_functions'))
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

mb = 2; Ib = 0.5*mb*0.1^2;

% arrange in vector
p   = [m1 m2 m3 m_motor I1 I2 I3 I_motor Ir N l_O_m1 l_A_m2 l_B_m3 l_OA l_AB l_BC g mb Ib]';

% define object inertia
mo = 1;
LLo = mo*eye(2);

%% test configuration
xb = 0; yb = 0; thb = 0;
th1 = -0.7;
th2 = 1.4;
th3 = -0.7;

q = [xb; yb; thb; th1; th2; th3];

% along a specified finger velocity direction
vf = [0;-1];

%% eval exit vel metrics
[phi_vf, Pvf, phi_vo, Pvo, phi_vol, Pvol] = eval_exit_vel_metrics(q,p,LLo);


%% evaluate along a circle of unit directions for finger velocity?

[V,D] = eig(Pvol);
% using parametric equations of ellipse
th_ellipse = atan2(V(2,1),V(1,1)); % angle between first eigenvector and positive x axis
gamma_ell = 0.5; %/(max([D(1,1),D(2,2)])); % ellipse scaling
l_x = gamma_ell*D(1,1); 
l_y = gamma_ell*D(2,2); 
jj = linspace(0, 2*pi, 100);
% make this belted ellipsoid?
x_ell = (l_x*cos(jj))*cos(th_ellipse) - (l_y*sin(jj))*sin(th_ellipse);
y_ell = (l_x*cos(jj))*sin(th_ellipse) + (l_y*sin(jj))*cos(th_ellipse);
% circle for reference? doesn't really make sense to add

%% plot results?
kp = keypoints_arm_floating(q,p);
figure(1); clf; hold on;
plot(kp(1,:),kp(2,:),'k','LineWidth',2);
plot(kp(1,:),kp(2,:),'ok','MarkerSize',5,'MarkerFaceColor','k');
plot([kp(1,end)+x_ell],[kp(2,end)+y_ell],'LineWidth',1.25,'Color',[0.85, 0.33, 0.]);
xlim([-0.5,2.0]); ylim([-1.5,1.5]);

%% evaluate across the entire workspace

% set link angles
link1_angles = deg2rad( -90:20:90 );% spaced to avoid singularities
link2_angles = deg2rad( -180:20:180 );
link3_angles = deg2rad( -90:20:90 );

n = length(link1_angles)*length(link2_angles)*length(link3_angles);
config_data = zeros(n,10); % matrix for output data
counter = 1;

% for loops to iterate
for ii = 1:length(link1_angles)
    link1_temp = link1_angles(ii);
    for jj = 1:length(link2_angles)
        link2_temp = link2_angles(jj);
        for kk = 1:length(link3_angles)
            link3_temp = link3_angles(kk);

            q_temp = [0, 0, 0, link1_temp, link2_temp, link3_temp]';

            [phi_vf, Pvf, phi_vo, Pvo, phi_vol, Pvol] = eval_exit_vel_metrics(q_temp,p,LLo);
            
            % evaluate H, J at q_temp
            Hf = H_arm_floating(q_temp,p);
            Jf = J_arm_floating(q_temp,p);

            % evaluate OSIMs
            LL_inv = Jf/Hf*Jf';
            LLf = inv(LL_inv);
            
            % generalized inertia ellipsoid
            gie1f = cond(LLf); 
            gie2f = sqrt(det(LLf)); 
            
            % evaluate fixed base
            H = A_arm(q_temp(4:6),p);
            J = jacobian_v_tip(q_temp(4:6),p);
            LL_inv = J/H*J';
            LL = inv(LL_inv);
            gie1 = cond(LL);
            gie2 = sqrt( det(LL) );
            
            % TODO: figure out what data I want to store at each configuration
            % store data
            config_data(counter,:) = [link1_temp; link2_temp; link3_temp; ...
                                      phi_vf; phi_vo; phi_vol; 
                                      gie1f; gie2f; gie1; gie2]';
            
            disp(counter);
            counter = counter + 1;
            
        end
    end
end

%% plot results over workspace
% figure(2); clf;
% scatter3(config_data(:,1),config_data(:,2),config_data(:,3),20,config_data(:,4),'filled');
% xlabel('q_1'); ylabel('q_2'); zlabel('q_3');
% cb1 = colorbar; 
% title('|\Phi_{vf}|');
% figure(3); clf; 
% scatter3(config_data(:,1),config_data(:,2),config_data(:,3),20,config_data(:,5),'filled');
% xlabel('q_1'); ylabel('q_2'); zlabel('q_3');
% cb2 = colorbar; 
% title('|\Phi_{vo}|');
figure(4); clf;
scatter3(config_data(:,1),config_data(:,2),config_data(:,3),20,config_data(:,6),'filled');
xlabel('q_1'); ylabel('q_2'); zlabel('q_3');
cb3 = colorbar; 
str = sprintf('|\\Phi_{vol}| with m_o=%.1f and m_b=%.1f',mo,mb);
title(str);

