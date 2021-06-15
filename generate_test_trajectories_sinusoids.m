% Generating random linear task-space trajectories for evaluating the effective mass
% optimization

% Andrew SaLoutos
% 4/16/2021

% modify to generate sinusoids starting at random points? traj = [x0,y0,a,T], run for a single period?
% add a random rotation angle theta? 


% MC arm
m1 = 0.195;             m2 = 0.262;
m3 = 0.053;             m_motor = 0.527;
I1 = 0.001170;          I2 = 0.001186;
I3 = 0.000096;          I_motor = 0.000508;
Ir = 0.000064;          N = 6;
l_O_m1 = 0.092;         l_A_m2 = 0.201;
l_B_m3 = 0.038;         l_OA = 0.2085;
l_AB = 0.265;           l_BC = 0.1225;
g = 9.81; % do I want gravity to start?
% parameters
p   = [m1 m2 m3 m_motor I1 I2 I3 I_motor Ir N l_O_m1 l_A_m2 l_B_m3 l_OA l_AB l_BC g]'; 



%% Setup parameters

N = 20; % number of trajectories

% choose to discretize arena space, pick random points, rather than generate random coordinates 
steps = 101;

% arena_x = linspace( 0.4, 1.2, steps);
% arena_y = linspace(-0.4, 0.4, steps);
% pt0_x_cands = linspace( 0.4, 0.9, steps); % restrict these so that trajectories end up inside arena?
% pt0_y_cands = linspace(-0.2, 0.2, steps);
% sin_steps = 21;
% amplitudes = linspace(0.1, 0.3, sin_steps);
% periods = linspace(0.1,0.5, sin_steps);
% thetas = linspace(-pi/2, pi/2, steps);

% for MC and UR3 arms
arena_x = linspace( 0.2, 0.5, steps);
arena_y = linspace(-0.15, 0.15, steps);
pt0_x_cands = linspace( 0.2, 0.35, steps); % restrict these so that trajectories end up inside arena?
pt0_y_cands = linspace(-0.1, 0.1, steps);
sin_steps = 21;
amplitudes = linspace(0.025, 0.1, sin_steps);
periods = linspace(0.1,0.2, sin_steps);
thetas = linspace(-pi/4, pi/4, steps);

trajectories = zeros(5,N); % each column is [x0;y0;T;a;th]
num_pts = 51; % number of pts in trajectory
% eventually save longer sets of points? or let optimization script run linspace() to fill in trajectories

%% Select trajectory endpoints

% for N trajectories...
pts = zeros(2,num_pts,N);
num_resamples = 0;
for ii=1:N
    
    valid_traj = 0;
    
    while (valid_traj==0)
    
        % assume a valid trajectory
        valid_traj = 1;
        
        % uniform random sample for first and second points
        pt0 = [arena_x(randi(steps));arena_y(randi(steps))];

        T = periods(randi(sin_steps));
        a = amplitudes(randi(sin_steps));
        th = thetas(randi(steps));

        % get points
        pts_x = linspace(0,T,num_pts);
        pts_y = a*sin((2*pi/T)*pts_x);
        pts_i = [pts_x;pts_y];
        % rotate by theta and translate by pt0
        Rth = [cos(th), -sin(th); sin(th), cos(th)];
        for jj=1:num_pts
            pts_i(:,jj) = pt0 + Rth*pts_i(:,jj);
            % check if new point is outside the arena
            if (pts_i(1,jj)<arena_x(1))||(pts_i(1,jj)>arena_x(end))||(pts_i(2,jj)<arena_y(1))||(pts_i(2,jj)>arena_y(end))
                valid_traj = 0;
                num_resamples = num_resamples + 1;
                continue
            end
        end
        
        
    end
    
    % save in pts array
    pts(:,:,ii) = pts_i;
    
    % save trajectory parameters
    trajectories(:,ii) = [pt0; T; a; th];
end

%% Plot random trajectories
figure(2); clf; hold on;
for ii=1:N
    % just use pts array    
    pts_plt = pts(:,:,ii);
    plot(pts_plt(1,:),pts_plt(2,:),'-','LineWidth',1.5);
end
% plot arena
x_a = [arena_x(1), arena_x(end), arena_x(end), arena_x(1), arena_x(1)];
y_a = [arena_y(1), arena_y(1), arena_y(end), arena_y(end), arena_y(1)];
plot(x_a,y_a,'r--');
% plot arm?
q = [-1.0,2.6,-1.2]';
kp = keypoints_arm(q,p);
rA = kp(:,1);
rB = kp(:,2);
rC = kp(:,3);
plot([0 rA(1)],[0 rA(2)],'Color',[0.3,0.3,0.3],'LineWidth',5);
plot([rA(1) rB(1)],[rA(2) rB(2)],'Color',[0.3,0.3,0.3],'LineWidth',5);
plot([rB(1) rC(1)],[rB(2) rC(2)],'Color',[0.3,0.3,0.3],'LineWidth',5);
plot([0 rA(1) rB(1)],[0 rA(2) rB(2)],'o','Color',[0.3,0.3,0.3],'MarkerSize',15,'MarkerFaceColor',[0.3,0.3,0.3]);

axis equal;
% xlim([arena_x(1)-0.4,arena_x(end)+0.4]); ylim([arena_y(1)-0.4,arena_y(end)+0.4]);
xlim([-0.05,arena_x(end)+0.05]); ylim([arena_y(1)-0.05,arena_y(end)+0.05]);
xlabel('X'); ylabel('Y'); title('20 Randomly Generated Sinusoidal Trajectories');

