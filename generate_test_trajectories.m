% Generating random linear task-space trajectories for evaluating the effective mass
% optimization

% Andrew SaLoutos
% 4/16/2021

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



% generates linear trajectories

%% Setup parameters

N = 30; % number of trajectories

% choose to discretize arena space, pick random points, rather than generate random coordinates 
steps = 50;
% arena_x = linspace( 0.4, 1.2, steps);
% arena_y = linspace(-0.4, 0.4, steps);
% min_traj_length = 0.2; % minimum length of trajectory

% for MC and UR3 arms
arena_x = linspace( 0.2, 0.5, steps);
arena_y = linspace(-0.15, 0.15, steps);
min_traj_length = 0.05;



%% Select trajectory endpoints

% trajectories = zeros(4,N); % each column is [x0;xf;y0;yf;]
% % eventually save longer sets of points? or let optimization script run linspace() to fill in trajectories
% 
% % for N trajectories...
% for ii=1:N
%     % uniform random sample for first and second points
%     pt1 = [arena_x(randi(steps));arena_y(randi(steps))];
%     pt2 = [arena_x(randi(steps));arena_y(randi(steps))];
%     while norm((pt2-pt1))<min_traj_length % ensure trajectory is long enough
%         pt2 = [arena_x(randi(steps));arena_y(randi(steps))];
%     end
%     % save data
%     trajectories(:,ii) = [pt1(1);pt2(1);pt1(2);pt2(2)];
% end

%% Plot random trajectories
figure(1); clf; hold on;
for ii=1:N
    plot(trajectories(1:2,ii),trajectories(3:4,ii),'o-','LineWidth',1.5);
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
xlabel('X'); ylabel('Y'); title('30 Randomly Generated Linear Trajectories');
