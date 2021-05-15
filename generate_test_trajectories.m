% Generating random linear task-space trajectories for evaluating the effective mass
% optimization

% Andrew SaLoutos
% 4/16/2021

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

trajectories = zeros(4,N); % each column is [x0;xf;y0;yf;]
% eventually save longer sets of points? or let optimization script run linspace() to fill in trajectories

% for N trajectories...
for ii=1:N
    % uniform random sample for first and second points
    pt1 = [arena_x(randi(steps));arena_y(randi(steps))];
    pt2 = [arena_x(randi(steps));arena_y(randi(steps))];
    while norm((pt2-pt1))<min_traj_length % ensure trajectory is long enough
        pt2 = [arena_x(randi(steps));arena_y(randi(steps))];
    end
    % save data
    trajectories(:,ii) = [pt1(1);pt2(1);pt1(2);pt2(2)];
end

%% Plot random trajectories
figure; hold on;
for ii=1:N
    plot(trajectories(1:2,ii),trajectories(3:4,ii),'o-','LineWidth',1.5);
end
% plot arena
x_a = [arena_x(1), arena_x(end), arena_x(end), arena_x(1), arena_x(1)];
y_a = [arena_y(1), arena_y(1), arena_y(end), arena_y(end), arena_y(1)];
plot(x_a,y_a,'r--');
axis equal;
% xlim([arena_x(1)-0.4,arena_x(end)+0.4]); ylim([arena_y(1)-0.4,arena_y(end)+0.4]);
xlim([arena_x(1)-0.1,arena_x(end)+0.1]); ylim([arena_y(1)-0.1,arena_y(end)+0.1]);
xlabel('X'); ylabel('Y'); title('Randomly Generated Linear Trajectories');
