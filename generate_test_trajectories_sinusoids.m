% Generating random linear task-space trajectories for evaluating the effective mass
% optimization

% Andrew SaLoutos
% 4/16/2021

% modify to generate sinusoids starting at random points? traj = [x0,y0,a,T], run for a single period?
% add a random rotation angle theta? 

%% Setup parameters

N = 30; % number of trajectories

% choose to discretize arena space, pick random points, rather than generate random coordinates 
steps = 101;
arena_x = linspace( 0.4, 1.2, steps);
arena_y = linspace(-0.4, 0.4, steps);

pt0_x_cands = linspace( 0.4, 0.9, steps); % restrict these so that trajectories end up inside arena?
pt0_y_cands = linspace(-0.2, 0.2, steps);

sin_steps = 21;
amplitudes = linspace(0.1, 0.3, sin_steps);
periods = linspace(0.1,0.5, sin_steps);
thetas = linspace(-pi/2, pi/2, steps);

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
figure; hold on;
for ii=1:N
    % just use pts array    
    pts_plt = pts(:,:,ii);
    plot(pts_plt(1,:),pts_plt(2,:),'-','LineWidth',1.5);
end
% plot arena
x_a = [0.4, 1.2, 1.2, 0.4, 0.4];
y_a = [-0.4, -0.4, 0.4, 0.4, -0.4];
plot(x_a,y_a,'r--');
axis equal;
xlim([arena_x(1)-0.4,arena_x(end)+0.4]); ylim([arena_y(1)-0.4,arena_y(end)+0.4]);
xlabel('X'); ylabel('Y'); title('Randomly Generated Sinusoidal Trajectories');

