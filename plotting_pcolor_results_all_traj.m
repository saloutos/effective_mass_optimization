%% pseudocolor plots across trials and cost functions

load('avg_data_linear_MC_2.mat');

avg_exit_vel

lin_ex_vel = [avg_exit_vel(1,:); avg_exit_vel(4,:)]
lin_meff = [avg_meff(1,:); avg_meff(4,:)];
lin_pe = [avg_p_error(1,:); avg_p_error(4,:)];
lin_rats = [meff_ratios(3,:); exit_vel_ratios(3,:)]; %solve_times(4,:)./solve_times(1,:)];

% build success table and print
success_table1 = zeros(2,6);
success_table1(2,1) = sum( (meff_ratios(3,:)<1) );
success_table1(2,2) = sum( (exit_vel_ratios(3,:)<1) );
success_table1(1,3:6) = [nanmean(solve_times(1,:),2), nanmean(avg_meff(1,:),2), nanmean(avg_exit_vel(1,:),2), nanmean(avg_p_error(1,:),2)];
success_table1(2,3:6) = [nanmean(solve_times(4,:),2), nanmean(avg_meff(4,:),2), nanmean(avg_exit_vel(4,:),2), nanmean(avg_p_error(4,:),2)];
disp(success_table1)

% get five best trajectories by exit vel ratio
[best_lin_traj,lin_inds] = mink(lin_rats(2,:),5)

load('avg_data_sinusoid_MC_2.mat');
sin_ex_vel = [avg_exit_vel(1,:); avg_exit_vel(4,:)];
sin_meff = [avg_meff(1,:); avg_meff(4,:)];
sin_pe = [avg_p_error(1,:); avg_p_error(4,:)];
sin_rats = [meff_ratios(3,:); exit_vel_ratios(3,:)];

% build success table and print
success_table2 = zeros(2,6);
success_table2(2,1) = sum( (meff_ratios(3,:)<1) );
success_table2(2,2) = sum( (exit_vel_ratios(3,:)<1) );
success_table2(1,3:6) = [nanmean(solve_times(1,:),2), nanmean(avg_meff(1,:),2), nanmean(avg_exit_vel(1,:),2), nanmean(avg_p_error(1,:),2)];
success_table2(2,3:6) = [nanmean(solve_times(4,:),2), nanmean(avg_meff(4,:),2), nanmean(avg_exit_vel(4,:),2), nanmean(avg_p_error(4,:),2)];
disp(success_table2)

% get five best trajectories by exit vel ratio
[best_sin_traj,sin_inds] = mink(sin_rats(2,:),5)

lin_good = [1:30]; %[1:7,9:23,25:30]; % entries that aren't NaN for either row
sin_good = [1:20];

conc_ex_vel = [lin_ex_vel(:,lin_good), sin_ex_vel(:,sin_good)]; % concatenate linear and sinusoidal data
conc_meff = [lin_meff(:,lin_good), sin_meff(:,sin_good)];
conc_pe = [lin_pe(:,lin_good), sin_pe(:,sin_good)];
conc_rats = [lin_rats(:,lin_good), sin_rats(:,sin_good)];


meff_change = (nanmean(conc_meff(2,:))-nanmean(conc_meff(1,:)))/nanmean(conc_meff(1,:)); % changes relative to plain case
ex_vel_change = (nanmean(conc_ex_vel(2,:))-nanmean(conc_ex_vel(1,:)))/nanmean(conc_ex_vel(1,:));
pe_change = (nanmean(conc_pe(2,:))-nanmean(conc_pe(1,:)))/nanmean(conc_pe(1,:));

changes = [meff_change, ex_vel_change, pe_change]

figure(6);
subplot(4,2,1);
sz = size(lin_meff);
C = [lin_meff, zeros(sz(1),1); zeros(1,sz(2)+1)];
pcolor(C);
colorbar;
axis ij;
yticks([1.5,2.5,3.5]);
yticklabels({'Plain','CA'});
xticks(5.5:5:30.5);
xticklabels({'5','10','15','20','25','30'});
title('Effective Mass');

subplot(4,2,3);
sz = size(lin_ex_vel);
C = [lin_ex_vel, zeros(sz(1),1); zeros(1,sz(2)+1)];
pcolor(C);
colorbar;
axis ij;
yticks([1.5,2.5,3.5]);
yticklabels({'Plain','CA'});
xticks(5.5:5:30.5);
xticklabels({'5','10','15','20','25','30'});
title('Exit Velocities');

subplot(4,2,5);
sz = size(lin_pe);
C = [lin_pe, zeros(sz(1),1); zeros(1,sz(2)+1)];
pcolor(C);
colorbar;
axis ij;
yticks([1.5,2.5,3.5]);
yticklabels({'Plain','CA'});
xticks(5.5:5:30.5);
xticklabels({'5','10','15','20','25','30'}); 
title('Position Errors');

subplot(4,2,7);
sz = size(lin_rats);
C = [lin_rats, zeros(sz(1),1); zeros(1,sz(2)+1)];
pcolor(C);
colorbar;
axis ij;
yticks([1.5,2.5,3.5,4.5]);
yticklabels({'m_{eff}','||v_{exit}||'});
xticks(5.5:5:30.5);
xticklabels({'5','10','15','20','25','30'});
title('Ratios, CA / Plain');

subplot(4,2,2);
sz = size(sin_meff);
C = [sin_meff, zeros(sz(1),1); zeros(1,sz(2)+1)];
pcolor(C);
colorbar;
axis ij;
yticks([1.5,2.5,3.5]);
yticklabels({'Plain','CA'});
xticks(5.5:5:30.5);
xticklabels({'5','10','15','20'});
title('Effective Mass');

subplot(4,2,4);
sz = size(sin_ex_vel);
C = [sin_ex_vel, zeros(sz(1),1); zeros(1,sz(2)+1)];
pcolor(C);
colorbar;
axis ij;
yticks([1.5,2.5,3.5]);
yticklabels({'Plain','CA'});
xticks(5.5:5:30.5);
xticklabels({'5','10','15','20'}); 
title('Exit Velocities');

subplot(4,2,6);
sz = size(sin_pe);
C = [sin_pe, zeros(sz(1),1); zeros(1,sz(2)+1)];
pcolor(C);
colorbar;
axis ij;
yticks([1.5,2.5,3.5]);
yticklabels({'Plain','CA'});
xticks(5.5:5:30.5);
xticklabels({'5','10','15','20'}); 
title('Position Errors');

subplot(4,2,8);
sz = size(sin_rats);
C = [sin_rats, zeros(sz(1),1); zeros(1,sz(2)+1)];
pcolor(C);
colorbar;
axis ij;
yticks([1.5,2.5,3.5,4.5]);
yticklabels({'m_{eff}','||v_{exit}||'});
xticks(5.5:5:30.5);
xticklabels({'5','10','15','20'}); 
title('Ratios, CA / Plain');






fs = 14;

figure(7);
subplot(4,1,1);
sz = size(conc_meff);
C = [conc_meff, zeros(sz(1),1); zeros(1,sz(2)+1)];
pcolor(C);
colorbar;
axis ij;
yticks([1.5,2.5]); yticklabels({'Plain','CA'});
xticks(5.5:5:45.5); xticklabels({'5','10','15','20','25','30','35','40','45','50'});
title('Effective Mass');
ax = gca;
ax.FontSize = fs;


subplot(4,1,2);
sz = size(conc_ex_vel);
C = [conc_ex_vel, zeros(sz(1),1); zeros(1,sz(2)+1)];
pcolor(C);
colorbar;
axis ij;
yticks([1.5,2.5]); yticklabels({'Plain','CA'});
xticks(5.5:5:45.5); xticklabels({'5','10','15','20','25','30','35','40','45','50'});
title('Exit Velocities');
ax = gca;
ax.FontSize = fs;

subplot(4,1,3);
sz = size(conc_pe);
C = [conc_pe, zeros(sz(1),1); zeros(1,sz(2)+1)];
pcolor(C);
colorbar;
axis ij;
yticks([1.5,2.5]); yticklabels({'Plain','CA'});
xticks(5.5:5:45.5); xticklabels({'5','10','15','20','25','30','35','40','45','50'});
title('Position Errors');
ax = gca;
ax.FontSize = fs;

subplot(4,1,4);
sz = size(conc_rats);
C = [conc_rats, zeros(sz(1),1); zeros(1,sz(2)+1)];
pcolor(C);
colorbar;
axis ij;
yticks([1.5,2.5]);
yticklabels({'m_{eff}','||v_{exit}||'});
xticks(5.5:5:45.5); xticklabels({'5','10','15','20','25','30','35','40','45','50'});
xlabel('Trajectory #');
title('Ratios, CA / Plain');
ax = gca;
ax.FontSize = fs;

sgtitle('Successful Trajectories, Linear (1-30) and Sinusoidal (31-50)');