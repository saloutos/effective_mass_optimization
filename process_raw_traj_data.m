function [data,col_time] = process_raw_traj_data(filename)
% process raw data
% return matrix of data w/o "Collision!" text, as well as first time of detected collision
%   Detailed explanation goes here

% data is [t, q1, q2, q3, dq1, dq2, dq3]
raw_data = readmatrix(filename);
col_time = 0;
data = [];
for ii=1:size(raw_data,1)
    if isnan(raw_data(ii,1))
        if col_time==0
            col_time = raw_data((ii-1),1)/1000.0;
        end
    else
        data = [data; raw_data(ii,:)/1000.0];
    end
end



end

