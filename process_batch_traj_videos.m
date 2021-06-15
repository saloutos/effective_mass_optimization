% process videos of collision trials

% Andrew SaLoutos
% 5/26/2021

% import video frames
% Note: video should be trimmed beforehand to just show trial

% diameter of finger dot is ~31.5
% diameter of mass dot is ~33.8
% roughly 2.45 pixels per mm

filepath = "Hardware Trajectories\Testing Videos\Linear18";

add_out_path = "\Processed";

extension1 = ".mp4";
extension2 = ".avi";

videonames = ["\L18_plain_pos1_trial1";
              "\L18_plain_pos1_trial2";
              "\L18_plain_pos1_trial3";
              "\L18_plain_pos1_trial4";
              "\L18_plain_pos1_trial5"
              "\L18_plain_pos2_trial1";
              "\L18_plain_pos2_trial2";
              "\L18_plain_pos2_trial3";
              "\L18_plain_pos2_trial4";
              "\L18_plain_pos2_trial5";
              "\L18_plain_pos3_trial1";
              "\L18_plain_pos3_trial2";
              "\L18_plain_pos3_trial3";
              "\L18_plain_pos3_trial4";
              "\L18_plain_pos3_trial5";
              "\L18_CA_pos1_trial1";
              "\L18_CA_pos1_trial2";
              "\L18_CA_pos1_trial3";
              "\L18_CA_pos1_trial4";
              "\L18_CA_pos1_trial5";
              "\L18_CA_pos2_trial1";
              "\L18_CA_pos2_trial2";
              "\L18_CA_pos2_trial3";
              "\L18_CA_pos2_trial4";
              "\L18_CA_pos2_trial5";
              "\L18_CA_pos3_trial1";
              "\L18_CA_pos3_trial2";
              "\L18_CA_pos3_trial3";
              "\L18_CA_pos3_trial4";
              "\L18_CA_pos3_trial5"]; 

% ["\S6_plain_pos1_trial1";
% "\S6_plain_pos1_trial2";
% "\S6_plain_pos1_trial3";
% "\S6_plain_pos1_trial4";
% "\S6_plain_pos1_trial5"
% "\S6_plain_pos2_trial1";
% "\S6_plain_pos2_trial2";
% "\S6_plain_pos2_trial3";
% "\S6_plain_pos2_trial4";
% "\S6_plain_pos2_trial5";
% "\S6_plain_pos3_trial1";
% "\S6_plain_pos3_trial2";
% "\S6_plain_pos3_trial3";
% "\S6_plain_pos3_trial4";
% "\S6_plain_pos3_trial5";
% "\S6_CA_pos1_trial1";
% "\S6_CA_pos1_trial2";
% "\S6_CA_pos1_trial3";
% "\S6_CA_pos1_trial4";
% "\S6_CA_pos1_trial5";
% "\S6_CA_pos2_trial1";
% "\S6_CA_pos2_trial2";
% "\S6_CA_pos2_trial3";
% "\S6_CA_pos2_trial4";
% "\S6_CA_pos2_trial5";
% "\S6_CA_pos3_trial1";
% "\S6_CA_pos3_trial2";
% "\S6_CA_pos3_trial3";
% "\S6_CA_pos3_trial4";
% "\S6_CA_pos3_trial5"]; 
          
          
% struct to save data?


video_data = struct('obj_centers',[],'obj_init',[],'obj_final',[],'obj_displacements',[],'time',[],'obj_total_disp',[]);
video_data = repmat(video_data,1,length(videonames));

total_disps = zeros(1,length(videonames));

figure(1);

for jj=1:length(videonames)
    
    fprintf('Processing video #%d\n',jj);
    
    videoname = join([filepath videonames(jj) extension1],'');
    v = VideoReader(videoname);

    nF = v.NumFrames;

%     figure(1);
    obj_centers = zeros(nF,2);
    obj_radii = zeros(nF,1);
    counter = 1;

    % struct to save cropped video with tracking 
    F1(nF) = struct('cdata',[],'colormap',[]);
    
    out_filename = join([filepath add_out_path videonames(jj) extension2],'');
    
    while hasFrame(v)
        vFrame = readFrame(v);
        vFrameCrop = vFrame(:,961:1920,:);

        [centers,radii,metrics] = imfindcircles(vFrameCrop,[16 20],'ObjectPolarity','bright','Sensitivity',0.97);

        % find neon green circle 
        grn_center = [0,0];
        grn_radius = 0;
        if (~isempty(radii))
            center_inds = floor(centers);
            for ii=1:length(radii)
                cc = center_inds(ii,:);
                red_val = vFrameCrop(cc(2),cc(1),1);
                grn_val = vFrameCrop(cc(2),cc(1),2);
                blu_val = vFrameCrop(cc(2),cc(1),3);
                if (red_val<215)&&(grn_val>235)&&(blu_val<205)
                    grn_radius = radii(ii,:);
                    grn_center = centers(ii,:);             
                end
            end

        end

        obj_centers(counter,:) = grn_center;
        obj_radii(counter) = grn_radius;

        % need to show image to get frame for videos
        imshow(vFrameCrop)
        if (counter>1)
            delete(h1);
        end
        h1 = viscircles(grn_center, grn_radius, 'Color','b');
%         pause(1/v.FrameRate);

        % save figure as frame
        
        axis off;
        set( gca, 'Position',[0 0 1 1]);
        
        F1(counter) = getframe(gca);

        counter = counter + 1;
    end
    
    %% save cropped video with tracked circles to new video    
    v1 = VideoWriter(out_filename);
    v1.FrameRate = 30; % optional: control frame rate
    open(v1)
    for ii=1:nF
        writeVideo(v1,F1(ii));
    end
    close(v1)
    
% %     %% preliminary plots
% % %     figure(2); 
% % %     plot(obj_centers(:,2)'/2.45,obj_centers(:,1)'/2.45);
% % %     xlabel('X (mm)'); ylabel('Y (mm)');
% % %     title('Object Center (tracked)');
% % 
% %     time_vec = 0:(1/v.FrameRate):((nF-1)/v.FrameRate);
% % %     figure(3); subplot(3,1,1);
% % %     plot(time_vec, obj_centers(:,2)'/2.45);
% % %     xlabel('Time (s)'); ylabel('X (mm)');
% % %     subplot(3,1,2);
% % %     plot(time_vec, obj_centers(:,1)'/2.45);
% % %     xlabel('Time (s)'); ylabel('Y (mm)');
% % 
% %     x_init = mean(obj_centers(1:15,2)/2.45);
% %     y_init = mean(obj_centers(1:15,1)/2.45);
% %        
% %     displacements = zeros(1,nF);
% %     for ii=1:length(displacements)
% %         dx = obj_centers(ii,2)/2.45-x_init;
% %         dy = obj_centers(ii,1)/2.45-y_init;
% %         displacements(ii) = sqrt(dx^2 + dy^2);
% %     end
% % %     subplot(3,1,3);
% % %     plot(time_vec, displacements);
% % %     xlabel('Time (s)'); ylabel('D (mm)');
% % 
% %     x_final = mean(obj_centers((nF-15):nF,2)/2.45);
% %     y_final = mean(obj_centers((nF-15):nF,1)/2.45);
% %     
% %     total_displacement = sqrt( (x_final-x_init)^2 + (y_final-y_init)^2 );
% %         
% %     total_disps(jj) = total_displacement;
% %     
% %     % save data
% %     video_data(jj).time = time_vec;
% %     video_data(jj).obj_centers = obj_centers/2.45; % save in mm
% %     video_data(jj).obj_init = [y_init, x_init]; % agree with pixel axes
% %     video_data(jj).obj_displacements = displacements; % also saved in mm
% %     video_data(jj).obj_final = [y_final, x_final];
% %     video_data(jj).obj_total_disp = total_displacement;
    

    
    
    
end
 
% total_disps


