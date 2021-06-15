% process videos of collision trials

% Andrew SaLoutos
% 5/26/2021

% import video frames
% Note: video should be trimmed beforehand to just show trial

% diameter of finger dot is ~31.5
% diameter of mass dot is ~33.8
% roughly 2.45 pixels per mm

% neon dot diameters is ~37.5 -> radius is ~18.75 pixels

videoname = 'Hardware Trajectories\Testing Videos\S6_CA_pos3_trial5.mp4';
v = VideoReader(videoname);

% check frames
% disp(v.NumFrames);
% 
% % get first frame
% vFrame = readFrame(v);
% % imshow(vFrame);
% 
% % crop frame
% vFrameCrop = vFrame(:,961:1920,:);
% imshow(vFrameCrop);
% 
% %%
% 
% % draw line to get diameter
% d = drawline;
% % calculate diameter
% pos = d.Position;
% diffPos = diff(pos);
% diameter = hypot(diffPos(1),diffPos(2))
% 
% %%
% 
% % see if circle color is darker or lighter
% grayFrame = rgb2gray(vFrameCrop);
% imshow(grayFrame);
% 
% %%
% 
% [centers,radii,metrics] = imfindcircles(vFrameCrop,[16 20],'ObjectPolarity','bright','Sensitivity',0.97);
% 
% cent_inds = floor(centers)
% radii
% metrics
% 
% for ii=1:length(radii)
%     rgb_val = vFrameCrop(cent_inds(ii,2),cent_inds(ii,1),:);
%     rgb_val = rgb_val(:)
% end
% 
% 
% imshow(vFrameCrop)
% delete(h);
% h = viscircles(centers,radii);

%%
% show frames
figure(1);

obj_centers = zeros(v.NumFrames,2);
obj_radii = zeros(v.NumFrames,1);
counter = 1;

while hasFrame(v)
% for ii=1:3
    vFrame = readFrame(v);
    vFrameCrop = vFrame(:,961:1920,:);
        
%     imshow(vFrame)
    [centers,radii,metrics] = imfindcircles(vFrameCrop,[15 20],'ObjectPolarity','bright','Sensitivity',0.97);
        
    grn_center = [0,0];
    grn_radius = 0;
    if (~isempty(radii))
        % find neon green and yellow circles
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
    
    imshow(vFrameCrop)
    if (counter>1)
%         delete(h1);
        delete(h);
    end
    h1 = viscircles(grn_center, grn_radius, 'Color','b');
%     h = viscircles(centers,radii);
    pause(1/v.FrameRate);
    counter = counter + 1;
end

%% preliminary plots
figure(2); 
plot(obj_centers(:,2)'/2.45,obj_centers(:,1)'/2.45);
xlabel('X (mm)'); ylabel('Y (mm)');
title('Object Center (tracked)');


nF = v.NumFrames;
time_vec = 0:(1/v.FrameRate):((nF-1)/v.FrameRate);

figure(3); subplot(3,1,1);
plot(time_vec, obj_centers(:,2)'/2.45);
xlabel('Time (s)'); ylabel('X (mm)');
subplot(3,1,2);
plot(time_vec, obj_centers(:,1)'/2.45);
xlabel('Time (s)'); ylabel('Y (mm)');

x_init = mean(obj_centers(1:15,2)/2.45);
y_init = mean(obj_centers(1:15,1)/2.45);

displacements = zeros(1,nF);
for ii=1:length(displacements)
    dx = obj_centers(ii,2)/2.45-x_init;
    dy = obj_centers(ii,1)/2.45-y_init;
    displacements(ii) = sqrt(dx^2 + dy^2);
end
subplot(3,1,3);
plot(time_vec, displacements);
xlabel('Time (s)'); ylabel('D (mm)');

x_final = mean(obj_centers((nF-15):nF,2)/2.45);
y_final = mean(obj_centers((nF-15):nF,1)/2.45);

total_displacement = sqrt( (x_final-x_init)^2 + (y_final-y_init)^2 )




