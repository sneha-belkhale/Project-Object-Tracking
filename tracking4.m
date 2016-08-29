
%now, we do the kalman filter multiple object tracking
%clc; % Clear command window.
%clear all; % Get rid of variables from prior run of this m-file.
%disp('tracking');
%workspace; % Show panel with all the variables.

%set(0,'DefaultFigureWindowStyle','docked') %dock the figures..just a personal preference you don't need this.
%base_dir ='C:\mandy\ucla\research\imaging\Tracking-Math\M\region1';

%cd(base_dir);
%load the tracking
load('vesicles_detections.mat')

%get frame list
%f_list =  dir('*jp2');

%% set our threshold parameters

MUNKRES_REJECT=20; 
COLLISION_REJECT=0; 
REMOVE_FRAME=5; 

%% define main variables for KALMAN FILTER!
dt = 1;  %our sampling rate
S_frame = 1; %starting frame 
E_frame = 14;%ending frame

%u = 0; % define acceleration magnitude to start
processNoise = 0.1; %process noise: the variability in how fast the vesicle is speeding up       % usually between 0.05-0.5
tkn_x = .1;  %measurement noise in the horizontal direction (x axis).
tkn_y = .1;  %measurement noise in the horizontal direction (y axis).
Ez = [tkn_x 0; 0 tkn_y];
Ex = [dt^4/4 0 dt^3/2 0; ...
    0 dt^4/4 0 dt^3/2; ...
    dt^3/2 0 dt^2 0; ...
    0 dt^3/2 0 dt^2].*processNoise^2; % Ex convert the process noise (stdv) into covariance matrix
P = Ex; % estimate of initial position variance (covariance matrix)


%% Define update equations in 2-D! (Coefficent matrices): A physics based model for where we expect the HEXBUG to be [state transition (state + velocity)] + [input control (acceleration)]
A = [1 0 dt 0; 0 1 0 dt; 0 0 1 0; 0 0 0 1]; %state update matrice
B = [(dt^2/2); (dt^2/2); dt; dt];
C = [1 0 0 0; 0 1 0 0];  %this is our measurement function C, that we apply to the state estimate Q to get our expect next/new measurement

%% initize result variables
Q_loc_meas = []; % the fly detecions  extracted by the detection algo

%% initialize estimation variables for two dimensions
Q= [X{S_frame} Y{S_frame} zeros(length(X{S_frame}),1) zeros(length(X{S_frame}),1)]'

Q_estimate = nan(4,100);
Q_estimate(:,1:size(Q,2)) = Q;  %estimate of initial location estimation of where the flies are(what we are updating)
Q_loc_estimateY = nan(100); %  position estimate
Q_loc_estimateX= nan(100); %  position estimate
P_estimate = P;  %covariance estimator
strk_trks = zeros(1,100);  %counter of how many strikes a track has gotten
nD = size(X{S_frame},1); %initize number of detections
nF =  find(isnan(Q_estimate(1,:))==1,1)-1;  %initize number of track estimates
u=zeros(nF,2);

%for each frame
for t = S_frame:E_frame
    
    % load the image
    img_tmp = double(imread(f_list(t).name));
    img = img_tmp(:,:,1);
    %get size
    [imsize_x, imsize_y,z]=size(img);
    % make the given detections matrix
    Q_loc_meas = [X{t} Y{t}];
    
    %% do the kalman filter
    % Predict next state of the flies with the last state and predicted motion.
    nD = size(X{t},1); %set new number of detections
   
    %% now we assign the detections to estimated track positions
    %make the distance (cost) matrice between all pairs rows = tracks, coln =
    %detections
    est_dist = pdist([Q_estimate(1:2,1:nF)'; Q_loc_meas]);
    est_dist = squareform(est_dist); %make square
    est_dist = est_dist(1:nF,nF+1:end) ; %limit to just the tracks to detection distances
    
  
    [asgn, cost] =munkres(est_dist); %do the assignment with munkres algorithm
    
    
    asgn = asgn';
    
    % check some special situations
    
    %check 1: is the detection far from the observation? if so, reject it.
    rej = zeros(size(asgn));  
    %for F = 1:size(asgn')
    for F = 1:nF
        if asgn(F) > 0
            rej(F) =  est_dist(F,asgn(F)) <MUNKRES_REJECT ; 
            %rej(F) =  est_dist(F,asgn(F)) <8 ; %5-10
        else
            rej(F) = 0;
        end
    end
    asgn = asgn.*rej;
    
    %For the estimates that did not get assigned, account for collisions by
    %checking if the new track is close enough to a previous estimate
    no_trk_list =  find(asgn==0);
    if ~isempty(no_trk_list)
    %find the nearest neighbors of the unassigned tracks 
    asgn_dup=knnsearch(Q_loc_meas,Q_estimate(1:2,no_trk_list)');
    %check 1: is the detection far from the observation? if so, reject it.
    rej = zeros(size(asgn_dup));  
    %for F = 1:size(asgn')
    for F = 1:size(asgn_dup)
        if asgn_dup(F) > 0
            rej(F) =  est_dist(no_trk_list(F),asgn_dup(F)) <COLLISION_REJECT ; 
            %rej(F) =  est_dist(F,asgn(F)) <8 ; %5-10
        else
            rej(F) = 0;
        end
    end
    asgn_dup = asgn_dup.*rej;
    %update our assignments 
    asgn(no_trk_list)=asgn_dup;
    end   
    
    %apply the assignment to the update
    k = 1;
    for F = 1:length(asgn)
        if asgn(F) > 0
            %Q_estimate(:,k) = Q_estimate(:,k) + K * (Q_loc_meas(asgn(F),:)' - C * Q_estimate(:,k));
            Q_estimate(1:2,k) = Q_loc_meas(asgn(F),:)';
        end
        k = k + 1;
    end

    %% Store data
    Q_loc_estimateX(t,1:nF) = Q_estimate(1,1:nF);
    Q_loc_estimateY(t,1:nF) = Q_estimate(2,1:nF);
    
    %ok, now that we have our assignments and updates, lets find the new detections and
    %lost trackings
    
    %find the new detections. basically, anything that doesn't get assigned
    %is a new tracking
    new_trk = Q_loc_meas(~ismember(1:size(Q_loc_meas,1),asgn),:)';
    
   %add new tracks to the end 
    if ~isempty(new_trk)
        Q_estimate(:,nF+1:nF+size(new_trk,2))=  [new_trk; zeros(2,size(new_trk,2))];
        nF = nF + size(new_trk,2);  % number of track estimates with new ones included
    end
    
   
    %give a strike to any tracking that didn't get matched up to a
    %detection
    no_trk_list =  find(asgn==0);
    
    if ~isempty(no_trk_list)
        strk_trks(no_trk_list) = strk_trks(no_trk_list) + 1;
        u(no_trk_list,:)=0;
    end
    
    %if a track has a strike greater than 30, delete the tracking. i.e.
    bad_trks = find(strk_trks >REMOVE_FRAME);                    %% 10-30  (40frames)
    Q_estimate(:,bad_trks) = NaN;
    %
    clf
    img = imread(f_list(t).name);
    img=img(:,:,1); %for jp3
    img=im2double(img); %for jp3    
    imshow(img);
    
    hold on;
    plot(X{t}(:),Y{t}(:),'or'); % the actual tracking
    T = size(Q_loc_estimateX,2);
    Ms = [3 5]; %marker sizes
    c_list = ['r' 'b' 'g' 'c' 'm' 'y'];
    for Dc = 1:nF
        if ~isnan(Q_loc_estimateX(t,Dc))
            Sz = mod(Dc,2)+1; %pick marker size
            Cz = mod(Dc,6)+1; %pick color
            tmX = Q_loc_estimateX(1:t,Dc);
            tmY = Q_loc_estimateY(1:t,Dc);
            plot(tmX,tmY,'.-','markersize',Ms(Sz),'color',c_list(Cz),'linewidth',1)
            axis on;
        end
    end
    t
    
end

fontsize =12;
for Dc = 1:nF
    if ~isnan(Q_loc_estimateX(41,Dc))                %frame number
    tmX = Q_loc_estimateX(41,Dc);
    tmY = Q_loc_estimateY(41,Dc);
    a = double(tmX);
    b = double(tmY);
   text(a,b,num2str(Dc),'FontSize', fontsize, 'FontWeight', 'Bold','Color', [1 1 1]);
    end
    axis on;
end



final_tracks=[];
for i=1:nF
tempx=Q_loc_estimateX(:,i);
tempx(isnan(tempx))=[];
tempy=Q_loc_estimateY(:,i);
tempy(isnan(tempy))=[];
final_tracks=[final_tracks;[ones(size(tempx))*i,tempx,tempy]];
end



save('position_estimates.mat', 'Q_loc_estimateX', 'Q_loc_estimateY')