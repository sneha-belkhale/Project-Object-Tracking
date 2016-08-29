
%now, we try the Bayesian approach
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
%how far away is an acceptable match?
MUNKRES_REJECT=20; 
%how far away is acceptable for two particles to be matched to the same
%detection?
COLLISION_REJECT=10; 
%how many frames of no matching can pass until we remove a detection
REMOVE_FRAME=5; 
%initial velocity of particles should always be 1 if we assume constant
INITIAL_VEL=1;
S_frame = 1; %starting frame 
E_frame = 14;%ending frame

%% initize result variables
Q_loc_meas = []; % the fly detecions  extracted by the detection algo

%% initialize estimation variables for two dimensions
Q= [X{S_frame} Y{S_frame} zeros(length(X{S_frame}),1) zeros(length(X{S_frame}),1)]'

Q_estimate = nan(4,100);
Q_estimate(:,1:size(Q,2)) = Q;  %estimate of initial location estimation of where the flies are(what we are updating)
Q_loc_estimateY = nan(100); %  position estimate
Q_loc_estimateX= nan(100); %  position estimate
strk_trks = zeros(1,100);  %counter of how many strikes a track has gotten
nD = size(X{S_frame},1); %initize number of detections
nF =  find(isnan(Q_estimate(1,:))==1,1)-1;  %initize number of track estimates
vel=ones(nF,1)*INITIAL_VEL;

%for each frame
for t = S_frame:E_frame
    
    % load the image
    img_tmp = double(imread(f_list(t).name));
    img = img_tmp(:,:,1);
    %get size
    [imsize_x, imsize_y,z]=size(img);
    % make the given detections matrix
    Q_loc_meas = [X{t} Y{t}];
    %set new number of detections
    nD = size(X{t},1); 
    
    %% now we assign the detections to estimated track positions 
    %make the probability (cost) matrice between all pairs rows = tracks, coln =
    %detections
    est_dist=zeros(nF,nD);
    for i=1:nF
        %calculate the new velocity matrix for each particle  
        K=[vel(i)^2,0;0,vel(i)^2];
        for j=1:nD
            %calculate the probability that the particles could be linked,
            %assuming gaussian distribution 
            est_dist(i,j)=1/sqrt((2*pi)^2*det(K)) * exp(-(Q_loc_meas(j,:)'-Q_estimate(1:2,i))'*inv(K)*(Q_loc_meas(j,:)'-Q_estimate(1:2,i))/2)*10000;

        end
    end

  
    [asgn, cost] =munkres(-est_dist); 
    asgn = asgn';
    
    %% check some special situations
    
    %check 1: is the detection far from the observation? if so, reject it.
    rej = zeros(size(asgn));  
    for F = 1:nF
        if asgn(F) > 0
            rej(F)=pdist([Q_loc_meas(asgn(F),:);Q_estimate(1:2,F)'])<MUNKRES_REJECT;
        else
            rej(F) = 0;
        end
    end
    asgn = asgn.*rej;
    
    %For the estimates that did not get assigned, account for collisions by
    %checking if the new track is close enough to a previous estimate
    no_trk_list =  find(asgn==0);
    if (~isempty(no_trk_list)&~isempty(Q_loc_meas))
        %find the nearest neighbors of the unassigned tracks 
        asgn_dup=knnsearch(Q_loc_meas,Q_estimate(1:2,no_trk_list)');
        rej = zeros(size(asgn_dup));  
   
        for F = 1:size(asgn_dup)
            if asgn_dup(F) > 0
                rej(F)=pdist([Q_loc_meas(asgn_dup(F),:);Q_estimate(1:2,no_trk_list(F))'])<COLLISION_REJECT;
            else
                rej(F) = 0;
            end
        end
        asgn_dup = asgn_dup.*rej;
        %update our assignments 
        asgn(no_trk_list)=asgn_dup;
    end   
    
    %% apply the assignment to the update
    k = 1;
    for F = 1:length(asgn)
        if asgn(F) > 0
            %Q_estimate(:,k) = Q_estimate(:,k) + K * (Q_loc_meas(asgn(F),:)' - C * Q_estimate(:,k));
            z=pdist([Q_loc_meas(asgn(F),:);Q_estimate(1:2,k)']);
            if (z>0)
                %average the past and current velocity
                vel(k)=(vel(k)+pdist([Q_loc_meas(asgn(F),:);Q_estimate(1:2,k)']))/2;
            end
            Q_estimate(1:2,k) = Q_loc_meas(asgn(F),:)';

        end
        k = k + 1;
    end
    
    
    %% Store data
    Q_loc_estimateX(t,1:nF) = Q_estimate(1,1:nF);
    Q_loc_estimateY(t,1:nF) = Q_estimate(2,1:nF);
    
    %% ok, now that we have our assignments and updates, lets find the new detections and lost trackings
    
    %find the new detections. basically, anything that doesn't get assigned
    %is a new tracking
    new_trk = Q_loc_meas(~ismember(1:size(Q_loc_meas,1),asgn),:)';
    
    if ~isempty(new_trk)
        Q_estimate(:,nF+1:nF+size(new_trk,2))=  [new_trk; zeros(2,size(new_trk,2))];
        nF = nF + size(new_trk,2);  % number of track estimates with new ones included
        vel=[vel;ones(1,size(new_trk,2))'*INITIAL_VEL];
    end

    %give a strike to any tracking that didn't get matched up to a
    %detection
    no_trk_list =  find(asgn==0);
    
    if ~isempty(no_trk_list)
        strk_trks(no_trk_list) = strk_trks(no_trk_list) + 1;
    end
    
    %if a track has a strike greater than 30, delete the tracking. i.e.
    bad_trks = find(strk_trks >REMOVE_FRAME);                   
    Q_estimate(:,bad_trks) = NaN;
    
    
    %% plot the resulting trajectories
    clf
    img = imread(f_list(t).name);
    img=img(:,:,1); %for jp3
    im=im2double(img); %for jp3
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

%% save the trajectories to a csv file
final_tracks=[];
for i=1:nF
    tempx=Q_loc_estimateX(:,i);
    tempx(isnan(tempx))=[];
    tempy=Q_loc_estimateY(:,i);
    tempy(isnan(tempy))=[];
    final_tracks=[final_tracks;[ones(size(tempx))*i,tempx,tempy]];
end
save('position_estimates.mat', 'Q_loc_estimateX', 'Q_loc_estimateY')