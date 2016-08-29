clc; % Clear command window.
clear all; % Get rid of variables from prior run of this m-file.
disp('DETECTION');
workspace; % Show panel with all the variables.

%set(0,'DefaultFigureWindowStyle','docked') %dock the figures..just a personal preference you don't need this.

base_dir = '/Users/snehabelkhale/Documents/MATLAB/tracking'; %change the directiory


cd(base_dir);

%have a folder in base_dir called 'jp5' to hold your .tiff frames
%%get listing of frames so that you can cycle through them easily. 

f_list =  dir('jp5');
f_list(1:3,:)=[];

%% manually change the vaule in line 30 to find an appropriate threshold value

hsizeh = 4;  %we might need to change these two values but so far they are good for all region analysis
sigmah = 0.43;
h = fspecial('log', hsizeh, sigmah);

%{
img =imread(f_list(1).name);
im=mat2gray(img);
subplot(1,2,1),imshow(im);
im1 = conv2(im,h,'same');
% imwrite(im1 , [f_list(1).name]);
im1 = im2bw(im1,0.2);              %0.2 is usually good for detection
subplot(1,2,2),imshow(im1);
%}

%% after finding the appropriate threshold value, 

for i = 1:length(f_list)
    img =imread(f_list(i).name);
    img=img(:,:,1);
    im=im2double(img);

    %if you want to compare LoG filter to Hough Transform, uncomment this section 

    %{
    [centersBright, radiiBright] = imfindcircles(im,[5 35],'ObjectPolarity','bright','EdgeThreshold',0.1);
    figure
    imshow(im)
    viscircles(centersBright, radiiBright,'EdgeColor','b');
    im=mat2gray(im);
    im = imsharpen(im,'Radius',2,'Amount',1); 
    %}

    %convolute image with the LoG Filter 
    im1 = conv2(im,h,'same');

    %figure
    %imshow(im1)

    im1 = im2bw(im1,0.15);  % plug in the threshold from above
 
    %use this section if the particles are not generally too close to eachother
    im1 = bwmorph(im1,'close');
    im1 = bwmorph(im1,'bridge'); 
    im1 = bwmorph(im1,'fill'); 
    im1 = bwmorph(im1,'close');  
    im1 = bwmorph(im1,'bridge'); 
    im1 = bwareaopen(im1, 4); %remove detections that are less than 4 pixels, because the size of the vesicle is around 10-15 pixels
    im1 = bwdist(im1) <= 2; %merge detections that are closer than 2 pixels
    im1 = imfill(im1,'holes');
    %

    %figure
    %imshow(im1)
    im1 = bwmorph(im1,'shrink',Inf); %shrink detection to one point
    %figure
    %imshow(im1)
 
    %%use extrema and extrema2 to find the local extrema and plot them
    [zmax,imax,zmin,imin] = extrema2(im1);
    [Y{i},X{i}] = ind2sub(size(im1),imax);
    %Y{i}=centersBright(:,2);  
    %X{i}=centersBright(:,1); 
 
    %remove points that are too close together, change this parameter
    %depending on data set
    %{
    est_dist = pdist([Y{i},X{i}]);
    est_dist = squareform(est_dist); %make square
    [row,col]=find(est_dist<5 & est_dist>0);
    row=unique(row);
    if ~isempty(row)
        row(1)=[];
    end
    col=unique(col);
    if ~isempty(col)
        col(1)=[];
    end
    Y{i}(row)=[];
    X{i}(col)=[];
    %}
    
 
    %for plotting, we can check the detection points by run the following
    %part
  
    figure,imagesc(im1)
        axis off
    figure,imshow(img)
    
    hold on
    
    for j = 1:length(X{i})
        plot(X{i}(j),Y{i}(j),'or')
    end
    axis off

end;

%SECTION FOR FINDING BOUNDARY FROM FIRST FRAME, and removing detections that are outside the boundary
%{
i=1;
img =imread(f_list(i).name);
img=img(:,:,1);
im=im2double(img);
BW = imfill(im,'holes');

level=0.3;
BW = im2bw(BW, level);
BW = bwdist(BW) <= 2;
figure
imshow(BW)
hold on
boundaries = bwboundaries(BW);
for k=1
   b = boundaries{k};
   plot(b(:,2),b(:,1),'g','LineWidth',3);
end
%%
%remove detections outside the boundary
for j=1:size(X,2)
    in=inpolygon(X{j}',Y{j}',boundaries{1}(:,2)',boundaries{1}(:,1)');
    X{j}=X{j}(in);
    Y{j}=Y{j}(in);
    %hold on
%plot(X{j}(in),Y{j}(in),'r+') % points inside
%plot(X{j}(~in),Y{j}(~in),'bo') % points outside
%hold off
end

%}
%%
save('vesicles_detections.mat',  'X','Y')   