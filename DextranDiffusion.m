%% Dextran perfusion


%% Load data
% Add bioformats
addpath('C:\Users\Rebek\OneDrive\Documenten\MATLAB\BioFormats');
javaaddpath('C:\Users\Rebek\OneDrive\Documenten\MATLAB\BioFormats\bioformats_package.jar');

% Read lif file with images of dextran diffusion
lifFilePath = 'C:\Users\Rebek\OneDrive\Documenten\MasterInternship\PermeabilityData\g6_day3_USMB TRITC.lif';

% Open the .lif file using Bio-Formats
reader = bfGetReader(lifFilePath);

%% Store all images in a list
% Access the metadata store to get information about the series (images)
metadata = reader.getMetadataStore();
numSeries = reader.getSeriesCount();

% Names of interest
targetNames = {'g6_TRITC_GFP_after USMB_t_0_10 min_Merged', ...
               'g6_TRITC_GFP_after USMB_t_10_20 min_Merged'};

% Initialize lists for storing images
imageList_0_10 = cell(83, 2);  % List for first series
imageList_10_20 = cell(83, 2); % List for second series

% Loop through series
for s = 1:numSeries
    reader.setSeries(s-1); % Zero-based indexing
    
    % Get the series name
    seriesName = char(metadata.getImageName(s-1));
    
    % Check if it's one of the target names
    if strcmp(seriesName, targetNames{1})
        targetList = imageList_0_10; % For first target series
    elseif strcmp(seriesName, targetNames{2})
        targetList = imageList_10_20; % For second target series
    else
        continue; % Skip if not an image of interest
    end
    
    % Get number of timepoints
    numTimepoints = reader.getSizeT();
    
    % Loop through all 83 timepoints
    for t = 1:numTimepoints
        % Read image from first channel (channel index 0 in zero-based)
        imgMatrix_channel1 = bfGetPlane(reader, 1 + (t-1) * reader.getSizeC());  % Channel 1 (first channel)
        
        % Read image from second channel (channel index 1 in zero-based)
        imgMatrix_channel2 = bfGetPlane(reader, 2 + (t-1) * reader.getSizeC());  % Channel 2 (second channel)
        
        % Store the first channel image in the first column and the second channel image in the second column
        targetList{t, 1} = imgMatrix_channel1; % First channel in column 1
        targetList{t, 2} = imgMatrix_channel2; % Second channel in column 2
    end
    
    % Assign the list back to the corresponding image list
    if strcmp(seriesName, targetNames{1})
        imageList_0_10 = targetList; % Store combined list for first target series
    else
        imageList_10_20 = targetList; % Store combined list for second target series
    end
end

%% Explore image 0-10
%Select an image of the imagelist to view
imgMatrix = imageList_0_10{1, 2};  % t = 0 min

% Display the TRITC image
figure;
imagesc(imgMatrix, [0 8200]);  % Show image with scaled colors
colormap([linspace(0,1,256)', zeros(256,1), zeros(256,1)]); % Red colormap
colorbar; % Show intensity scale
axis image; % Keep aspect ratio
title('Dextran diffusion t=0 min');


% Display the GFP image
% imgMatrix = imageList_0_10{8, 1};  
% figure
% imagesc(imgMatrix, [0 9000]);  % Show image with scaled colors
% colormap([zeros(256,1), linspace(0,1,256)', zeros(256,1)]); % Green colormap
% colorbar; % Show intensity scale
% axis image; % Keep aspect ratio
% title('Microvessel t=0 min');

%% Explore image 10-20
%Select an image of the imagelist to view
imgMatrix = imageList_10_20{83, 2};  % t = 20 min

% Display the image
figure;
imagesc(imgMatrix, [0 8200]);  % Show image with scaled colors
colormap([linspace(0,1,256)', zeros(256,1), zeros(256,1)]); % Red colormap
colorbar; % Show intensity scale
axis image; % Keep aspect ratio
title('Dextran diffusion t=20 min');



%% Get a cross section image 
% 0 minute is frame 1 of _0_10
img_t1 = imageList_0_10{1, 2};  

% select the middle of the vessel
crssec_t1 = img_t1(:, 975); %ncolums / 2

% get rownumber as axis
[nrow, ncol] = size(crssec_t1);

% Translate to real distance
height_micron = 6586;
height_micron_pixel = height_micron / nrow;
array_height = (1:nrow)*height_micron_pixel;

% Plot the fluorescent cross section
figure
plot(array_height, crssec_t1(:, 1), 'r' )
xlim([0 6300])
ylim([0 6100])
hold on

% 20 minutes is frame 83 of _10_20
img_t20 = imageList_10_20{83, 2};  

% select the middle of the vessel
crssec_t15 = img_t20(:, 975);

% Add it to the plot
plot(array_height, crssec_t15(:,1), 'b')
xline(1037*height_micron_pixel, 'g', 'LineWidth', 1)
xline(1391*height_micron_pixel, 'g', 'LineWidth', 1)
legend('t = 0 min',   't = 20 min', 'FontSize', 11)
xlabel('μm', 'FontSize', 14)
ylabel('Fluorescent intensity (a.u.)', 'FontSize', 14)
title('Cross-diameter TRITC-dextran distribution', 'FontSize', 14)
set(gca, 'FontSize', 13)



%% Compute Permeability in the middle
% The values in these first 2 lines are based on dividing the entire FOV in
% 5 regions. For 'middle' the middle region is analyzed. For left and right
% the outer regions are analyzed
middle_width = 890;
middle_height = 1266;
w = 356;  % Width of FOV
h = 2500;  % Height of FOV
x1 = middle_width - w/2;  % X start (column index)
y1 = middle_height - h/2;  % Y start (row index)
x2 = x1 + w - 1;  % X end
y2 = y1 + h - 1;  % Y end

% View Image at t=0
imgMatrix = imageList_0_10{1, 2};  % t = 0 min
figure;
imagesc(imgMatrix, [0 8200]);  % Show image with scaled colors
colormap([linspace(0,1,256)', zeros(256,1), zeros(256,1)]); % Red colormap
colorbar; 
axis image;
title('Dextran diffusion t=0 min');
hold on;

% Draw FOV rectangle
rectangle('Position', [x1, y1, w, h], 'EdgeColor', 'y', 'LineWidth', 1);

% Extract FOV and compute mean intensity
FOV_region_0_10 = imgMatrix(y1:y2, x1:x2);
mean_intensity_0 = mean(FOV_region_0_10(:));
I1 = mean_intensity_0;

fprintf('Mean intensity (0-10 min): %.2f\n', mean_intensity_0);


% View Image at t=20
imgMatrix = imageList_10_20{83, 2};  % t = 20 min
figure;
imagesc(imgMatrix, [0 8200]); 
colormap([linspace(0,1,256)', zeros(256,1), zeros(256,1)]); % Red colormap
colorbar; 
axis image;
title('Dextran diffusion t=20 min');
hold on;

% Draw FOV rectangle
rectangle('Position', [x1, y1, w, h], 'EdgeColor', 'y', 'LineWidth', 1);

% Extract FOV and compute mean intensity
FOV_region_10_20 = imgMatrix(y1:y2, x1:x2);
mean_intensity_20 = mean(FOV_region_10_20(:));
I2 = mean_intensity_20;

fprintf('Mean intensity (10-20 min): %.2f\n', mean_intensity_20);

Ib = 0;
dt = 20*60; %seconds
d = 870e-4; % cm 
p_middle = 1/(I1-Ib) * ((I2-I1)/dt) * (d/4);

%% Compute Permeability on the left
middle_width = 179;
middle_height = 1266;
w = 356;  % Width of FOV
h = 2500;  % Height of FOV
x1 = middle_width - w/2;  % X start (column index)
y1 = middle_height - h/2;  % Y start (row index)
x2 = x1 + w - 1;  % X end
y2 = y1 + h - 1;  % Y end

% View Image at t=0
imgMatrix = imageList_0_10{1, 2};  % t = 0 min
figure;
imagesc(imgMatrix, [0 8200]);  % Show image with scaled colors
colormap([linspace(0,1,256)', zeros(256,1), zeros(256,1)]); % Red colormap
colorbar; 
axis image;
title('Dextran diffusion t=0 min');
hold on;

% Draw FOV rectangle
rectangle('Position', [x1, y1, w, h], 'EdgeColor', 'y', 'LineWidth', 1);

% Extract FOV and compute mean intensity
FOV_region_0_10 = imgMatrix(y1:y2, x1:x2);
mean_intensity_0 = mean(FOV_region_0_10(:));
I1 = mean_intensity_0;

fprintf('Mean intensity (0-10 min): %.2f\n', mean_intensity_0);


% View Image at t=20
imgMatrix = imageList_10_20{83, 2};  % t = 20 min
figure;
imagesc(imgMatrix, [0 8200]); 
colormap([linspace(0,1,256)', zeros(256,1), zeros(256,1)]); % Red colormap
colorbar; 
axis image;
title('Dextran diffusion t=20 min');
hold on;

% Draw FOV rectangle
rectangle('Position', [x1, y1, w, h], 'EdgeColor', 'y', 'LineWidth', 1);

% Extract FOV and compute mean intensity
FOV_region_10_20 = imgMatrix(y1:y2, x1:x2);
mean_intensity_20 = mean(FOV_region_10_20(:));
I2 = mean_intensity_20;

fprintf('Mean intensity (10-20 min): %.2f\n', mean_intensity_20);

Ib = 0;
dt = 20*60; %seconds
d = 870e-4; % cm 
p_left = 1/(I1-Ib) * ((I2-I1)/dt) * (d/4);

%% Compute Permeability on the right
middle_width = 1600;
middle_height = 1266;
w = 356;  % Width of FOV
h = 2500;  % Height of FOV
x1 = middle_width - w/2;  % X start (column index)
y1 = middle_height - h/2;  % Y start (row index)
x2 = x1 + w - 1;  % X end
y2 = y1 + h - 1;  % Y end

% View Image at t=0
imgMatrix = imageList_0_10{1, 2};  % t = 0 min
figure;
imagesc(imgMatrix, [0 8200]);  % Show image with scaled colors
colormap([linspace(0,1,256)', zeros(256,1), zeros(256,1)]); % Red colormap
colorbar; 
axis image;
title('Dextran diffusion t=0 min');
hold on;

% Draw FOV rectangle
rectangle('Position', [x1, y1, w, h], 'EdgeColor', 'y', 'LineWidth', 1);

% Extract FOV and compute mean intensity
FOV_region_0_10 = imgMatrix(y1:y2, x1:x2);
mean_intensity_0 = mean(FOV_region_0_10(:));
I1 = mean_intensity_0;

fprintf('Mean intensity (0-10 min): %.2f\n', mean_intensity_0);


% View Image at t=20
imgMatrix = imageList_10_20{83, 2};  % t = 20 min
figure;
imagesc(imgMatrix, [0 8200]); 
colormap([linspace(0,1,256)', zeros(256,1), zeros(256,1)]); % Red colormap
colorbar; 
axis image;
title('Dextran diffusion t=20 min');
hold on;

% Draw FOV rectangle
rectangle('Position', [x1, y1, w, h], 'EdgeColor', 'y', 'LineWidth', 1);

% Extract FOV and compute mean intensity
FOV_region_10_20 = imgMatrix(y1:y2, x1:x2);
mean_intensity_20 = mean(FOV_region_10_20(:));
I2 = mean_intensity_20;

fprintf('Mean intensity (10-20 min): %.2f\n', mean_intensity_20);

Ib = 0;
dt = 20*60; %seconds
d = 870e-4; % cm 
p_right = 1/(I1-Ib) * ((I2-I1)/dt) * (d/4);