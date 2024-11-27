% Developed by Yuntao Li, OMNI Lab, Northeastern University

% This is data analysis code for GCaMP vs. Pupil diameter

% The code presents trial-by-trial data, ROI on cortex, average data, and
% calculates rising time and onset time

% For further information, please reach out to: li.yunt@northeastern.edu
%%
clc; clear; close all

%% load data %%%%%%%%%%
% load path and find the number of images
[filename, pathname] = uigetfile('*.csv', 'Select One or More Files');

%% read data and divide into datasets based on every 30 seconds elapsed time
% Import data
csv = readtable(fullfile(pathname, filename));

% Extract relevant columns
data = table2array(csv(:, 2:5));
D = table2array(csv(:, 5));

% Convert to numeric arrays
acc = str2double(table2array(csv(:, 2)));
diameter = str2double(D);

% Extract elapsed time and pulse data
elapsetime = table2array(csv(:, 3));
pulse = table2array(csv(:, 4));

% Find indices where pulse data is 1
pulse_data = find(pulse == 1);

% Define the threshold values for elapsed time
thresholds = 0:30:max(elapsetime-1);

% Initialize cell array to store indices
indices = cell(length(thresholds), 1);

% Find the indices where elapsed time exceeds each threshold value
for i = 1:length(thresholds)
    indices{i} = find(elapsetime > thresholds(i), 1);
end

% Display the indices
for i = 1:length(thresholds)
    fprintf('Elapsed time exceeds %d seconds at index %d\n', thresholds(i), indices{i});
end

% Initialize cell array to store datasets
datasets = cell(length(thresholds), 1);

% Find the indices where elapsed time exceeds each threshold value
for i = 1:length(thresholds) - 1
    start_index = indices{i};
    end_index = indices{i + 1} - 1;
    datasets{i} = csv(start_index:end_index, :);
end

% Handle the remaining data
start_index = find(elapsetime > thresholds(end), 1);
datasets{end} = csv(start_index:end, :);

%% Plot data

% Initialize variables to store global minimum and maximum values
globalMinDia = inf;
globalMaxDia = -inf;

% Create new figure for OISI trials data
load intensityMatrix
fps = 5; 
dt = 1/fps;
time = (0:dt:30-dt);

% Determine the number of curves
numCurves = size(intensityMatrix, 2);

for i = 1:length(datasets)
    % Get current dataset
    dataset = datasets{i};
    
    % Extract data from dataset
    elapsedTime = dataset{:, 3};
    acceleration = dataset{:, 2}; % 9.75: empirical acceleration constant
    diameter = dataset{:, 5};
    pulse = dataset{:, 4};

    % Find the index where pulse signal turns from 0 to 1
    pulse_start_index = find(diff(pulse) == 1, 1) + 1;

    % Calculate the time point corresponding to 5 seconds into the dataset
    time_at_5s = elapsedTime(pulse_start_index) + 5;

    % Find the index corresponding to time_at_5s
    [~, idx_5s] = min(abs(elapsedTime - time_at_5s));

    % Calculate the average diameter from the start of the dataset to 5 seconds
    avg_diameter = mean(diameter(1:pulse_start_index));
    
    avg_basline_acceleration = 9.80665;% compensate the measurement of gravity in Z-diretion

    acceleration_nobaseline = acceleration - avg_basline_acceleration;

    % Calculate the average acceleration from the start of the dataset to 5 seconds
    avg_acceleration = mean(acceleration_nobaseline);
    std_acceleration = std(acceleration_nobaseline);

    % Define the threshold
    threshold = 0.2;
    
    abs_acceleration = abs(acceleration_nobaseline - avg_acceleration);

    % Find the points outside the range (greater than 0 + 0.3)
    outside_range_indices = find((abs_acceleration > threshold));
    


    % Initialize a new array of zeros
    new_array = zeros(size(abs_acceleration));
    new_value = 1;
    new_array(outside_range_indices) = new_value;

    % Define the maximum allowable gap to merge clusters (fps/4)
    max_gap = 15;
    neighborhood_size = 15;


    % Loop through the new_array to find isolated '1's
    for i_acc = 1:length(new_array)
        if new_array(i_acc) == new_value
            % Define the left and right boundaries for checking neighbors
            left_bound = max(1, i_acc - neighborhood_size);  % Ensure we don't go out of bounds on the left
            right_bound = min(length(new_array), i_acc + neighborhood_size);  % Ensure we don't go out of bounds on the right
            
            % Get the number of '1's in the neighborhood excluding the current index
            left_neighbors = sum(new_array(left_bound:i_acc-1) > 0);
            right_neighbors = sum(new_array(i_acc+1:right_bound) > 0);
            
            % If there is only one neighbor total (either on the left or right), remove the current '1'
            if left_neighbors + right_neighbors <=6
                new_array(i_acc) = 0;
            end
        end
    end

    % Find indices where new_array equals 1
    one_indices = find(new_array == new_value);
    
    % Loop through the indices and merge close clusters
    i_acc = 1;  % Start from the first '1'
    while i_acc < length(one_indices)
        % Check the distance between the current index and the next index
        if one_indices(i_acc + 1) - one_indices(i_acc) <= max_gap
            % If the gap is within the max_gap, fill the gap by setting values to 1
            new_array(one_indices(i_acc):one_indices(i_acc + 1)) = new_value;
        end
        i_acc = i_acc + 1;  % Move to the next index
    end

    figure('Position', [300, 400, 400, 200]);bar(elapsedTime,abs_acceleration,'b');
    ylim([0, 2]);
    hold on
    plot(elapsedTime,new_array);
    %plot(elapsedTime(outside_range_indices), abs_acceleration(outside_range_indices), 'ro', 'MarkerFaceColor', 'r', 'DisplayName', 'Out of Range');
   
    new_array = 10.*new_array;

    % calculate the absolution value of relative change of accleration
    relative_diameter = (diameter - avg_diameter)./avg_diameter*100;

    % butterworhth filter
    nyquist_freq = 0.5*length(relative_diameter)/30;
    normalized_cutoff_freq = 1 / nyquist_freq;
    [b,a] = butter(4, normalized_cutoff_freq, 'low');
    relative_diameter = filtfilt(b,a,relative_diameter);
    

    % % Scale acceleration data to [0 10]
    % accRange = max(relative_acceleration) - min(relative_acceleration);
    % diaRange = max(relative_diameter) - min(relative_diameter);
    
    % Update global minimum and maximum values
    globalMinDia = min(globalMinDia, min(relative_diameter));
    globalMaxDia = max(globalMaxDia, max(relative_diameter));

    % Define the desired width and height for the figure window
    width = 400;  % Width in pixels
    height = 200; % Height in pixels
    
    % Create the figure and set the position of the figure window
    figure('Position', [300, 100, width, height]);
    
    % Plot diameter against elapsed time
    yyaxis left
    h1 = plot(elapsedTime, relative_diameter, 'k','LineWidth',1.2);
    % Highlight the points outside the range in red

    %ylabel('Diameter (pxls)','Color','g');
    
    % Set the color of numbers on the y-axis for diameter
    set(gca, 'YColor', 'k');
    
    % Set y-axis limits for diameter
    if globalMinDia > 0
        ylim([globalMinDia * 0.95, globalMaxDia * 1.05]);
    else
        ylim([-abs(globalMinDia) * 1.05, globalMaxDia * 1.05]);
    end

    % Highlight pulse events within shaded area
    for j = 1:length(pulse)
        if pulse(j) == 1
            % Find start and end times of pulse event
            pulse_start = elapsedTime(j);
            pulse_end = elapsedTime(min(j + 1, length(elapsedTime)));
            
            % Draw shaded rectangle with partial transparency using 'patch'
            v1 = [pulse_start globalMinDia * 0.95;pulse_end globalMinDia * 0.95;pulse_end globalMaxDia*1.05;pulse_start globalMaxDia*1.05];
            f1 = [1 2 3 4];
            %x1 = [pulse_start, pulse_end, pulse_end, pulse_start];
            %y1 = [0, 0, globalMaxDia*1.05, globalMaxDia*1.05];
            p1 = patch('Faces', f1, 'Vertices', v1, 'EdgeColor', 'none', 'FaceColor', [1 0 0], 'FaceAlpha', 0.3);
        end
    end

    % Plot acceleration against elapsed time as the shaded area
    for j = 1:length(abs_acceleration)
        % Find start and end times of acceleration event
        delta_eplapsedTime = 30/length(elapsedTime);
        acceleration_start = elapsedTime(j) - delta_eplapsedTime;
        acceleration_end = elapsedTime(j) + delta_eplapsedTime;
        
        % Draw shaded rectangle with partial transparency using 'patch'
        %v2 = [acceleration_start globalMinDia * 0.95;acceleration_end globalMinDia * 0.95;acceleration_end globalMinDia * 0.95 + relative_acceleration(j);acceleration_start globalMinDia * 0.95 + relative_acceleration(j)];
        % v2 = [acceleration_start 0;acceleration_end 0;acceleration_end relative_acceleration(j);acceleration_start relative_acceleration(j)];
        % f2 = [1 2 3 4];
        % %x2 = [acceleration_start, acceleration_end, acceleration_end, acceleration_start];
        % %y2 = [0, 0, acceleration_scaled(j), acceleration_scaled(j)];
        % if relative_acceleration(j) > 1
        %     p2 = patch('Faces', f2, 'Vertices', v2, 'EdgeColor', 'none', 'FaceColor', [1 0 0], 'FaceAlpha', 0.7);
        % else
        %     p2 = patch('Faces', f2, 'Vertices', v2, 'EdgeColor', 'none', 'FaceColor', [0 0 1], 'FaceAlpha', 0.7);
        % end
 
    end
    
    % Create a new vector of time values
    yyaxis right
    numPoints = length(elapsedTime);
    newTimeVector = linspace(0, max(time), numPoints);
    
    % Interpolate the intensity data corresponding to the new time values
    newIntensityMatrix = zeros(numPoints, size(intensityMatrix(:,i), 2)); % Initialize interpolated intensity matrix
    for k = 1:size(intensityMatrix(:,i), 2)
        newIntensityMatrix(:, k) = interp1(time, intensityMatrix(:,i), newTimeVector, 'linear');
    end

    % Plot the current curve
    h2 = plot(elapsedTime, newIntensityMatrix, 'r', 'LineWidth',1.2);
    %ylabel('Norm Intensity (a.u.)','Color','r');

    % Set y-axis limits for OISI intensity
    ylim([min(newIntensityMatrix) - 0.005, max(newIntensityMatrix) + 0.005]);

    % Set the color of numbers on the y-axis for diameter
    set(gca, 'YColor', 'r');

    % Set x-axis limits
    xlim([min(elapsedTime), max(elapsedTime)]);
    
    hold off
    % Add labels and title
    %xlabel('Elapsed Time (s)');
    %title(['Trial ' num2str(i)]);
    %set(gca,'fontsize',10);

    % save pictures (comment when dont need)
    %print('-painters', '-dsvg', num2str(i));
end

%% plot ROI with activation area
load data_fractional_change_all
load ROI
I=medfilt2(data_fractional_change_all,[3,3]);
figure;
imagesc(I);
%imcontrast();
axis image; 
axis off;
hold on
set(gca,'fontsize',12)
clim([-0.03 0.03])

rect_pos = ROI.Position;
rectangle('Position',rect_pos,'EdgeColor','r','LineWidth',2)

%% average plot
D1 = datasets{1};
D1 = D1{:, 5};

D2 = datasets{2};
D2 = D2{:, 5};

D3 = datasets{3};
D3 = D3{:, 5};

D4 = datasets{4};
D4 = D4{:, 5};

D5 = datasets{5};
D5 = D5{:, 5};

% Concatenate matrices into a single matrix
min_length = min([length(D1), length(D2), length(D3), length(D4), length(D5)]);
D1 = D1(1:min_length);
D2 = D2(1:min_length);
D3 = D3(1:min_length);
D4 = D4(1:min_length);
D5 = D5(1:min_length);

% Concatenate the matrices into a single vector
concatenated_vector = [D1, D2, D5];

% Calculate the average diameter from the start of the dataset to 5 seconds
avg_diameter = mean(concatenated_vector(1:180,:));


% calculate the absolution value of relative change of accleration
relative_diameter = (concatenated_vector - avg_diameter)./avg_diameter*100;

% butterworhth filter
nyquist_freq = 0.5*length(relative_diameter)/90;
normalized_cutoff_freq = 1 / nyquist_freq;
[b,a] = butter(4, normalized_cutoff_freq, 'low');
relative_diameter = filtfilt(b,a,relative_diameter);
save('relative_diameter.mat','relative_diameter');

mean_Diameter = mean(relative_diameter,2);

intensityMatrix_selected = [intensityMatrix(:,1),intensityMatrix(:,2),intensityMatrix(:,5)];
save('intensityMatrix_selected.mat','intensityMatrix_selected');
mean_Intensity = mean(intensityMatrix_selected,2);

figure('Position', [300, 100, 400, 200]);
% Plot mean_Intensity
yyaxis left
numPoints = length(mean_Diameter);
newTimeVector = linspace(0, max(time), numPoints)';

% Interpolate the intensity data corresponding to the new time values
Interepolate_newIntensity = zeros(numPoints, size(mean_Intensity(:,1), 2)); % Initialize interpolated intensity matrix
Interepolate_newIntensity(:, 1) = interp1(time, mean_Intensity(:,1), newTimeVector, 'linear');
plot(Interepolate_newIntensity, 'g'); % Blue line

% Hold on to the existing plot
hold on;

% Plot mean_D
yyaxis right
plot(mean_Diameter, 'k'); % Red line

%% calculate the rising time and onset time
% butterworhth filter
nyquist_freq = 0.5*length(elapsedTime)/30;
normalized_cutoff_freq = 1 / nyquist_freq;
[b,a] = butter(4, normalized_cutoff_freq, 'low');

RisingandFalling_time = zeros(5,4);
onsettime = zeros(5,2);

for i = 1:length(datasets)
    % Get current dataset
    dataset = datasets{i};
    
    % Extract data from dataset
    elapsedTime = dataset{:, 3};
    acceleration = dataset{:, 2}; % 9.75: empirical acceleration constant
    diameter = dataset{:, 5};
    pulse = dataset{:, 4};

    % Find the index where pulse signal turns from 0 to 1
    pulse_start_index = find(diff(pulse) == 1, 1) + 1;

    % apply filter to diameter
    diameter_fil = filtfilt(b,a,diameter);
    
    % diameter curve fitting
    diameter_fit = fit(elapsedTime,diameter_fil,'gauss6');
    diameter_fit = diameter_fit(elapsedTime);
    figure;
    plot(elapsedTime, diameter_fil,'b');
    hold on

    % diameter rising and onset time
    [pks,locs,w,p,wxPk] = findpeaks2(diameter_fit,elapsedTime,'NPEAKS',1,'SORTSTR','descend','Annotate','extents');
    findpeaks2(diameter_fit,elapsedTime,'NPEAKS',1,'SORTSTR','descend','Annotate','extents');
    diameter_10_index = find(diameter_fit(pulse_start_index:end)>diameter_fit(pulse_start_index)+((pks - diameter_fit(pulse_start_index))*0.1),1);
    diameter_90_index = find(diameter_fit(pulse_start_index:end)>diameter_fit(pulse_start_index)+((pks - diameter_fit(pulse_start_index))*0.9),1);
    diameter_rising_time = elapsedTime(pulse_start_index+diameter_90_index) - elapsedTime(pulse_start_index+diameter_10_index);
    diameter_falling_time = wxPk(2) - locs;
    diameter_onsettime = elapsedTime(pulse_start_index+diameter_10_index) - elapsedTime(pulse_start_index);
    hold on  % to plot on the current figure
    plot(elapsedTime(pulse_start_index+diameter_10_index),diameter_fit(pulse_start_index+diameter_10_index),'r*')  % adds a red asterisk at the point (x_pos,y_pos)
    plot(wxPk(2),diameter_fit(find(elapsedTime > (wxPk(2)),1)-1),'r*')  % adds a red asterisk at the point (x_pos,y_pos)
    plot(elapsedTime(pulse_start_index+diameter_90_index),diameter_fit(pulse_start_index+diameter_90_index),'r*')  % adds a red asterisk at the point (x_pos,y_pos)
    plot(elapsedTime(pulse_start_index),diameter_fit(pulse_start_index),'b*')  % adds a red asterisk at the point (x_pos,y_pos)
    hold off
    legend('filtered signal','fitted signal', 'Peak','width','Peak height','Rising point','falling point');
    
    % create new time vector for OISI data
    numPoints = length(elapsedTime);
    newTimeVector = linspace(0, max(time), numPoints);
    
    % Interpolate the intensity data corresponding to the new time values
    newIntensityMatrix = zeros(numPoints, size(intensityMatrix(:,i), 2)); % Initialize interpolated intensity matrix
    for k = 1:size(intensityMatrix(:,i), 2)
        newIntensityMatrix(:, k) = interp1(time, intensityMatrix(:,i), newTimeVector, 'linear');
    end
    
    % OISI curve fitting
    upsidedown_newIntensityMatrix = 1./newIntensityMatrix;
    Intensity_fit = fit(elapsedTime,upsidedown_newIntensityMatrix,'gauss8');
    Intensity_fit = Intensity_fit(elapsedTime);
    figure;
    plot(elapsedTime, upsidedown_newIntensityMatrix,'b');
    hold on

    % OISI rising and onset time
    [pks,locs,w,p,wxPk] = findpeaks2(Intensity_fit,elapsedTime,'NPEAKS',1,'SORTSTR','descend','Annotate','extents');
    findpeaks2(Intensity_fit,elapsedTime,'NPEAKS',1,'SORTSTR','descend','Annotate','extents');
    onsetIntensity10 = Intensity_fit(pulse_start_index) + (pks - Intensity_fit(pulse_start_index))*0.1;
    onsetIntensity90 = Intensity_fit(pulse_start_index) + (pks - Intensity_fit(pulse_start_index))*0.9;
    Intensity_10_index = find(Intensity_fit(pulse_start_index:end) > onsetIntensity10, 1) + pulse_start_index;
    Intensity_90_index = find(Intensity_fit(pulse_start_index:end) > onsetIntensity90, 1) + pulse_start_index;
    Intensity_rising_time = elapsedTime(Intensity_90_index) - elapsedTime(Intensity_10_index);
    Intensity_falling_time = wxPk(2) - locs;
    Intensity_onsettime = elapsedTime(Intensity_10_index) - elapsedTime(pulse_start_index);
    plot(elapsedTime(Intensity_10_index),Intensity_fit(Intensity_10_index),'r*')  % adds a red asterisk at the point (x_pos,y_pos)
    plot(elapsedTime(Intensity_90_index),Intensity_fit(Intensity_90_index),'r*')  % adds a red asterisk at the point (x_pos,y_pos)
    plot(wxPk(2),Intensity_fit(find(elapsedTime > (wxPk(2)),1)-1),'r*')  % adds a red asterisk at the point (x_pos,y_pos)
    plot(elapsedTime(pulse_start_index),Intensity_fit(pulse_start_index),'b*')
    hold off
    legend('filtered signal','fitted signal', 'Peak','width','Peak height','Rising point','falling point');
    
    RisingandFalling_time(i,:) = [diameter_rising_time, diameter_falling_time, Intensity_rising_time, Intensity_falling_time];
    onsettime(i,:) = [diameter_onsettime, Intensity_onsettime]; 
end
