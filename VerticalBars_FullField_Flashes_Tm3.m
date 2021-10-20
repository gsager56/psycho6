%% PARAMETERS YOU NEED TO SPECIFY
clear; clc;
param.movieChan = 1; % which channels do you want to analyze in the movie
param.probe_epochs = 1 : 7;
param.epochs_for_selectivity = [];%[1, 4:7];
param.interleave_epochs = [8];

% parameters for the ROi selection
param.corr_thresh = 0.5; % correlation threshold
param.probe_correlation = true; % do you want to use probe correlations for ROI selection
param.probe_correlation_type = 'hard code indices'; % type of correlation comparison for the 
param.corr_idxs{1} = 205 : 567; % indices of first probe I want to use for correlation
param.corr_idxs{2} = 1821 : 2183; % indices of 2nd probe I want to use for correlation
param.num_corr_epochs = round(length(param.epochs_for_selectivity) / 2); % number of epochs in param.epochs_for_selectivity above correlation threshold

param.force_alignment = false; % force alignment calculation?
param.manual_roi = true;

% how would you like to do ROI selection?

%% load in the experimental details

cellType='Tm3';
%stim='GreyBlink_Probe_LongGreyBlink';
stim='VerticalBars_FullField_Flashes';
sensor='gc6f';
%sensor='GC6f';
flyEye='right';
surgeon='Allison';

dataPath=GetPathsFromDatabase(cellType,stim,sensor,flyEye,surgeon);

sensor='GC6f';

dataPath_2=GetPathsFromDatabase(cellType,stim,sensor,flyEye,surgeon);
dataPath{3,1} = dataPath_2{1};
param.analyze_these = 3;%[2 : length(dataPath)];

%% loop through each recording to be analyzed

for i_ex = param.analyze_these
    
    %% load in experimental parameters and raw movie
    [exp, raw_movie] = get_exp_details( dataPath{i_ex}, param);

    %% register the relaxed mean image with each frame to get the alignment
    % find dimensions of image
    [n_rows, n_cols, n_t] = size(raw_movie);
    
    % average intensity of each pixel
    mean_movie = mean(raw_movie,3);
    
    alignment_folder = fullfile(dataPath{i_ex}, 'psycho6');
    
    if isfolder(alignment_folder) && ~param.force_alignment
        % download what you need to do alignment from your last calculation
        load( fullfile( alignment_folder, 'alignment_info.mat' ) );
        
        % only save the the rows and columns that stayed in frame
        good_cols = [abs(min(dx)) + 1 : n_cols - max(dx)]; % columns that stayed in frame
        good_rows = [max(dy) + 1 : n_rows - abs(min(dy))]; % rows that stayed in frame
        
        % initialize the aligned_movie
        aligned_movie = zeros( length(good_rows), length(good_cols), n_t);
        
        % iterate through time
        for i_t = 1 : n_t
            aligned_movie(:,:,i_t) = raw_movie(good_rows - dy(i_t) , good_cols + dx(i_t) , i_t);
        end
    else
        % do the alignment
        
        [dx, dy, aligned_movie, tissue_mask] = movie_alignment_registration( raw_movie, dataPath{i_ex} );
        
        % only save the the rows and columns that stayed in frame
        good_cols = [abs(min(dx)) + 1 : n_cols - max(dx)]; % columns that stayed in frame
        good_rows = [max(dy) + 1 : n_rows - abs(min(dy))]; % rows that stayed in frame

        aligned_movie = aligned_movie(good_rows, good_cols, :);
    end
    
    % redefine n_rows and n_cols to be size of cropped movie
    [n_rows, n_cols, ~] = size(aligned_movie);
    
    %% Subtract out background in aligned_movie
    
    % initiliaze the filtered movie
    filtered_movie = zeros(size(aligned_movie));
    
    % define background as below median, temporal mean intensity of the 0 pixels in
    % tissue_mask
    mean_movie = mean(aligned_movie, 3); % recalculate mean movie on the aligned movie
    bkg_mask = (mean_movie < median(mean_movie(~tissue_mask))) & ~tissue_mask;
    
    % loop through time
    for i = 1 : n_t
        this_frame = aligned_movie(:,:,i);
        filtered_movie(:,:,i) = aligned_movie(:,:,i) - mean(this_frame(bkg_mask));
    end
    
    mean_movie = mean(filtered_movie, 3); % recalculate mean movie on the fully filtered movie
    
    disp('Finished Processing Raw Movie')
    
    %% Compute correlation image from filtered movie

    corr_img = neighboring_pixel_correlation( filtered_movie );
    
    % dot correlation image with the tissue mask, since we assume 0
    % correlation for background
    corr_img_filt = corr_img .* tissue_mask;
    
    disp('Finished Calculating Correlation Image')
    %% Get ROIs from correlation image
    
    if param.manual_roi
        % manually define the ROIs
        roi_init = draw_rois( corr_img ); % no reason to do ROI selection for this
    else
        % try to programmatically define the ROIs
        roi_init = WatershedImage(corr_img); % psycho5 function that uses watershed, then tries to fill in the borders
    end
    
    roi_final = probe_correlation( filtered_movie, param, roi_init, corr_img );
    num_rois = max(roi_final,[],'all'); % recalculate number of ROIs
    
    disp('Finished Calculating ROIs')

    %% Get the response of each ROI over time 
    roi_dff = zeros( n_t, num_rois );

    % loop through each ROI
    for i = 1 : num_rois
        num_pixels = sum(roi_final == i, 'all'); % size of ROI
        
        % intensity of this ROI over time
        intensity_trace = (sum(sum((roi_final == i) .* filtered_movie, 1), 2)) ...
                              ./ (num_pixels);
                          
        % Delta F over F of this fly
        F0 = mean( intensity_trace( ismember(exp.epochVal, param.interleave_epochs) ) );
        roi_dff(:,i) = ( intensity_trace - F0 ) / F0;
    end
    
    epoch_trace = exp.epochVal; % rename this to make it easier to remember
    
end

%% Plot Stuff That Every User Probably Wants

% number of fly I'm on
fly_num = num2str( find(i_ex == analyze_these)  );

% plot the dx and dy
figure; hold on

subplot(2, 1, 1)
% subplot of dx change
plot(dx, 'linewidth', 1)
xlabel('time')
ylabel('dx (pixels)')
title([fly_num, '; Displacement in "x" direction'])
set(gca, 'fontsize', 20)

subplot(2, 1, 2)
% subplot of dy change
plot(dy, 'linewidth', 1)
xlabel('time')
ylabel('dy (pixels)')
title([fly_num, '; Displacement in "y" direction'])
set(gca, 'fontsize', 20)

% plot responses during 1st and 2nd probe
figure;
y_max = max(roi_dff( [param.corr_idxs{1}, param.corr_idxs{2}] , :), [], 'all') * 1.1; % 10% bigger than largest response
y_min = min(roi_dff( [param.corr_idxs{1}, param.corr_idxs{2}] , :), [], 'all');
y_min = y_min * (1 - sign(y_min)*0.1);
for i_p = 1 : 2

    % find where there's an epoch change
    epoch_change = find( diff( epoch_trace( param.corr_idxs{i_p} ) ) );

    % loop through both probes
    subplot(1,2,i_p)
    plot( roi_dff( param.corr_idxs{i_p} , :) ); hold on;
    for i_epoch = 1 : length(epoch_change)
        plot( ones(2,1).*epoch_change(i_epoch), [y_min; y_max], '--', 'color', 0.5 * ones(1,3) )
    end
    xlabel('time')
    ylabel('$\Delta F / F$', 'interpreter', 'latex')
    title(['Probe ', num2str(i_p)])
    subtitle('temporal trace used for ROI selection')
    set(gca, 'FontSize', 25)
    ylim( [ y_min, y_max] )
end

% plot total response
figure; hold on;
plot(roi_dff);
y_min = -0.5; y_max = 9.5;
epoch_change = find( diff( epoch_trace ) );
for i_epoch = 1 : length(epoch_change)
    plot( ones(2,1).*epoch_change(i_epoch), [y_min; y_max], '--', 'color', 0.5 * ones(1,3) )
end
xlabel('time');
ylabel('$\Delta F / F$', 'interpreter', 'latex')
title(['Fly ', fly_num ,' Response'])
ylim([y_min, y_max])
xlim([1, size(roi_dff,1)])
set(gca, 'FontSize', 25)

% plot the correlation image with labeled ROIs
figure;
imagesc(mean_movie); hold on;
colormap(gray)
axis off
for i_roi = 1 : num_rois
    visboundaries(bwboundaries(roi_final==i_roi), 'LineStyle', '--','LineWidth', 0.1, 'color', 'red');
end
title(['Fly ', fly_num ,' Mean Movie with ROIs'])
set(gca, 'FontSize', 25)

clc;
disp(['Finished Processing Fly ', fly_num, ' / ', num2str(length(analyze_these)) ]) 

