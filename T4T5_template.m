%% PARAMETERS YOU NEED TO SPECIFY
clear; clc;
param.movieChan = 1; % which channels do you want to analyze in the movie
param.probe_epochs = 1 : 12;
param.interleave_epochs = 13;

% parameters for the ROi selection
param.corr_thresh = 0.2; % correlation threshold
param.probe_correlation = true; % do you want to use probe correlations for ROI selection
param.probe_correlation_type = 'hard code indices'; % type of correlation comparison for the 
param.corr_idxs{1} = [637, 1689, 634, 1689]; % indices of first probe I want to use for correlation
param.corr_idxs{2} = [5280, 6332, 5282, 6337]; % indices of 2nd probe I want to use for correlation

param.force_alignment = false; % force alignment calculation?
param.force_roi_selection = false; % force calculation of ROIs?
param.manual_roi = true; % manually define the ROIs

%% load in the experimental details

cellType='T4T5';
stim='MovingEdges_MovingSinusoid45d_MovingBars60d';
sensor='GACh';
flyEye='right';
surgeon='Annie';

dataPath=GetPathsFromDatabase(cellType,stim,sensor,flyEye,surgeon);

param.analyze_these = [1 : length(dataPath)];

%%
resp = cell( length(param.analyze_these), 1 );
for i_ex = param.analyze_these
    
    
    %% load in experimental parameters and raw movie
    [exp, raw_movie] = get_exp_details( dataPath{i_ex}, param);
    % number of fly I'm on
    param.fly_num = find(i_ex == param.analyze_these);
    clc;
    % in case you want to plot the epochs...
    % figure; plot( exp.epochVal ); ylabel('epoch index'); xlabel('time index'); title(['Fly', num2str(param.fly_num)])

    %% register the relaxed mean image with each frame to get the alignment
    % find dimensions of image
    [n_rows, n_cols, n_t] = size(raw_movie);
    
    storage_folder = fullfile(dataPath{i_ex}, 'psycho6');
    
    if isfolder(storage_folder) && ~param.force_alignment
        % download what you need to do alignment from your last calculation
        load( fullfile( storage_folder, 'alignment_info.mat' ) );
        
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
    
    roi_mask_file = fullfile(storage_folder, 'roi_masks.mat');
    if ~param.force_roi_selection && isfile( roi_mask_file )
        load( roi_mask_file )
    else
        if param.manual_roi
            % manually define the ROIs
            roi_init = draw_rois( corr_img ); % no reason to do ROI selection for this
        else
            % try to programmatically define the ROIs
            roi_init = WatershedImage(corr_img); % psycho5 function that uses watershed, then tries to fill in the borders
        end
        roi_final = probe_correlation( filtered_movie, param, roi_init, corr_img );
        
        save(roi_mask_file, 'roi_init', 'roi_final');
    end
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
    
    %% Plot Stuff That Every User Probably Wants
    
    if num_rois > 0
    
        % plot the dx and dy
        MakeFigure; hold on

        subplot(2, 1, 1)
        % subplot of dx change
        plot(exp.time ./ 60, dx, 'linewidth', 1)
        xlabel('time (minutes)')
        ylabel('dx (pixels)')
        title([num2str( param.fly_num ), '; Displacement in "x" direction'])
        set(gca, 'fontsize', 20)

        subplot(2, 1, 2)
        % subplot of dy change
        plot(exp.time ./ 60, dy, 'linewidth', 1)
        xlabel('time (minutes)')
        ylabel('dy (pixels)')
        title([num2str( param.fly_num ), '; Displacement in "y" direction'])
        set(gca, 'fontsize', 20)

        % plot responses during 1st and 2nd probe
        MakeFigure;
        probe_idxs{1} = param.corr_idxs{1}( 1 + 2*(param.fly_num-1) ) : param.corr_idxs{1}( 2 + 2*(param.fly_num-1) );
        probe_idxs{2} = param.corr_idxs{2}( 1 + 2*(param.fly_num-1) ) : param.corr_idxs{2}( 2 + 2*(param.fly_num-1) );
        y_max = max(roi_dff( [probe_idxs{1}, probe_idxs{2}] , :), [], 'all') * 1.1; % 10% bigger than largest response
        y_min = min(roi_dff( [probe_idxs{1}, probe_idxs{2}] , :), [], 'all');
        y_min = y_min * (1 - sign(y_min)*0.1);
        for i_p = 1 : 2

            % find where there's an epoch change
            epoch_change = find( diff( epoch_trace( probe_idxs{i_p} ) ) );

            % loop through both probes
            subplot(1,2,i_p)
            plot( exp.time(probe_idxs{i_p}) ./ 60, roi_dff( probe_idxs{i_p} , :) ); hold on;
            for i_epoch = 1 : length(epoch_change)
                plot( ones(2,1).* exp.time( epoch_change(i_epoch) + probe_idxs{i_p}(1) )  ./ 60, [y_min; y_max], '--', 'color', 0.5 * ones(1,3) )
            end
            xlabel('time (minutes)')
            ylabel('$\Delta F / F$', 'interpreter', 'latex')
            title(['Probe ', num2str(i_p)])
            subtitle('temporal trace used for ROI selection')
            set(gca, 'FontSize', 25)
            ylim( [ y_min, y_max] )
            x_min = min(exp.time( probe_idxs{i_p} )  ./ 60);
            x_max = max(exp.time( probe_idxs{i_p} )  ./ 60);
            xlim([x_min, x_max])
        end

        % plot total response
        MakeFigure; hold on;
        plot(exp.time ./ 60, roi_dff);

        y_max = max(roi_dff, [], 'all') * 1.1; % 10% bigger than largest response
        y_min = min(roi_dff, [], 'all');
        y_min = y_min * (1 - sign(y_min)*0.1);

        epoch_change = find( diff( epoch_trace ) );
        for i_epoch = 1 : length(epoch_change)
            plot( ones(2,1).* exp.time( epoch_change(i_epoch) ) ./ 60, [y_min; y_max], '--', 'color', 0.5 * ones(1,3) )
        end
        xlabel('time (minutes)');
        ylabel('$\Delta F / F$', 'interpreter', 'latex')
        title(['Fly ', num2str( param.fly_num ) ,' Response'])
        ylim([y_min, y_max])
        xlim([0, exp.time(end)./60])
        set(gca, 'FontSize', 25)

        % plot the correlation image with labeled ROIs
        MakeFigure;
        imagesc(corr_img); hold on;
        caxis([0,1])
        colorbar
        colormap(gray)
        axis off
        for i_roi = 1 : num_rois
            visboundaries(bwboundaries(roi_final==i_roi), 'LineStyle', '--','LineWidth', 0.1, 'color', 'red');
        end
        title(['Fly ', num2str(param.fly_num) ,' Correlation Image with ROIs'])
        set(gca, 'FontSize', 25)
        
        % plot the mean movie with labeled ROIs
        MakeFigure;
        imagesc(mean_movie); hold on;
        colormap(gray)
        axis off
        for i_roi = 1 : num_rois
            visboundaries(bwboundaries(roi_final==i_roi), 'LineStyle', '--','LineWidth', 0.1, 'color', 'red');
        end
        title(['Fly ', num2str(param.fly_num) ,' Mean Movie with ROIs'])
        set(gca, 'FontSize', 25)
    end
    
    clc;
    disp(['Finished Processing Fly ', num2str(param.fly_num), ' / ', num2str(length(param.analyze_these)) ]) 
end

%% loop through each recording to be analyzed

% for i_ex = analyze_these
%     
%     %% load in experimental parameters and raw movie
%     [exp, raw_movie] = get_exp_details( dataPath{analyze_these(i_ex)}, param);
%     
%     %% Compute the correlation image during the stimulus
%     % find dimensions of image
%     [n_rows, n_cols, n_t] = size(raw_movie);
% 
%     % initialize the correlation image
%     corr_img = zeros(n_rows, n_cols);
% 
%     % matrix indices of all you partners
%     rows_j = {2:n_rows, 2:n_rows, 1:n_rows, 1:(n_rows-1), 1:(n_rows-1), 1:(n_rows-1), 1:n_rows, 2:n_rows};
%     cols_j = {1:n_cols, 2:n_cols, 2:n_cols, 2:n_cols, 1:n_cols, 1:(n_cols-1), 1:(n_cols-1), 1:(n_cols-1)};
%     rows_i = {1:(n_rows-1), 1:(n_rows-1), 1:n_rows, 2:n_rows, 2:n_rows, 2:n_rows, 1:n_rows, 1:(n_rows-1)};
%     cols_i = {1:n_cols, 1:(n_cols-1), 1:(n_cols-1), 1:(n_cols-1), 1:n_cols, 2:n_cols, 2:n_cols, 2:n_cols};
% 
%     % Find correlation of each pixel with neighboring pixels
%     for i = 1 : length(rows_j) % iterate through each of 8 neighbors
%         center = raw_movie(rows_i{i}, cols_i{i}, :);
%         neigh = raw_movie(rows_j{i}, cols_j{i}, :);
% 
%         % compute correlations
%         r = ( n_t * sum( center.*neigh, 3) - sum(center,3).*sum(neigh,3) ) ./ ...
%             sqrt( ( n_t*sum(center.^2,3) - sum(center,3).^2  ) .* ( n_t*sum(neigh.^2,3) - sum(neigh,3).^2 ) );
% 
%         % add to the correaltion image
%         corr_img(rows_i{i}, cols_i{i}) = corr_img(rows_i{i}, cols_i{i}) + r;
%     end
% 
%     % divide correlation image by number of neighboring pixels to get average
%     % correlation
%     corr_norm = zeros(size(corr_img)) + 8; % most pixels have eight neighbors
% 
%     % borders only have 5 neighbors
%     corr_norm(:,1)=5; corr_norm(:,end)=5; corr_norm(1,:)=5; corr_norm(end,:)=5;
% 
%     % 4 edges only have three neighbors
%     corr_norm(1,1)=3; corr_norm(1,end)=3; corr_norm(end,1)=3; corr_norm(end,end)=3;
% 
%     % divide correlation image by the respective normalization and make
%     % nonnegative
%     corr_img = abs( corr_img ./ corr_norm );
%     
%     disp('Finished Calculating Correlation Image')
% 
%     %% register the correlation image with each frame to get the alignment
%     
%     alignment_folder = fullfile(dataPath{i_ex}, 'concise_psycho5');
%     
%     if isfolder(alignment_folder)
%         % download what you need to do alignment from your last calculation
%         load( fullfile( alignment_folder, 'alignment_info.mat' ) );
%         
%         % only save the the rows and columns that stayed in frame
%         good_cols = [abs(min(dx)) + 1 : n_cols - max(dx)]; % columns that stayed in frame
%         good_rows = [max(dy) + 1 : n_rows - abs(min(dy))]; % rows that stayed in frame
%         
%         % initialize the aligned_movie
%         aligned_movie = zeros( length(good_rows), length(good_cols), n_t);
%         
%         % iterate through time
%         for i_t = 1 : n_t
%             cols = good_cols + dx(i_t);
%             rows = good_rows - dy(i_t);
%             
%             aligned_movie(:,:,i_t) = raw_movie(rows, cols, i_t);
%         end
%     else
%         % do the alignment
%         disp('Doing alignment, this will take on the order of 10s of minutes')
%         
%         % get the relaxed image
%         tissue_mask = ~RelaxationNetwork(corr_img, 5, 'correlation_img');
% 
%         dx = zeros(1,n_t);
%         dy = zeros(1,n_t);
% 
%         transformType = 'translation';
%         optimizer = registration.optimizer.RegularStepGradientDescent;
%         metric = registration.metric.MattesMutualInformation;
% 
%         aligned_movie = zeros(size(raw_movie));
% 
%         median_filt_raw_movie = medfilt3( raw_movie ); % median filter raw movie to help with identification
% 
%         for i = 1 : n_t % iterate through each stimulus frame
% 
%             moving(:,:) = median_filt_raw_movie(:,:,i); % image to be registered
%             aligned_movie(:,:,i) = imregister(moving,tissue_mask * max(moving,[],'all'),transformType,optimizer,metric);
% 
%             % compute displacement in x
%             if all(aligned_movie(:,1,i)==0)
%                 % the structure moved to the left
%                 this_dx = 0;
%                 while all(aligned_movie(:,this_dx+1,i)==0)
%                     this_dx = this_dx + 1; % this counts the displacement in x
%                 end
% 
%                 this_dx = -this_dx; % define left as negative
%             elseif all(aligned_movie(:,end,i)==0)
%                 % this structure moved to the right
%                 this_dx = n_cols;
%                 while all(aligned_movie(:,this_dx-1,i)==0)
%                     this_dx = this_dx - 1; % this counts the displacement in x
%                 end
% 
%                 this_dx = n_cols - this_dx;
%             else
%                 % there was no movement
%                 this_dx = 0;
%             end
% 
%             % compute displacement in y
%             if all(aligned_movie(1,:,i)==0)
%                 % the structure moved up (assuming the 1st row is the top)
%                 this_dy = 0;
%                 while all(aligned_movie(this_dy+1,:,i)==0)
%                     this_dy = this_dy + 1; % this counts the displacement in x
%                 end
% 
%             elseif all(aligned_movie(end,:,i)==0)
%                 % the structure moved down (assuming the 1st row is the top)
%                 this_dy = n_rows;
%                 while all(aligned_movie(this_dy-1,:,i)==0)
%                     this_dy = this_dy - 1; % this counts the displacement in x
%                 end
% 
%                 this_dy = -(n_rows - this_dy); % define down as negative
%             else
%                 % there was no movement
%                 this_dy = 0;
%             end
% 
%             % update dx and dy vectors
%             dx(i) = this_dx;
%             dy(i) = this_dy;
% 
%             if rem(i,round(n_t/10))==0
%                 disp([num2str(round(100 * i / n_t)), '% done with alignment'])
%             end
%             
%             
%         end
% 
%         % save the displacement in x and y
%         mkdir(fullfile(dataPath{i_ex}, 'concise_psycho5'));
%         save(fullfile(dataPath{i_ex}, 'concise_psycho5', 'alignment_info.mat'), ...
%             'dx', 'dy', 'tissue_mask');
%         
%         % only save the the rows and columns that stayed in frame
%         good_cols = [abs(min(dx)) + 1 : n_cols - max(dx)]; % columns that stayed in frame
%         good_rows = [max(dy) + 1 : n_rows - abs(min(dy))]; % rows that stayed in frame
% 
%         aligned_movie = aligned_movie(good_rows, good_cols, :);
%     end
%     
%     corr_img = corr_img(good_rows, good_cols);
%     
%     if size(tissue_mask,1) ~= length(good_rows)
%         tissue_mask = tissue_mask(good_rows, good_cols);
%     end
%     
%     %% Subtract out background in aligned_movie
%     
%     % initiliaze the filtered movie
%     filtered_movie = zeros(size(aligned_movie));
%     
%     % define background as below median correlation of the 0 pixels in
%     % tissue_mask
%     bkg_mask = (corr_img < median(corr_img(~tissue_mask))) & ~tissue_mask;
%     
%     % loop through time
%     for i = 1 : n_t
%         this_frame = aligned_movie(:,:,i);
%         filtered_movie(:,:,i) = aligned_movie(:,:,i) - mean(this_frame(bkg_mask));
%     end
%     
%     disp('Finished Processing Raw Movie')
%     
%     %% Get ROIs from correlation image
%     roi_init = WatershedImage(corr_img); % psycho5 function that uses watershed, then tries to fill in the borders
%     roi_idx_init = unique(roi_init); % list of unique, initial ROI
% 
%     % initialize the roi final mask
%     roi_final = zeros( size(roi_init) );
% 
%     num_rois = 0; % keeps track of final ROI index
%     for i_roi = 1 : length(roi_idx_init) % iterate through each roi
%         % mask for this ROI
%         this_roi_mask = roi_init == roi_idx_init(i_roi);
%         
%         if (sum(this_roi_mask .* corr_img, 'all') / sum(this_roi_mask,'all')) < 0.1
%             % this means the average correlation between pixels in this ROI
%             % is less than 0.1, so there's no way this is a good ROI. For
%             % the sake of time, let's not bother calculating the
%             % correlation of this ROI's activity during the 1st and 2nd
%             % probe
%         else
%         
%             % mean ROI intensity
%             intensity_trace = sum( sum( (filtered_movie .* this_roi_mask) ./ sum(this_roi_mask, 'all'), 1), 2);
%             intensity_trace = reshape(intensity_trace, [1, length(intensity_trace)] ); % make intensity trace 1D
% 
%             % convert to delta F over F where F0 is the average intensity over
%             % the interleaves. This is a valid assumption for minimal
%             % photobleaching
%             F0 = mean( intensity_trace( ismember(exp.epochVal, param.interleave_epochs) ) );
%             dff_trace = (intensity_trace - F0) ./ F0;
% 
% 
%             % get correlation of responses with 1st probe and 2nd probe
% 
%             % loop through each probe epoch that you want to use for ROI
%             % selection
% 
%             % find where last param file epoch is presented and use this to
%             % find the 1st and 2nd probe presentation
%             probe_cutoff = find(exp.epochVal == max(exp.epochVal), 1);
% 
%             % indices that correspond to the 1st probe
%             probe_1.idx = 1 : find( ismember(exp.epochVal(1 : probe_cutoff), param.probe_epochs), 1, 'last' );
% 
%             % indices that correspond to the 2nd probe
%             probe_2.idx = (find( ismember( exp.epochVal(probe_cutoff : end), param.probe_epochs ), 1) + probe_cutoff - 1) : n_t;
%             
%             probe_corr = zeros(1, length(param.epochs_for_selectivity)); % initialize this variable
% 
%             for i_epoch = 1 : length(param.epochs_for_selectivity)
%                 this_epoch = param.epochs_for_selectivity(i_epoch);
%                 probe_1.epoch = exp.epochVal(probe_1.idx) == this_epoch;
% 
%                 probe_1.epoch_start = find(diff( probe_1.epoch ) == 1);
%                 probe_1.epoch_end = find(diff( probe_1.epoch ) == -1);
% 
%                 if probe_1.epoch(1) == this_epoch
%                     % using diff will miss the first idx
%                     probe_1.epoch_start = [1; probe_1.epoch_start];
%                 elseif exp.epochVal(probe_1.idx(end)) == this_epoch
%                     % using diff will miss the last idx
%                     probe_1.epoch_end = [probe_1.epoch_end; probe_1.idx(end)];
%                 end
% 
%                 if length(probe_1.epoch_start) ~= length(probe_1.epoch_end)
%                     error('ERROR: ugh, dimensions mismatch')
%                 end
% 
%                 % probe 2 trace
%                 probe_2.full_trace = dff_trace( find(exp.epochVal(probe_2.idx) == this_epoch) + probe_2.idx(1) - 1 );
%                 
%                 
% 
%                 for i_p = 1 : length( probe_1.epoch_start )
% 
%                     probe_1.full_trace = dff_trace( probe_1.epoch_start(i_p) : probe_1.epoch_end(i_p) );
% 
%                     min_length = min( [length(probe_1.full_trace), length(probe_2.full_trace)] );
% 
%                     if length(probe_1.full_trace) == length(probe_2.full_trace)
%                         % we're good
%                         probe_1.trace = probe_1.full_trace;
%                         probe_2.trace = probe_2.full_trace;
%                     elseif (abs( length(probe_1.full_trace) - length(probe_2.full_trace) ) / min_length) < 0.1
%                         % there's about a 5% difference in the length of the
%                         % two vectors, so it's probably ok to just cutoff the
%                         % longer trace
%                         probe_1.trace = probe_1.full_trace(1 : min_length);
%                         probe_2.trace = probe_2.full_trace(1 : min_length);
%                     else
%                         % there's a big difference in the dimensionality
%                         error(['ERROR: this epoch of probe was presented for\n',...
%                             ' significantly different time lengths in the beginning and end'])
%                     end
% 
%                     correlation = corrcoef( probe_1.trace, probe_2.trace );
%                     probe_corr(i_epoch) = probe_corr(i_epoch) + correlation(1,2);
%                 end
%                 probe_corr(i_epoch) = probe_corr(i_epoch) / length( probe_1.epoch_start ); % average correlation this ROI
%             end
%             
%             if sum(probe_corr > param.corr_thresh) > param.num_corr_epochs
%                 % WE HAVE A GOOD ROI, YAYYYY!!
%                 num_rois = num_rois + 1;
%                 
%                 roi_final( this_roi_mask ) = num_rois;
%             end
%         end
%     end
%     
%     disp('Finished Calculating ROIs')
% 
%     %% Get the response of each ROI over time 
%     roi_dff = zeros( n_t, num_rois );
% 
%     % loop through each ROI
%     for i = 1 : num_rois
%         num_pixels = sum(roi_final == i, 'all'); % size of ROI
%         
%         % intensity of this ROI over time
%         intensity_trace = (sum(sum((roi_final == i) .* filtered_movie, 1), 2)) ...
%                               ./ (num_pixels);
%                           
%         % Delta F over F of this fly
%         F0 = mean( intensity_trace( ismember(exp.epochVal, param.interleave_epochs) ) );
%         roi_dff(:,i) = ( intensity_trace - F0 ) / F0;
%     end
%     
%     %% Plot Stuff That Every User Probably Wants
%     
%     % plot the dx and dy
%     figure; hold on
% 
%     subplot(2, 1, 1)
%     % subplot of dx change
%     plot(dx, 'linewidth', 1)
%     xlabel('time')
%     ylabel('dx (pixels)')
%     title('Displacement in "x" direction')
%     set(gca, 'fontsize', 20)
% 
%     subplot(2, 1, 2)
%     % subplot of dy change
%     plot(dy, 'linewidth', 1)
%     xlabel('time')
%     ylabel('dy (pixels)')
%     title('Displacement in "y" direction')
%     set(gca, 'fontsize', 20)
%     
%     % plot total response
%     figure;
%     plot(roi_dff);
%     xlabel('time');
%     ylabel('$\Delta F / F$', 'interpreter', 'latex')
%     title(['Fly ', num2str(i_ex) ,' Response'])
%     set(gca, 'FontSize', 25)
%     
%     % plot the correlation image with labeled ROIs
%     figure;
%     imagesc(corr_img); hold on;
%     colormap(gray)
%     axis off
%     visboundaries(bwboundaries(roi_final>0), 'LineStyle', '--','LineWidth', 0.1, 'color', 'red');
%     title(['Fly ', num2str(i_ex) ,' Correlation Image with ROIs'])
%     set(gca, 'FontSize', 25)
%     
%     clc;
%     disp(['Finished Processing Fly ', num2str( find(i_ex == analyze_these)  ), ' / ', num2str(length(analyze_these)) ]) 
% end


