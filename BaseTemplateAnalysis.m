%% load in the experimental details
clear; clc;
% parameters about your experiment
cellType='Mi1';
stim='Flash_32frames';
sensor='GC6f'; 
flyEye='right';
surgeon='Harsh';

dataPath=GetPathsFromDatabase(cellType,stim,sensor,flyEye,surgeon);

%% PARAMETERS YOU NEED TO SPECIFY
param.stim = stim; % don't change this line
param.movieChan = 1; % which channels do you want to analyze in the movie
param.probe_epochs = 'do it for me';
param.interleave_epochs = 8;

% what flies do you want to analyze?
param.analyze_these = [3];
resp = cell( length(param.analyze_these), 1 );

% parameters for the ROI selection
param.corr_thresh = 0.9; % correlation threshold
param.probe_correlation = true; % do you want to use probe correlations for ROI selection
param.probe_correlation_type = 'probe indices'; % how to roi correlations from the 1st and 2nd probes
param.mean_amplification_thresh = 1.2; % this gets rid of ROIs with bad F0 fits when computing delta F over F

% parameters for what code will compute
param.force_alignment = false; % force alignment calculation?
param.force_roi_selection = true; % force calculation of ROIs?
param.manual_roi = false;
param.group_method = 'none'; % how to group selected ROIs
param.frac_interleave = 0.5; % fraction of the interleave to use for computing F0 for DF/F

%% Run Analysis
for i_ex = param.analyze_these
    
    %% load in experimental parameters and raw movie
    [exp_info, raw_movie, param] = get_exp_details( dataPath{i_ex}, param);
    % number of fly I'm on
    
    param.fly_num = find(i_ex == param.analyze_these);
    clc;
    
    % in case you want to plot the epochs...
    % figure; plot( exp_info.epochVal ); ylabel('epoch index'); xlabel('time index'); title(['Fly', num2str(param.fly_num)])

    %% register the relaxed mean image with each frame to get the alignment
    % find dimensions of image
    [n_rows, n_cols, n_t] = size(raw_movie);
    
    storage_folder = fullfile(dataPath{i_ex}, 'psycho6');
    medfilt_size = 5;
    
    if isfolder(storage_folder) && ~param.force_alignment
        % download what you need to do alignment from your last calculation
        load( fullfile( storage_folder, 'alignment_info.mat' ) );
        
        dx_filt = medfilt1(dx,medfilt_size);
        dy_filt = medfilt1(dy,medfilt_size);
        
        % only save the the rows and columns that stayed in frame
        good_cols = [abs(min(dx_filt)) + 1 : n_cols - max(dx_filt)]; % columns that stayed in frame
        good_rows = [max(dy_filt) + 1 : n_rows - abs(min(dy_filt))]; % rows that stayed in frame
        tissue_mask = tissue_mask(good_rows, good_cols);
        
        % initialize the aligned_movie
        aligned_movie = zeros( length(good_rows), length(good_cols), n_t);
        
        % iterate through time
        for i_t = 1 : n_t
            aligned_movie(:,:,i_t) = raw_movie(good_rows - dy_filt(i_t) , good_cols + dx_filt(i_t) , i_t);
        end
    else
        % do the alignment
        neuron_hist_thresh = 0.95;
        init_tissue_prob = linear_tissue_prob_HighNeuronThresh(mean(raw_movie,3), neuron_hist_thresh);
        %init_tissue_prob = linear_tissue_prob(mean(raw_movie,3)); % try this if not many pixels are responding
        init_tissue_prob = sqrt(init_tissue_prob); % I've found adding this nonlinearity usually helps out
        tissue_mask = ~RelaxationNetworkSimulation(init_tissue_prob);
        
        [dx, dy, aligned_movie, tissue_mask] = movie_alignment_registration( raw_movie, tissue_mask );
        
        % save the displacement in x and y
        mkdir(fullfile(dataPath{i_ex}, 'psycho6'));
        save(fullfile(dataPath{i_ex}, 'psycho6', 'alignment_info.mat'), 'dx', 'dy', 'tissue_mask');
        
        dx_filt = medfilt1(dx,medfilt_size);
        dy_filt = medfilt1(dy,medfilt_size);
        
        % only save the the rows and columns that stayed in frame
        good_cols = [abs(min(dx_filt)) + 1 : n_cols - max(dx_filt)]; % columns that stayed in frame
        good_rows = [max(dy_filt) + 1 : n_rows - abs(min(dy_filt))]; % rows that stayed in frame

        aligned_movie = aligned_movie(good_rows, good_cols, :);
        tissue_mask = tissue_mask(good_rows, good_cols);
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
    
    disp('Finished Calculating Correlation Image')
    %% Get ROIs from correlation image
    
    roi_mask_file = fullfile(storage_folder, 'roi_masks.mat');
    if ~param.force_roi_selection && isfile( roi_mask_file )
        load( roi_mask_file )
    else
        if param.manual_roi
            % manually define the ROIs
            roi_extract = draw_rois( corr_img );
        else
            % try to programmatically define the ROIs
            roi_extract = WatershedImage(corr_img); % psycho5 function that uses watershed, then tries to fill in the borders
        end
        probe_idxs = get_probe_idxs(exp_info.epochVal, param);
        roi_select = probe_correlation( filtered_movie, param, roi_extract, corr_img, probe_idxs );
        
        if strcmp(param.group_method, 'none')
            % do nothing
        elseif strcmp(param.group_method, 'touching')
            % group rois that are touching
            roi_select = bwlabel( roi_select>0 );
        elseif strcmp(param.group_method, 'manual')
            % allow user to manually group ROIs
            roi_select = manually_group_rois(roi_select, mean_movie);
        end
        save(roi_mask_file, 'roi_extract', 'roi_select');
    end

    %% Get the response of each ROI over time 
    
    %roi_dff_method = 'mean_resp'; % sets F0 to the mean interleave response
    roi_dff_method = 'exponential'; % fits exponential to each ROI
    [roi_dff, roi_final, mean_inter] = roi_dff_calc( param, exp_info, roi_select, filtered_movie, roi_dff_method);
    num_rois = max(roi_final,[],'all');
    epoch_trace = exp_info.epochVal; % rename this to make it easier to remember
    disp('Finished Calculating ROIs')
    %% save important output to resp cell array
    resp{param.fly_num}.dff = roi_dff; % delta F over F for each "good" ROI
    resp{param.fly_num}.epoch_trace = epoch_trace; % epochs presented over time
    resp{param.fly_num}.roi_final = roi_final; % final ROI mask
    
    %% Plot Stuff That Every User Probably Wants
    
    if num_rois > 0
    
        % plot the dx and dy
        MakeFigure; hold on
        subplot(2, 1, 1)
        % subplot of dx change
        plot(exp_info.time ./ 60, dx_filt, 'linewidth', 1)
        xlabel('time (minutes)')
        ylabel('dx (pixels)')
        title(['Fly ', num2str( param.fly_num ), '; Median Filtered Displacement in "x" direction'])
        set(gca, 'fontsize', 20)

        subplot(2, 1, 2)
        % subplot of dy change
        plot(exp_info.time ./ 60, dy_filt, 'linewidth', 1)
        xlabel('time (minutes)')
        ylabel('dy (pixels)')
        title(['Fly ', num2str( param.fly_num ), '; Median Filtered Displacement in "y" direction'])
        set(gca, 'fontsize', 20)

        % plot responses during 1st and 2nd probe
        MakeFigure;
        y_max = max(roi_dff( [probe_idxs{1}, probe_idxs{2}] , :), [], 'all') * 1.1; % 10% bigger than largest response
        y_min = min(roi_dff( [probe_idxs{1}, probe_idxs{2}] , :), [], 'all');
        y_min = y_min * (1 - sign(y_min)*0.1);
        for i_p = 1 : 2
            % find where there's an epoch change
            epoch_change = find( diff( epoch_trace( probe_idxs{i_p} ) ) );

            % loop through both probes
            subplot(1,2,i_p)
            plot( exp_info.time(probe_idxs{i_p}) ./ 60, roi_dff( probe_idxs{i_p} , :) ); hold on;
            for i_epoch = 1 : length(epoch_change)
                plot( ones(2,1).* exp_info.time( epoch_change(i_epoch) + probe_idxs{i_p}(1) )  ./ 60, [y_min; y_max], '--', 'color', 0.5 * ones(1,3) )
            end
            xlabel('time (minutes)')
            ylabel('$\Delta F / F$', 'interpreter', 'latex')
            title(['Fly ', num2str(param.fly_num), ' Probe ', num2str(i_p)])
            subtitle('temporal trace used for ROI selection')
            set(gca, 'FontSize', 25)
            ylim( [ y_min, y_max] )
            x_min = min(exp_info.time( probe_idxs{i_p} )  ./ 60);
            x_max = max(exp_info.time( probe_idxs{i_p} )  ./ 60);
            xlim([x_min, x_max])
        end

        % plot response over entire recording
        MakeFigure; hold on;
        plot(exp_info.time ./ 60, roi_dff);

        y_max = max(roi_dff, [], 'all') * 1.1; % 10% bigger than largest response
        y_min = min(roi_dff, [], 'all');
        y_min = y_min * (1 - sign(y_min)*0.1);

        epoch_change = find( diff( epoch_trace ) );
        for i_epoch = 1 : length(epoch_change)
            plot( ones(2,1).* exp_info.time( epoch_change(i_epoch) ) ./ 60, [y_min; y_max], '--', 'color', 0.5 * ones(1,3) )
        end
        xlabel('time (minutes)');
        ylabel('$\Delta F / F$', 'interpreter', 'latex')
        title(['Fly ', num2str( param.fly_num ) ,' Response'])
        ylim([y_min, y_max])
        xlim([0, exp_info.time(end)./60])
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
        
        % plot the region circling neural tissue
        MakeFigure;
        imagesc(mean_movie); hold on;
        colormap(gray)
        axis off
        visboundaries(bwboundaries(tissue_mask), 'LineStyle', '--','LineWidth', 0.1, 'color', 'red');
        title(['Fly ', num2str(param.fly_num) ,' Tissue Mask Over Mean Movie'])
        set(gca, 'FontSize', 25)
        
        % plot interleave fit divided by interleave response (mean)
        MakeFigure; hold on;
        plot(mean_inter.times / 60, mean_inter.fit ./ mean_inter.resp, 'LineWidth', 1)
        h = plot(mean_inter.times / 60, 1, 'LineWidth', 1);
        legend(h, 'Ideal Case')
        xlabel('Time (minutes)')
        ylabel('Factor Enhancement')
        title(['Fly ', num2str(param.fly_num) ,' Interleave Fit / Interleave Response'])
        set(gca, 'FontSize', 25)
    end
    
    clc;
    disp(['Finished Processing Fly ', num2str(param.fly_num), ' / ', num2str(length(param.analyze_these)) ]) 
end
