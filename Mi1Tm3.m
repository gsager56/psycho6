clear all;
%% PARAMETERS YOU NEED TO SPECIFY
clear; clc;
param.movieChan = [1,2]; % for two channels
param.probe_epochs = 1 : 2;
param.interleave_epochs = 3; 
%% load in the experimental details
cellType='Mi1Tm3';
stim='HarshFlash32Frames180framesinterleave';
sensor='GC6f;jRGECO1a'; 
flyEye='right'; 
surgeon='Harsh';
dataPath=GetPathsFromDatabase(cellType,stim,sensor,flyEye,surgeon);
param.analyze_these = [25,28];
%%
% parameters for the ROi selection
param.corr_thresh = 0.2; % correlation threshold1
param.probe_correlation = true; % do you want to use probe correlations for ROI selection
param.probe_correlation_type = 'hard code indices'; % type of correlation comparison for the 
param.corr_idxs{1} = [611 789]; % indices of first probe I want to use for correlation
param.corr_idxs{2} = [5716 5894]; % indices of 2nd probe I want to use for correlation
param.force_alignment = false; % force alignment calculation?yy
param.force_roi_selection =true ; % 
param.manual_roi = true; % manually define the ROIs 

%% 
resp = cell( length(param.analyze_these), 1 );

for i_ex = param.analyze_these
    %% load in experimental parameters and raw movie
    [expt, raw_movie_all] = get_exp_details( dataPath{i_ex}, param);
 
    %figure; plot( expt.epochVal ); ylabel('epoch index'); xlabel('time index'); title(['Fly', num2str(param.fly_num)])
    
    % number of fly I'm on
    param.fly_num = find(i_ex == param.analyze_these);
    clc;
    %%
    for i_channel = 1: size(raw_movie_all,4)
        raw_movie = raw_movie_all(:,:,:,i_channel);

        %% register the relaxed mean image with each frame to get the alignment
        % find dimensions of image
        [n_rows, n_cols, n_t] = size(raw_movie);
        storage_folder = fullfile(dataPath{i_ex},'psycho6',num2str(i_channel));

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
            [dx, dy, aligned_movie, tissue_mask] =  movie_alignment_registration_Two_channels(raw_movie, dataPath{i_ex},i_channel );

            % only save the the rows and columns that stayed in frame
            good_cols = [abs(min(dx)) + 1 : n_cols - max(dx)]; % columns that stayed in frame
            good_rows = [max(dy) + 1 : n_rows - abs(min(dy))]; % rows that stayed in frame

            aligned_movie = aligned_movie(good_rows, good_cols, :);
        end

        % redefine n_rows and n_cols to be size of cropped movie
        [n_rows, n_cols, ~] = size(aligned_movie);
        disp('Finished Processing Raw Movie')
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
        disp(strcat('Total ROIs found:',num2str(num_rois)));
        

        %% Get the response of each ROI over time
        roi_dff{i_channel} = zeros( n_t, num_rois );
        
        %         for i=1:num_rois
        % %             num_pixels = sum(roi_final == i, 'all'); % size of ROI
        % %            % intensity of this ROI over time
        % %             intensity_trace = (sum(sum((roi_final== i) .* filtered_movie, 1), 2)) ...
        % %                 ./ (num_pixels);
        % %
        %         end
%% Get the response of each ROI over time
% roi_dff_method = 'mean_resp'; % sets F0 to the mean interleave response
roi_dff_method = 'exponential'; % fits exponential to each ROI
roi_dff{i_channel} = roi_dff_calc( param, expt, roi_final, filtered_movie, roi_dff_method);
epoch_trace = expt.epochVal; % rename this to make it easier to remember
time(:,1)=expt.time;

        %% save important output to resp cell array
        resp{param.fly_num}.dff{i_channel} = roi_dff{i_channel};
        resp{param.fly_num}.epoch_trace = epoch_trace;
        resp{param.fly_num}.time = time;
        resp{param.fly_num}.mean_movie{i_channel} = mean_movie;
        resp{param.fly_num}.roi_final{i_channel} =roi_final;


        %% Plot Stuff That Every User Probably Wants

        if num_rois > 0

            % plot the dx and dy
            MakeFigure; hold on

            subplot(2, 1, 1)
            % subplot of dx change
            plot(expt.time ./ 60, dx, 'linewidth', 1)
            xlabel('time (minutes)')
            ylabel('dx (pixels)')
            title([num2str( param.fly_num ), '; Displacement in "x" direction'])
            set(gca, 'fontsize', 20)

            subplot(2, 1, 2)
            % subplot of dy change
            plot(expt.time ./ 60, dy, 'linewidth', 1)
            xlabel('time (minutes)')
            ylabel('dy (pixels)')
            title([num2str( param.fly_num ), '; Displacement in "y" direction'])
            set(gca, 'fontsize', 20)

            % plot responses during 1st and 2nd probe
            MakeFigure;
            probe_idxs{1} = param.corr_idxs{1}( 1 + 2*(param.fly_num-1) ) : param.corr_idxs{1}( 2 + 2*(param.fly_num-1) );
            probe_idxs{2} = param.corr_idxs{2}( 1 + 2*(param.fly_num-1) ) : param.corr_idxs{2}( 2 + 2*(param.fly_num-1) );
            y_max = max(roi_dff{i_channel}( [probe_idxs{1}, probe_idxs{2}] , :), [], 'all') * 1.1; % 10% bigger than largest response
            y_min = min(roi_dff{i_channel}( [probe_idxs{1}, probe_idxs{2}] , :), [], 'all');
            y_min = y_min * (1 - sign(y_min)*0.1);
            for i_p = 1 : 2

                % find where there's an epoch change
                epoch_change = find( diff( epoch_trace( probe_idxs{i_p} ) ) );

                % loop through both probes
                subplot(1,2,i_p)
                plot( expt.time(probe_idxs{i_p}) ./ 60, roi_dff{i_channel}( probe_idxs{i_p} , :) ); hold on;
                for i_epoch = 1 : length(epoch_change)
                    plot( ones(2,1).* expt.time( epoch_change(i_epoch) + probe_idxs{i_p}(1) )  ./ 60, [y_min; y_max], '--', 'color', 0.5 * ones(1,3) )
                end
                xlabel('time (minutes)')
                ylabel('$\Delta F / F$', 'interpreter', 'latex')
                title(['Probe ', num2str(i_p)])
                subtitle('temporal trace used for ROI selection')
                set(gca, 'FontSize', 25)
                ylim( [ y_min, y_max] )
                x_min = min(expt.time( probe_idxs{i_p} )  ./ 60);
                x_max = max(expt.time( probe_idxs{i_p} )  ./ 60);
                xlim([x_min, x_max])
            end

            % plot total response
            MakeFigure; hold on;
            plot(expt.time ./ 60, roi_dff{i_channel});

            y_max = max(roi_dff{i_channel}, [], 'all') * 1.1; % 10% bigger than largest response
            y_min = min(roi_dff{i_channel}, [], 'all');
            y_min = y_min * (1 - sign(y_min)*0.1);

            epoch_change = find( diff( epoch_trace ) );
            for i_epoch = 1 : length(epoch_change)
                plot( ones(2,1).* expt.time( epoch_change(i_epoch) ) ./ 60, [y_min; y_max], '--', 'color', 0.5 * ones(1,3) )
            end
            xlabel('time (minutes)');
            ylabel('$\Delta F / F$', 'interpreter', 'latex')
            title(['Fly ', num2str( param.fly_num ) ,' Response'])
            ylim([y_min, y_max])
            xlim([0, expt.time(end)./60])
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

        clearvars -except resp param expt raw_movie_all i_ex dataPath
    end
    disp(['Finished Processing Fly ', num2str(param.fly_num), ' / ', num2str(length(param.analyze_these)) ]) 
end

selpath = uigetdir;
name='resp.mat';
filename=strcat(selpath,'/fly_num',num2str(param.fly_num),name);
save(filename, 'resp')
% 
% name='stiminfo.mat';
% save(filename, 'resp')









