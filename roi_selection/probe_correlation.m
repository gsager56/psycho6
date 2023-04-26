function roi_final = probe_correlation( movie, param, roi_init, corr_img, probe_idxs, tissue_mask)
% this roi selection function picks ROIs whose response is highly
% correlated from the first and second probe presentation

    probe_1_idxs = probe_idxs{1};
    probe_2_idxs = probe_idxs{2};
    roi_idx_init = unique( roi_init( roi_init>0 ) ); % list of unique, initial ROI

    % initialize the roi final mask
    roi_final = zeros( size(roi_init) );

    num_rois = 0; % keeps track of final ROI index
    c=0;
    corr_vec = zeros(length(unique(roi_init)), 1);
    for i_roi = 1 : length(roi_idx_init) % iterate through each roi
        % mask for this ROI
        this_roi_mask = roi_init == roi_idx_init(i_roi);

        if (sum(this_roi_mask .* corr_img, 'all') / sum(this_roi_mask,'all')) <= 0
            % this means the average correlation between pixels in this ROI
            % is less than 0, so there's no way this is a good ROI. For
            % the sake of time, let's not bother calculating the
            % correlation of this ROI's activity during the 1st and 2nd
            % probe
        elseif strcmp(param.probe_correlation_type, 'in tissue')
            % keep ROIs in the neural tissue
            frac_tissue = mean( tissue_mask(this_roi_mask) );
            if frac_tissue >= 0.5
                roi_size = sum(this_roi_mask, 'all'); % size of ROI in pixels

                % get intensity trace during 1st and second probe
                probe_1.trace = sum( sum( (movie(:,:,probe_1_idxs) .* this_roi_mask) ./ roi_size, 1), 2);
                probe_2.trace = sum( sum( (movie(:,:,probe_2_idxs) .* this_roi_mask) ./ roi_size, 1), 2);

                % calculate correlation coefficient
                correlation = corrcoef( probe_1.trace, probe_2.trace );

                % see if there's a significant correlation
                c=c+1;
                corr_vec(c) = correlation(1,2);

                if correlation(1,2) >= 0
                    % WE HAVE A GOOD ROI, YAYYYY!!
                    num_rois = num_rois + 1;
                    roi_final( this_roi_mask ) = num_rois;
                end
            end
            
        elseif strcmp(param.probe_correlation_type, 'probe indices')
            % use the indices specified by the user to compute which ROIs
            % are responding to the probe
            roi_size = sum(this_roi_mask, 'all'); % size of ROI in pixels

            % get intensity trace during 1st and second probe
            probe_1.trace = sum( sum( (movie(:,:,probe_1_idxs) .* this_roi_mask) ./ roi_size, 1), 2);
            probe_2.trace = sum( sum( (movie(:,:,probe_2_idxs) .* this_roi_mask) ./ roi_size, 1), 2);

            % calculate correlation coefficient
            correlation = corrcoef( probe_1.trace, probe_2.trace );

            % see if there's a significant correlation
            c=c+1;
            corr_vec(c) = correlation(1,2);

            if correlation(1,2) >= param.corr_thresh
                % WE HAVE A GOOD ROI, YAYYYY!!
                num_rois = num_rois + 1;
                roi_final( this_roi_mask ) = num_rois;
            end
        elseif strcmp(param.probe_correlation_type, 'epoch wise')
            % does ROI selection based on correlations of inidividual epoch
            % presentations. This is generally not recommended, but it is
            % sometimes useful
            
            % mean ROI intensity
            intensity_trace = sum( sum( (movie .* this_roi_mask) ./ sum(this_roi_mask, 'all'), 1), 2);
            intensity_trace = reshape(intensity_trace, [1, length(intensity_trace)] ); % make intensity trace 1D

            % convert to delta F over F where F0 is the average intensity over
            % the interleaves. This is a valid assumption for minimal
            % photobleaching
            F0 = mean( intensity_trace( ismember(exp.epochVal, param.interleave_epochs) ) );
            dff_trace = (intensity_trace - F0) ./ F0;

            % get correlation of responses with 1st probe and 2nd probe

            % loop through each probe epoch that you want to use for ROI
            % selection

            % find where last param file epoch is presented and use this to
            % find the 1st and 2nd probe presentation
            probe_cutoff = find(exp.epochVal == max(exp.epochVal), 1);

            % indices that correspond to the 1st probe
            probe_1.idx = 1 : find( ismember(exp.epochVal(1 : probe_cutoff), param.probe_epochs), 1, 'last' );

            % indices that correspond to the 2nd probe
            probe_2.idx = (find( ismember( exp.epochVal(probe_cutoff : end), param.probe_epochs ), 1) + probe_cutoff - 1) : n_t;

            probe_corr = zeros(1, length(param.epochs_for_selectivity)); % initialize this variable

            for i_epoch = 1 : length(param.epochs_for_selectivity)
                this_epoch = param.epochs_for_selectivity(i_epoch);
                probe_1.epoch = exp.epochVal(probe_1.idx) == this_epoch;

                probe_1.epoch_start = find(diff( probe_1.epoch ) == 1);
                probe_1.epoch_end = find(diff( probe_1.epoch ) == -1);

                if probe_1.epoch(1) == this_epoch
                    % using diff will miss the first idx
                    probe_1.epoch_start = [1; probe_1.epoch_start];
                elseif exp.epochVal(probe_1.idx(end)) == this_epoch
                    % using diff will miss the last idx
                    probe_1.epoch_end = [probe_1.epoch_end; probe_1.idx(end)];
                end

                if length(probe_1.epoch_start) ~= length(probe_1.epoch_end)
                    error('ERROR: ugh, dimensions mismatch')
                end

                % probe 2 trace
                probe_2.full_trace = dff_trace( find(exp.epochVal(probe_2.idx) == this_epoch) + probe_2.idx(1) - 1 );

                for i_p = 1 : length( probe_1.epoch_start )

                    probe_1.full_trace = dff_trace( probe_1.epoch_start(i_p) : probe_1.epoch_end(i_p) );

                    min_length = min( [length(probe_1.full_trace), length(probe_2.full_trace)] );

                    if length(probe_1.full_trace) == length(probe_2.full_trace)
                        % we're good
                        probe_1.trace = probe_1.full_trace;
                        probe_2.trace = probe_2.full_trace;
                    elseif (abs( length(probe_1.full_trace) - length(probe_2.full_trace) ) / min_length) < 0.1
                        % there's about a 5% difference in the length of the
                        % two vectors, so it's probably ok to just cutoff the
                        % longer trace
                        probe_1.trace = probe_1.full_trace(1 : min_length);
                        probe_2.trace = probe_2.full_trace(1 : min_length);
                    else
                        % there's a big difference in the dimensionality
                        error(['ERROR: this epoch of probe was presented for\n',...
                            ' significantly different time lengths in the beginning and end'])
                    end

                    correlation = corrcoef( probe_1.trace, probe_2.trace );
                    probe_corr(i_epoch) = probe_corr(i_epoch) + correlation(1,2);
                end
                probe_corr(i_epoch) = probe_corr(i_epoch) / length( probe_1.epoch_start ); % average correlation this ROI
            end

            if sum(probe_corr > param.corr_thresh) > param.num_corr_epochs
                % WE HAVE A GOOD ROI, YAYYYY!!
                num_rois = num_rois + 1;
                roi_final( this_roi_mask ) = num_rois;
            end
        end
    end
end