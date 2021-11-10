function roi_dff = roi_dff_calc( param, exp_info, roi_final, filtered_movie, roi_dff_method)
% param is the structure of params given by user about the experiment
% epoch_trace is the temporal trace of epoch indices
% roi_final is a mask of all the ROIs
% filtered_movie is the aligned and background subtracted movie

[n_rows, n_cols, n_t] = size(filtered_movie);
num_rois = max(roi_final,[],'all');

% intialize matrices
roi_dff = zeros( n_t, num_rois );
roi_trace = roi_dff;

% compute mean intensity temporal trace of each ROI
% loop through each ROI
for i_roi = 1 : num_rois
    num_pixels = sum(roi_final == i_roi, 'all'); % size of ROI

    % intensity of this ROI over time
    roi_trace(:,i_roi) = (sum(sum((roi_final == i_roi) .* filtered_movie, 1), 2)) ./ (num_pixels);
end

% boolean vector of time indices belonging to the interleave
inter_bool = ismember(exp_info.epochVal, param.interleave_epochs);
inter_times = exp_info.time( inter_bool )'; % time that interleave is presented

if strcmp(roi_dff_method, 'mean_resp')
    % calculate F0 as the mean response during the iterleaves. This method
    % assumes no photobleaching and that F0 is not changing in time. For
    % most experiments, these are probably OK assumptions. The simplicity
    % of this method is a major positive.
    
    F0 = mean( roi_trace(inter_bool, :), 1);
    roi_dff = (roi_trace - F0) ./ F0;
elseif strcmp(roi_dff_method, 'exponential')
    % fit an exponential to the interleave responses
    % x(:,1) is the time vector
    % x(:,2) is the vector of responses in intensity
    exp_func = @(p, x) ( p(1).*exp(-x./p(2)) );
    for i_roi = 1 : num_rois
        p0 = [mean(roi_trace(inter_bool,i_roi)), 60]; % initial parameter guess
        best_p = lsqcurvefit( exp_func, p0, inter_times, roi_trace(inter_bool,i_roi));
        clc;
        F0 = exp_func(best_p, exp_info.time )';
        roi_dff(:,i_roi) = (roi_trace(:,i_roi) - F0) ./ F0;
    end
end


end