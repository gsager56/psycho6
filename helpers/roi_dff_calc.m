function [roi_dff, roi_final, mean_inter] = roi_dff_calc( param, exp_info, roi_select, filtered_movie, roi_dff_method)
% param is the structure of params given by user about the experiment
% epoch_trace is the temporal trace of epoch indices
% roi_final is a mask of all the ROIs
% filtered_movie is the aligned and background subtracted movie

n_t = size(filtered_movie,3);
num_rois = max(roi_select,[],'all');

% intialize matrices
roi_dff = zeros( n_t, num_rois );
roi_raw_intensity = roi_dff;

% compute mean intensity temporal trace of each ROI
% loop through each ROI
for i_roi = 1 : num_rois
    num_pixels = sum(roi_select == i_roi, 'all'); % size of ROI

    % intensity of this ROI over time
    roi_raw_intensity(:,i_roi) = (sum(sum((roi_select == i_roi) .* filtered_movie, 1), 2)) ./ (num_pixels);
end

% boolean vector of time indices belonging to the interleave
init_inter_bool = ismember(exp_info.epochVal, param.interleave_epochs);
inter_bool = false( size(init_inter_bool) );
inter_epoch_ids = bwlabel(init_inter_bool);
for i_inter_epoch = 1 : max(inter_epoch_ids)
    this_epoch_bool = i_inter_epoch == inter_epoch_ids;
    idxs = find(this_epoch_bool, ceil(sum(this_epoch_bool)*param.frac_interleave), 'last');
    inter_bool(idxs) = true;
end

inter_times = exp_info.time( inter_bool )'; % time that interleave is presented
if strcmp(roi_dff_method, 'mean_resp')
    % calculate F0 as the mean response during the iterleaves. This method
    % assumes no photobleaching and that F0 is not changing in time. For
    % most experiments, these are probably OK assumptions. The simplicity
    % of this method is a major positive.
    
    F0 = mean( roi_raw_intensity(inter_bool, :), 1);
    roi_dff = (roi_raw_intensity - F0) ./ F0;
    inter_fits = repmat(F0,length(roi_raw_intensity),1);
elseif strcmp(roi_dff_method, 'exponential')
    % fit an exponential to the interleave responses
    % x(:,1) is the time vector
    % x(:,2) is the vector of responses in intensity
    exp_func = @(p, x) ( p(1).*exp(-x./p(2)) );
    inter_fits = zeros(length(inter_bool),num_rois);
    for i_roi = 1 : num_rois
        p0 = [mean(roi_raw_intensity(inter_bool,i_roi)), 60]; % initial parameter guess
        best_p = lsqcurvefit( exp_func, p0, inter_times, roi_raw_intensity(inter_bool,i_roi));
        clc;
        F0 = exp_func(best_p, exp_info.time )';
        inter_fits(:,i_roi) = F0;
        roi_dff(:,i_roi) = (roi_raw_intensity(:,i_roi) - F0) ./ F0;
    end
end

% Remove ROIs with poor interleave fits

% find mean interleave fit divided by mean interleave response for each
% interleave presentation
inter_epoch_ids = bwlabel(inter_bool);
mean_inter_resp = zeros(max(inter_epoch_ids),num_rois);
mean_inter_fit = zeros(max(inter_epoch_ids),num_rois);
mean_inter_times = zeros(max(inter_epoch_ids),1);
for i_inter_epoch = 1 : max(inter_epoch_ids)
    mean_inter_times(i_inter_epoch) = mean(exp_info.time(i_inter_epoch==inter_epoch_ids));
    mean_inter_resp(i_inter_epoch,:) = mean(roi_raw_intensity(i_inter_epoch==inter_epoch_ids,:),1);
    mean_inter_fit(i_inter_epoch,:) = mean(inter_fits(i_inter_epoch==inter_epoch_ids,:),1);
end

inter_factor = mean_inter_fit ./ mean_inter_resp;
inter_factor( inter_factor < 1 ) = 1 ./ inter_factor( inter_factor < 1 );
remove_rois_bool = mean(inter_factor,1) > param.mean_amplification_thresh;
remove_rois = find(remove_rois_bool);

if ~isempty(remove_rois)
    disp( ['Removing ', num2str(length(remove_rois)),'/', num2str(num_rois),' ROIs due to poor interleave fits'] )
    % there are ROIs outside the amplification threshold that need to
    % be removed
    roi_counter = 0;
    roi_final = zeros( size(roi_select) );
    for i_roi = 1 : num_rois
        if ~any(i_roi==remove_rois)
            % this ROI is not to be removed, so add it
            roi_counter = roi_counter + 1;
            roi_final( i_roi == roi_select ) = roi_counter;
        end
    end

    roi_dff = roi_dff(:,~remove_rois_bool);
    %roi_raw_intensity = roi_raw_intensity(:,~remove_rois_bool);
    %inter_fits = inter_fits(:,~remove_rois_bool);
    %inter_factor = inter_factor(:,~remove_rois_bool);
    mean_inter_fit = mean_inter_fit(:,~remove_rois_bool);
    mean_inter_resp = mean_inter_resp(:,~remove_rois_bool);

else
    % there are no ROIs to remove
    roi_final = roi_select;
end

mean_inter.resp = mean_inter_resp;
mean_inter.fit = mean_inter_fit;
mean_inter.times = mean_inter_times;

end