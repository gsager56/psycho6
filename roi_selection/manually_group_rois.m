function roi_select = manually_group_rois(roi_select_init, display_img)

%figure; imagesc(img); colormap(gray)
num_rois = 0;
roi_select = zeros( size(roi_select_init) );
keep_going = true;

% force image to be between 0 and 1
display_img = display_img - min(display_img,[],'all');
display_img = display_img ./ max(display_img,[],'all');

roi_display = roi_select_init;
num_rois = 0;
while keep_going
    figure;
    imagesc(display_img); hold on;
    colormap(gray)
    axis off
    all_rois = unique( roi_display(roi_display > 0) );
    
    for i_roi = 1 : length(all_rois)
        visboundaries(bwboundaries(roi_display==all_rois(i_roi)), 'LineStyle', '--','LineWidth', 0.1, 'color', 'red');
        hold on;
    end
    title('Draw Polygon Around ROIs You Want To Group')
    set(gca,'FontSize',20)
    BW = roipoly;
    close all;
    
    % confirm answer with user
    figure;
    chosen_roi_img = display_img;
    chosen_rois = unique(roi_display( BW ));
    chosen_rois = chosen_rois( chosen_rois>0 );
    for i_roi = 1 : length(chosen_rois)
        chosen_roi_img(roi_display==chosen_rois(i_roi)) = 1;
    end
    imagesc(chosen_roi_img); hold on;
    colormap(gray)
    axis off
    for i_roi = 1 : length(all_rois)
        visboundaries(bwboundaries(roi_display==all_rois(i_roi)), 'LineStyle', '--','LineWidth', 0.1, 'color', 'red');
        hold on;
    end
    title('Confirm Your Selection In The Command Window')
    set(gca,'FontSize',20)
    confirm_group = input( 'Would you like to group these ROIs (true/false)? ' );
    close all;
    
    if confirm_group
        num_rois = num_rois + 1;
        for i_roi = 1 : length(chosen_rois)
            mask = roi_display==chosen_rois(i_roi);
            roi_select(mask) = num_rois;
            roi_display(mask) = 0; % remove those ROIs from future consideration
        end
    end
    if any(roi_display > 0,'all')
        keep_going = input( 'Do you want to keep going (true/false)? ' );
    else
        disp('All ROIs have been grouped')
    end
end

if any(roi_display>0,'all')
    % there are still rois left
    unchosen_rois = unique(roi_display(roi_display>0));
    for i_roi = 1 : length(unchosen_rois)
        num_rois = num_rois + 1;
        roi_select(roi_display==unchosen_rois(i_roi)) = num_rois;
    end
end