
function roi_mask = draw_rois( img )
% given an input image, this function returns the ROI mask. The user draws
% polygons around the ROIs they want. It does this until the user specifies
% otherwise. It also displays a black color for pixels already labeled by
% the user.

%figure; imagesc(img); colormap(gray)
num_rois = 0;
roi_mask = zeros( size(img) );
keep_going = true;

% force image to be between 0 and 1
img = img - min(img,[],'all');
img = img ./ max(img,[],'all');

while keep_going
    num_rois = num_rois + 1;
    BW = roipoly(img .* (roi_mask==0));
    close all;
    roi_mask( BW ) = num_rois;
    keep_going = input( 'Do you want to define another ROI (true/false)? ' );
end
        
end