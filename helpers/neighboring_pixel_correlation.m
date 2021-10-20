function corr_img = neighboring_pixel_correlation( movie )
% given some input movie, this function finds the temporal correlation of
% each pixels with its respective neighboring pixels in the movie. These
% correlations are averaged over each possible neighbor. This function does
% consider that boundary pixels do not have 8 neighbors

    [n_rows, n_cols, n_t] = size(movie);
    
    % initialize the correlation image
    corr_img = zeros(n_rows, n_cols);

    % matrix indices of all you partners
    rows_j = {2:n_rows, 2:n_rows, 1:n_rows, 1:(n_rows-1), 1:(n_rows-1), 1:(n_rows-1), 1:n_rows, 2:n_rows};
    cols_j = {1:n_cols, 2:n_cols, 2:n_cols, 2:n_cols, 1:n_cols, 1:(n_cols-1), 1:(n_cols-1), 1:(n_cols-1)};
    rows_i = {1:(n_rows-1), 1:(n_rows-1), 1:n_rows, 2:n_rows, 2:n_rows, 2:n_rows, 1:n_rows, 1:(n_rows-1)};
    cols_i = {1:n_cols, 1:(n_cols-1), 1:(n_cols-1), 1:(n_cols-1), 1:n_cols, 2:n_cols, 2:n_cols, 2:n_cols};

    % Find correlation of each pixel with neighboring pixels
    for i = 1 : length(rows_j) % iterate through each of 8 neighbors
        center = movie(rows_i{i}, cols_i{i}, :);
        neigh = movie(rows_j{i}, cols_j{i}, :);

        % compute correlations
        r = ( n_t * sum( center.*neigh, 3) - sum(center,3).*sum(neigh,3) ) ./ ...
            sqrt( ( n_t*sum(center.^2,3) - sum(center,3).^2  ) .* ( n_t*sum(neigh.^2,3) - sum(neigh,3).^2 ) );

        % add to the correaltion image
        corr_img(rows_i{i}, cols_i{i}) = corr_img(rows_i{i}, cols_i{i}) + r;
    end

    % divide correlation image by number of neighboring pixels to get average
    % correlation
    corr_norm = zeros(size(corr_img)) + 8; % most pixels have eight neighbors

    % borders only have 5 neighbors
    corr_norm(:,1)=5; corr_norm(:,end)=5; corr_norm(1,:)=5; corr_norm(end,:)=5;

    % 4 edges only have three neighbors
    corr_norm(1,1)=3; corr_norm(1,end)=3; corr_norm(end,1)=3; corr_norm(end,end)=3;

    % divide correlation image by the respective normalization
    corr_img = corr_img ./ corr_norm ;
    
end