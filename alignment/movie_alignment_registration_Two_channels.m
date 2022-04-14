    
function [dx, dy, aligned_movie, tissue_mask] = movie_alignment_registration_Two_channels( movie, dataPath, channel )
    % do the alignment
    disp('Doing alignment, this will take on the order of 10s of minutes')
    
    [n_rows, n_cols, n_t] = size( movie );
    
    % get the relaxed image
    mean_movie = mean( movie, 3);
    tissue_mask = ~RelaxationNetwork(mean_movie, 5, 'mean_movie');
    
    dx = zeros(1,n_t);
    dy = zeros(1,n_t);
    
    transformType = 'translation';
      optimizer = registration.optimizer.RegularStepGradientDescent;
    metric = registration.metric.MattesMutualInformation;
    
    aligned_movie = zeros(size(movie));
    
    median_filt_movie = medfilt3( movie ); % median filter raw movie to help with identification
    
    for i_t = 1 : n_t % iterate through each stimulus frame
        
        moving(:,:) = median_filt_movie(:,:,i_t); % image to be registered
        aligned_movie(:,:,i_t) = imregister(moving,tissue_mask * max(moving,[],'all'),transformType,optimizer,metric);
        
        % compute displacement in x
        if all(aligned_movie(:,1,i_t)==0)
            % the structure moved to the left
            this_dx = 0;
            while all(aligned_movie(:,this_dx+1,i_t)==0)
                this_dx = this_dx + 1; % this counts the displacement in x
            end

            this_dx = -this_dx; % define left as negative
        elseif all(aligned_movie(:,end,i_t)==0)
            % this structure moved to the right
            this_dx = n_cols;
            while all(aligned_movie(:,this_dx-1,i_t)==0)
                this_dx = this_dx - 1; % this counts the displacement in x
            end

            this_dx = n_cols - this_dx;
        else
            % there was no movement
            this_dx = 0;
        end

        % compute displacement in y
        if all(aligned_movie(1,:,i_t)==0)
            % the structure moved up (assuming the 1st row is the top)
            this_dy = 0;
            while all(aligned_movie(this_dy+1,:,i_t)==0)
                this_dy = this_dy + 1; % this counts the displacement in x
            end

        elseif all(aligned_movie(end,:,i_t)==0)
            % the structure moved down (assuming the 1st row is the top)
            this_dy = n_rows;
            while all(aligned_movie(this_dy-1,:,i_t)==0)
                this_dy = this_dy - 1; % this counts the displacement in x
            end

            this_dy = -(n_rows - this_dy); % define down as negative
        else
            % there was no movement
            this_dy = 0;
        end

        % update dx and dy vectors
        dx(i_t) = this_dx;
        dy(i_t) = this_dy;

        if rem(i_t,round(n_t/10))==0
            disp([num2str(round(100 * i_t / n_t)), '% done with alignment'])
        end
    end
    
    % only save the the rows and columns that stayed in frame
    good_cols = [abs(min(dx)) + 1 : n_cols - max(dx)]; % columns that stayed in frame
    good_rows = [max(dy) + 1 : n_rows - abs(min(dy))]; % rows that stayed in frame
    
    % only save the good part of tissue mask
    tissue_mask = tissue_mask(good_rows, good_cols);

    % save the displacement in x and y
    mkdir(fullfile(dataPath, 'psycho6',num2str(channel)));
%     save(fullfile(dataPath, 'psycho6', 'alignment_info.mat'), ...
%         'dx', 'dy', 'tissue_mask');
     save(fullfile(dataPath, 'psycho6',num2str(channel), 'alignment_info.mat'), ...
        'dx', 'dy', 'tissue_mask');
end
    