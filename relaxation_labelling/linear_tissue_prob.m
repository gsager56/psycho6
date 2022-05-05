function tissue_prob = linear_tissue_prob_HighNeuronThresh(img)
% this function takes in an image and computes an initial probability that
% a pixel belongs to neural tissue. THIS FUNCTION ASSUMES bright pixels are more
% likely to be neural tissue than dark pixels. Therefore, this function
% tries to find a good initial intensity threshold for the neural tissue
% and background. It then assumes the probability pixels between these
% thresholds is neural tissue grows linearly with intensity. The user may
% want to add a nonlinearity

    assert( size(img,3) == 1, 'ERROR: this function assumes a 2D image' )
    img= img ./ max(img, [], 'all'); % normalize img
    
    % all other methods require a neuron and background threshold
    % start by finding the two peaks of intensities
    values= img(:);
    [counts, edges]= histcounts(log10(values));
    [~, idx, ~, prom] = findpeaks(counts);
    prom_thresh = 300; % prominance threshold

    if sum( prom > prom_thresh ) == 2
        % the two peaks correspond to the background and neural tissue

        % find the peak at smaller intensity
        lrg_prom_idx = find(prom > prom_thresh); % indices of the large prominance values

        % peaks for the background and neural tissue
        bkg.I_thresh = 10^edges(idx(lrg_prom_idx(1)));
    elseif sum( prom > prom_thresh ) == 1
        % it could only find a peak for the background or neural tissue
        lrg_prom_idx = find(prom > prom_thresh); % indices of the large prominance values

        % this is the large prominance threshold
        big_thresh = 10^edges(idx(lrg_prom_idx));
        if big_thresh > 0.4
            % this is almost certainly the neural tissue
            % let's find the distribution peak of pixels below the mean
            % image minus standard deviation

            init_bkg_thresh = mean(img, 'all') - ( std(img,0,'all')/2 );
            values= img( img < init_bkg_thresh );
            [counts, edges]= histcounts(log10(values));
            [~, idx, ~, prom] = findpeaks(counts);

            % largest prominance peak corresponds probably to the background
            [~, bkg_prom_idx] = max(prom);
            bkg.I_thresh = 10^edges(idx(bkg_prom_idx));
        else
            % this peak almost certainly corresponds to the background
            bkg.I_thresh = big_thresh;
            init_neuron_thresh = mean(img, 'all') + 2*std(img,0,'all');

            values= img( img > init_neuron_thresh );
            [counts, edges]= histcounts(log10(values));
            [~, idx, ~, prom] = findpeaks(counts);

            % largest prominance peak corresponds probably to the neural tissue
            % since neuron is so weak, let's use this peak as the threshold
            [~, neuron_prom_idx] = max(prom);
            neuron.I_thresh = 10^edges(idx(neuron_prom_idx));
        end
    else
        % something unexpected happened
        error('ERROR: expected 2 peaks in the distribution of correlations')
    end

    % find probabilities each pixel is part of a neuron
    slope= 1 / (neuron.I_thresh-bkg.I_thresh); % slope of linear equation
    b= -slope*bkg.I_thresh; % y-intercept of linear equation
    tissue_prob= slope.*img + b;
    tissue_prob(tissue_prob<0)= 0; % probabilities less than 0 are called 0
    tissue_prob(tissue_prob>1)= 1; % probabilities greater than 1 are called 1
