function relaxed= RelaxationNetworkSimulation(init_tissue_prob)
% This function takes an img and performs the relaxation labelling
% alrgorithm to segment the object pixels from the background
% written by Garrett on 4-27-2021

    p1{1} = init_tissue_prob; % initial probabiligy of being label 1 - i.e. neuronal tissue
    p0{1}= 1 - init_tissue_prob; % initial probabilities of being label 0 - i.e. background
    
    delta=0.1; % step size
    niters = 5;
    for count= 1:niters

        s_i_1= support(p1{count}, 1); % finds support for i'th pixel being label 1
        s_i_0= support(p1{count}, 0); % finds support for i'th pixel being label 0

        p1{count+1} = p1{count} + delta * s_i_1;
        p1_g1_idxs= p1{count+1} > 1; % indices with probability greater than 1
        p1_l0_idxs= p1{count+1} < 0; % indices with probability less than 0
        p1{count+1}(p1_g1_idxs)= 1;
        p1{count+1}(p1_l0_idxs)= 0;

        p0{count+1}= p0{count} + delta * s_i_0;
        p0_g1_idxs= p0{count+1} > 1; % indices with probability greater than 1
        p0_l0_idxs= p0{count+1} < 0; % indices with probability less than 0
        p0{count+1}(p0_g1_idxs)= 1;
        p0{count+1}(p0_l0_idxs)= 0;

        error(count)= sum( abs(p1{count+1} - p1{count}), 'all') + ...
        sum( abs(p0{count+1} - p0{count}), 'all');
        % note, both sums should be equivalent

        if ~all( (p0{count+1}+p1{count+1}) - 1 < 10*eps(1), 'all')
            % The sum of the 2 probabilities should be within 10x machine epsilon
            % away from 1
            error('ERROR: Probabilities do not sum to 1')
        end

    end
    
    relaxed= p0{end}>0.5; % this is the background mask
    
    % visualize the boundaries of neural tissue
%     MakeFigure;
%     imagesc(img); hold on;
%     neural_bound= bwboundaries(round(p1{end}));
%     visboundaries(neural_bound, 'LineStyle', '--','LineWidth', 0.1, 'color', 'red')
%     if nargin==3
%         % relaxation labeling was done assuming ROI locations
%         title('Neural Tissue Labeled by Relaxation Network on ROI mask'); 
%     elseif nargin==2
%         % relaxation labeling was done on the raw recording
%         title('Neural Tissue Labeled by Relaxation Network on Raw Recording'); 
%     end
%     colormap gray;
%     set(gca,'fontsize',15)
%     axis off
end

function s_i= support(img, label)

    [row_max, col_max]= size(img);

    row_matrix= repmat( [1:row_max]', 1      , col_max );
    col_matrix= repmat( [1:col_max] , row_max, 1);

    % 8 touching neighbors
    d_row= [0,1,1,1,0,-1,-1,-1];
    d_col= [1,1,0,-1,-1,-1,0,1];
    
    % initialize support matrix
    s_i= zeros(size(img));

    for in= 1:8
        % iterate through the 8 possible neighbors

        % calculate which pixels have the in'th neighbor - there's a problem
        % with pixels on the border of the image
        allow= (d_row(in)+row_matrix >= 1) & (d_row(in)+row_matrix <= row_max) & ...
                (d_col(in)+col_matrix >= 1) & (d_col(in)+col_matrix <= col_max);
        % binary vectors of allowed row and column indices
        row_allow= row_matrix(allow(:,1),1);
        col_allow= col_matrix(1,allow(1,:));
        
        s_i(row_allow,col_allow)= s_i(row_allow,col_allow) + ...
            support_neighbor(img(row_allow,col_allow),img(row_allow+d_row(in), col_allow+d_col(in)),label);
    end
end

function s_ij= support_neighbor(img_i, img_j, label)
    if ~isequal(size(img_i), size(img_j))
        error('ERROR: Dimensions do not agree. Check your image shifts.')
    end

    p_j1= img_j; % probability neighbor is label 1
    p_j0= 1 - img_j; % probability neighbor is label 0
    
    % r_ij=1 means compatible
    % r_ij=-1 means incompatible
    r_ij_1= 2*(1-abs(label - 1)) - 1;
    r_ij_0= 2*(1-abs(label - 0)) - 1;

    % compute the support output for this step
    s_ij= (r_ij_1 .* p_j1) + (r_ij_0 .* p_j0);
end