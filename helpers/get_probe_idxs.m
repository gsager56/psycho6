function probe_idxs = get_probe_idxs(epoch_trace, param)

    probe_epochs = param.probe_epochs;
    probe_bool = ismember( epoch_trace, probe_epochs );
    CC = bwconncomp( probe_bool );
    
    assert( CC.NumObjects == 2, 'There should only be two times the probe is presented' )
    
    if length(CC.PixelIdxList{1}) >= length(CC.PixelIdxList{2})
        % the first probe has more than or the same number of frames as the
        % second probe
        long_probe = 1; short_probe = 2;
    else
        long_probe = 2; short_probe = 1;
    end
    long_probe_idxs = CC.PixelIdxList{long_probe};
    short_probe_idxs = CC.PixelIdxList{short_probe};
    short_probe_trace = epoch_trace( short_probe_idxs );
    
    % find alignment that maximizes the cross-correlation between the two
    % probes
    Corr = zeros( length(long_probe_idxs) - length(short_probe_idxs) + 1, 1 );
    for delay = 1 : length(Corr)
        % DO NOT use matlab's xcorr function to do this!!!!
        % xcorr weights large values more than small values, which is not
        % what we want here
        Corr(delay) = mean(epoch_trace(long_probe_idxs( delay : length(short_probe_idxs)+delay-1 ) ) == short_probe_trace);
    end
    [~, delay] = max(Corr);
    cropped_long_trace = epoch_trace(long_probe_idxs( delay : length(short_probe_idxs)+delay-1 ) );
    
    % find first index where the two signals are equal
    start_idx = find(cropped_long_trace == short_probe_trace, 1, 'first');
    
    % find last index where the two signals are equal
    end_idx = find(cropped_long_trace == short_probe_trace, 1, 'last');
    
    probe_idxs = cell(2,1);
    probe_idxs{1} = start_idx + delay - 1 : end_idx + delay - 1;
    probe_idxs{2} = start_idx : end_idx;
    
    probe_idxs{1} = probe_idxs{1} + long_probe_idxs(1) - 1;
    probe_idxs{2} = probe_idxs{2} + short_probe_idxs(1) - 1;
    
    probe_match = mean(epoch_trace(probe_idxs{1}) == epoch_trace(probe_idxs{2}));
    assert( probe_match > 0.95, 'Could not find good alignment between the first and second probe' )
    
    
    if false
        % here's code to plot the aligned first and second probes
        figure;
        subplot(2,1,1)
        plot(epoch_trace(probe_idxs{1}))
        title('probe 1')

        subplot(2,1,2)
        plot(epoch_trace(probe_idxs{2}))
        title('probe 2')
    end
    

