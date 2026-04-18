function [Obj_Mask_combined,cell_type,size_arr,Mask] = getHealthTypeAndSize(I_stain,Mask)

    [Stain_mask,Size_mask] = MakeHealthMask(I_stain);   

    % ── SOMA objects (original logic preserved) ───────────────────────────────
    [Stain_labeled, ~] = bwlabel(Stain_mask);
    Obj_Mask_soma_raw = Mask .* Stain_labeled;
    unique_soma       = unique(Obj_Mask_soma_raw);
    unique_soma       = unique_soma(2:end);          % drop background
    N_soma            = length(unique_soma);
    [~, Obj_soma] = ismember(Obj_Mask_soma_raw, unique_soma);
    % Obj_soma: 0 = background, 1..N_soma = soma labels

    % ── NEURITE ring objects ──────────────────────────────────────────────────
    dilation_radius   = 20;                          % tune to your resolution
    SE                = strel('disk', dilation_radius);
    [Stain_labeled, ~] = bwlabel(Stain_mask);
    Stain_dilated     = imdilate(logical(Stain_mask), SE);
    Neurite_ring      = Stain_dilated & ~logical(Stain_mask) ...
                                      & logical(Mask);
    % Assign each ring pixel to its nearest soma label (nearest-neighbor)
    [~, nearest_idx]  = bwdist(logical(Stain_mask));
    nearest_label     = Stain_labeled(nearest_idx);
    Neurite_labeled   = double(nearest_label) .* double(Neurite_ring);
    
    unique_neurite    = unique(Neurite_labeled);
    unique_neurite    = unique_neurite(2:end);        % drop background
    N_neurite         = length(unique_neurite);
    
    [~, Obj_neurite]  = ismember(Neurite_labeled, unique_neurite);
    % Obj_neurite: 0 = background, 1..N_neurite = neurite labels
    
    % ── Merge into one sequential label map ───────────────────────────────────
    % Somas  → labels 1        .. N_soma
    % Neurites → labels N_soma+1 .. N_soma+N_neurite
    Obj_Mask_combined = Obj_soma + (Obj_neurite > 0) .* (Obj_neurite + N_soma);
    
    % ── Update outputs ─────────────────────────────────────────────────────────
    Mask_temp                  = Obj_Mask_combined > 0;
    Mask   = Mask_temp;
    
    % cell_type: 1 = soma, 2 = neurite  (one entry per object, same order as labels)
    cell_type = [ones(N_soma, 1); 2 * ones(N_neurite, 1)];


    %%% Size of selected somas only
    valid_soma = Obj_soma > 0;
    stats = regionprops(Obj_soma(valid_soma), Size_mask(valid_soma), 'MeanIntensity');
    soma_size_arr = [stats.MeanIntensity]';
    Size_mask_soma = Size_mask.*valid_soma;

    %%% Size of all stains
    stats = regionprops(Stain_labeled, Size_mask, 'MeanIntensity');
    stain_area  = [stats.MeanIntensity]';
    % For each neurite object label 1..N_neurite, get its parent soma area

    % Build neurite size mask:    
    % every neurite pixel gets the area of its parent soma
    Size_mask_neurite = zeros(size(Size_mask));
    valid_neurite = Obj_neurite > 0;
    Size_mask_neurite(valid_neurite) = stain_area(Neurite_labeled(valid_neurite));

    stats = regionprops(Obj_neurite(valid_neurite), Size_mask_neurite(valid_neurite), 'MeanIntensity');
    neurite_size_arr = [stats.MeanIntensity]';

    % Combine into one final size mask aligned with Obj_Mask_combined
    Size_mask = Size_mask_neurite + Size_mask_soma;

    size_arr = [soma_size_arr;neurite_size_arr];

    function [Stain_mask, Size_mask] = MakeHealthMask(I_stain)
        
        scale = 0.561;
        darkCurrent = 205;
        sigmaBgInput = 0.7;      % light smoothing only for background estimation
        bgRadius1 = 70;          % first-pass morphological background scale
        bgSmooth1 = 4;           % smooth the first background estimate
    
        clipPct2 = 90;           % clip bright objects before second-pass background estimate
        bgSigma2 = 30;           % second-pass haze scale
    
        lowPct1 = 1.0;           % subtract this percentile after step 1, set 0 to disable
        lowPct2 = 1.0;           % subtract this percentile after step 2, set 0 to disable
    
        Iraw = double(I_stain);
    
        I = Iraw - darkCurrent;
        I(I < 0) = 0;
    
        % Step 1: broad background removal by morphological opening
        IbgInput = imgaussfilt(I, sigmaBgInput);
    
        ds = 4;
    
        Ismall = imresize(IbgInput, 1/ds, 'bicubic');
        bg1_small = imopen(Ismall, strel('disk', round(bgRadius1/ds), 0));
        bg1_small = imgaussfilt(bg1_small, bgSmooth1/ds);
        
        bg1 = imresize(bg1_small, size(IbgInput), 'bicubic');
    
        % bg1 = imopen(IbgInput, strel('disk', bgRadius1, 0));
        % bg1 = imgaussfilt(bg1, bgSmooth1);
    
        Istep1 = I - bg1;
        Istep1(Istep1 < 0) = 0;
    
        if lowPct1 > 0
            vals1 = Istep1(Istep1 > 0);
            p1 = prctile(vals1, lowPct1);
            Istep1 = Istep1 - p1;
            Istep1(Istep1 < 0) = 0;
        else
            p1 = 0;
        end
    
        % Step 2: residual haze removal from clipped image
        valsClip = Istep1(Istep1 > 0);
        clipHi = prctile(valsClip, clipPct2);
        tmp = min(Istep1, clipHi);
    
        bg2 = imgaussfilt(tmp, bgSigma2);
    
        Icorr = Istep1 - bg2;
        Icorr(Icorr < 0) = 0;
    
        if lowPct2 > 0
            vals2 = Icorr(Icorr > 0);
            p2 = prctile(vals2, lowPct2);
            Icorr = Icorr - p2;
            Icorr(Icorr < 0) = 0;
        else
            p2 = 0;
        end
    
        %%% Thresholding (Test Donna's code again)
        vals = Icorr(Icorr>10);
        floorVals = vals(vals <= prctile(vals, 70));
        mu_floor = median(floorVals);
        sigma_floor = 1.4826 * mad(floorVals, 1);
        Tauto = mu_floor + 6 * sigma_floor;
        I_thresh = bwareaopen(Icorr > Tauto, 100, 8);
    
        %%% Separating closely spaced neurons here
       
        % Stain_mask = ...
    
        % --- Parameters (tune these) ---
        sigmaSeeds   = 3;     % Gaussian smoothing before peak detection; larger = fewer peaks
        hMaxVal      = [];    % h-maxima suppression threshold; [] = auto (15% of 95th pct)
        minSeedDist  = 8;     % minimum allowed distance (px) between two seed points
        % --------------------------------
        
        % 1. Smooth the corrected image for robust, noise-free peak detection
        Ismooth = imgaussfilt(Icorr, sigmaSeeds);
        
        seeds = false(size(I_thresh));
        CC    = bwconncomp(I_thresh);
        
        hFrac    = 0.10;   % h = this fraction of each blob's local intensity range; LOWER = more splits
        hFloor   = 2.0;    % absolute minimum h to avoid seeding pure noise blobs
        
        % for k = 1:CC.NumObjects
        %     blobMask = false(size(I_thresh));
        %     blobMask(CC.PixelIdxList{k}) = true;
        % 
        %     % Local intensity range for THIS blob only
        %     blobVals  = Ismooth(blobMask);
        %     localMax  = max(blobVals);
        %     localMin  = min(blobVals);
        % 
        %     hAdaptive = max(hFrac * (localMax - localMin), hFloor);
        % 
        %     % Extended maxima scoped to this blob
        %     blobSeeds = imextendedmax(Ismooth .* double(blobMask), hAdaptive) & blobMask;
        %     seeds     = seeds | blobSeeds;
        % end
    
        
        
        % 4. Suppress seeds that are too close together (keeps the dominant peak)
        %    Dilate seeds, then re-find maxima of smoothed image within dilated regions
        if minSeedDist > 0
            seedCC  = bwconncomp(seeds);
            seedImg = zeros(size(seeds));
            for k = 1:seedCC.NumObjects
                idx           = seedCC.PixelIdxList{k};
                [~, bestIdx]  = max(Ismooth(idx));
                seedImg(idx(bestIdx)) = 1;
            end
            seeds = logical(seedImg);
        end
        
        % 5. Dilate seeds slightly so watershed has a stable starting region
        seeds = imdilate(seeds, strel('disk', 2));
        
        % 6. Build the watershed surface: inverted distance transform of the mask.
        %    Distance transform naturally creates ridges between touching blobs.
        D = -bwdist(~I_thresh);
        
        % 7. Impose minima at seed locations AND at background (outside mask).
        %    imimposemin forces watershed to only split at YOUR chosen seeds.
        bgMarkers   = ~imdilate(I_thresh, strel('disk', 3));   % just outside mask edge
        allMarkers  = seeds | bgMarkers;
        D_imposed   = imimposemin(D, allMarkers);
        
        % 8. Run watershed and mask back to original thresholded region
        L = watershed(D_imposed);
        Stain_mask = I_thresh & (L > 0);   % L==0 are watershed boundary lines
        % imagesc(L)
    
    
        %%%
        stats = regionprops(Stain_mask, I_stain, 'Area', 'Perimeter', 'Centroid', 'MeanIntensity', 'PixelIdxList');
        areas = [stats.Area]*scale^2; % Area in pixels
        
        Size_mask = zeros(size(Stain_mask));
        for idx_obj = 1:length(stats)
            Size_mask(stats(idx_obj).PixelIdxList) = stats(idx_obj).Area;
        end
    
    
        % figure;
        % t = tiledlayout(1,5, 'Padding', 'compact', 'TileSpacing', 'compact');
        % 
        % ax1 = nexttile;
        % imagesc(ax1, I_stain);
        % axis(ax1, 'image');
        % colorbar(ax1);
        % title(ax1, 'Raw');
        % 
        % ax2 = nexttile;
        % imagesc(ax2, Istep1);
        % axis(ax2, 'image');
        % colorbar(ax2);
        % title(ax2, 'After stage 1');
        % 
        % ax3 = nexttile;
        % imagesc(ax3, Icorr);
        % axis(ax3, 'image');
        % colorbar(ax3);
        % title(ax3, 'After stage 2');
        % 
        % ax4 = nexttile;
        % imagesc(ax4, I_thresh);
        % axis(ax4, 'image');
        % colorbar(ax4);
        % title(ax4, 'Candidate mask overlay');
        % hold(ax4, 'on');
        % 
        % 
        % ax5 = nexttile;
        % imagesc(ax5, Stain_mask);
        % axis(ax5, 'image');
        % colorbar(ax5);
        % title(ax5, 'Candidate mask overlay');
        % hold(ax5, 'on');
        % 
        % 
        % linkaxes([ax1, ax2, ax3, ax4,ax5], 'xy');
        
        % figure;
        % imshow(I_stain, [200 600]);
        % colormap(gca, gray);
        % axis image;
        % hold on;
        % overlay = cat(3, ones(size(Stain_mask)), zeros(size(Stain_mask)), zeros(size(Stain_mask)));
        % h = imshow(overlay);
        % set(h, 'AlphaData', 0.30 * double(Stain_mask));
        % title('Intensity image with red mask overlay');
    
    end


end

