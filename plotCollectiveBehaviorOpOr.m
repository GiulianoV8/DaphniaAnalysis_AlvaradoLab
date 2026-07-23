function [all_Op, all_Or, sampledFrames, sampledIDs] = plotCollectiveBehaviorOpOrDensity(base, inputDir, totalN, subGroupK, numDraws, quarantine_coords, frame_rate, video_length, selectionMethod)
%PLOTCOLLECTIVEBEHAVIOROPOR Plot the joint density of Op and Or.
% Op and Or are calculated only from Daphnia that are not quarantined
% and have a finite, nonzero, threshold-valid velocity in that frame.
%
% NOTE: pick quarantine_coords well outside your arena's real coordinate
% range. If a genuine detection ever lands exactly on quarantine_coords,
% it will be wrongly treated as a missing detection and excluded.
%
% INPUTS
%   base               Filename prefix before "_daphnia#.csv"
%   inputDir           Folder containing the trajectory CSV files
%   totalN             Total number of tracked Daphnia
%   subGroupK          Number of Daphnia used in each Op/Or calculation
%   numDraws           Number of random valid frame-subgroup samples
%   quarantine_coords  [x, y] coordinates assigned to invalid positions
%   frame_rate         Video frame rate in frames per second
%   video_length       Video duration in seconds
%   selectionMethod    If "random", select Daphnia by random. If "KNN", select
%                      by K-nearest neighbors. If "KFN", select by
%                      K-Farthest Neighbors.
%
% OUTPUTS
%   all_Op             numDraws-by-1 polarization measurements
%   all_Or             numDraws-by-1 rotation measurements
%   sampledFrames      Reference-frame index used for each observation
%   sampledIDs         Daphnia row indices used for each observation
%
% EXAMPLE
%   [Op, Or] = plotCollectiveBehaviorOpOr( ...
%       'AItracking20minutes', ...
%       'C:\Users\vrspr\Videos\Organism_Motion_Data\AItracking20minutes_csv', ...
%       30, 3, 500, [-1000, -1000], 60, 1200);

    % Check parameter validity
    validateOpOrInputs(totalN, subGroupK, numDraws, quarantine_coords, selectionMethod);

    function validateOpOrInputs(totalN, subGroupK, numDraws, quarantine_coords, selectionMethod)
    %VALIDATEOPORINPUTS Basic sanity checks on the function's arguments.
        if totalN < 1 || totalN ~= floor(totalN)
            error('totalN must be a positive integer.');
        end
        if subGroupK < 2 || subGroupK ~= floor(subGroupK)
            error('subGroupK must be an integer of at least 2.');
        end
        if subGroupK > totalN
            error('subGroupK (%d) cannot exceed totalN (%d).', subGroupK, totalN);
        end
        if numDraws < 1 || numDraws ~= floor(numDraws)
            error('numDraws must be a positive integer.');
        end
        if numel(quarantine_coords) ~= 2
            error('quarantine_coords must be a two-element vector: [x, y].');
        end
        if selectionMethod ~= "random" && selectionMethod ~= "KNN" && selectionMethod ~= "KFN"
            error("selection method must be 'random', 'KNN', or 'KFN'");
        end
    end

    % Obtain corrected frame-aligned matrices from quarantine.m
    [ref_time, ~, ~, V_total, Xmat, Ymat, Vx, Vy] = quarantine(base, inputDir, totalN, quarantine_coords, frame_rate, video_length);

    ref_time = ref_time(:).';
    timeLength = numel(ref_time);
    
    % Confirm matrix formats
    checkMatrixSizes(Xmat, Ymat, Vx, Vy, V_total, totalN, timeLength);

    function checkMatrixSizes(Xmat, Ymat, Vx, Vy, V_total, totalN, timeLength)
    %CHECKMATRIXSIZES Confirm quarantine.m returned N-by-T matrices as expected.
        expectedSize = [totalN, timeLength];
        matrixNames = {'Xmat', 'Ymat', 'Vx', 'Vy', 'V_total'};
        matrices = {Xmat, Ymat, Vx, Vy, V_total};
    
        for matrixNumber = 1:numel(matrices)
            if ~isequal(size(matrices{matrixNumber}), expectedSize)
                error('%s must be an N-by-T matrix returned by quarantine.m.', ...
                      matrixNames{matrixNumber});
            end
        end
    end

    % Find frames with valid subgroups
    [usableDaphnia, validFrames] = findUsableFrames(Xmat, Ymat, Vx, Vy, V_total, quarantine_coords, subGroupK);

    function [usableDaphnia, validFrames] = findUsableFrames(Xmat, Ymat, Vx, Vy, V_total, quarantine_coords, subGroupK)
    %FINDUSABLEFRAMES Identify, for every Daphnia and frame, whether it is valid
    % and frames where at least subGroupK such Daphnia are available.
    
        quarantineTolerance = 1e-4;
        
        % Check if Daphnia quarantined
        inQuarantine = ...
            abs(Xmat - quarantine_coords(1)) <= quarantineTolerance & ...
            abs(Ymat - quarantine_coords(2)) <= quarantineTolerance;
        
        % Identify valid Daphnia
        usableDaphnia = ...
            ~inQuarantine & ...
            isfinite(Xmat) & isfinite(Ymat) & ...
            isfinite(Vx) & isfinite(Vy) & isfinite(V_total) & ...
            V_total > eps;
    
        % Get frames where Op and Or can be calculated with subgroup of K
        numberUsablePerFrame = sum(usableDaphnia, 1);
        validFrames = find(numberUsablePerFrame >= subGroupK);
    
        fprintf('Total usable Daphnia-frame velocities: %d\n', nnz(usableDaphnia));
        fprintf('Maximum usable Daphnia in one frame: %d\n', max(numberUsablePerFrame));
        fprintf('Frames containing at least K=%d usable Daphnia: %d\n', ...
                subGroupK, numel(validFrames));
    
        if isempty(validFrames)
            error([ ...
                'No frame contains at least %d valid moving Daphnia. ' ...
                'Check the quarantine coordinates, frame rate, video length, ' ...
                'and the velocity threshold in quarantine.m.'], subGroupK);
        end
    end

    % Calculate Op and Or
    [all_Op, all_Or, sampledFrames, sampledIDs] = sampleOpOr(Xmat, Ymat, Vx, Vy, V_total, usableDaphnia, validFrames, subGroupK, numDraws, selectionMethod);

    plotDensityMap(all_Op, all_Or, subGroupK, numel(all_Op));
end

function [all_Op, all_Or, sampledFrames, sampledIDs] = sampleOpOr(Xmat, Ymat, Vx, Vy, V_total, usableDaphnia, validFrames, subGroupK, numDraws, method)
%SAMPLEOPOR For all valid frames compute Op and Or for a subgroup of K
%valid Daphnia, where the subgroup is the k-nearest neighbors of a
%randomly chosen valid "center" Daphnia in that frame.

    all_Op = NaN(numDraws, 1);
    all_Or = NaN(numDraws, 1);
    sampledFrames = NaN(numDraws, 1);
    sampledIDs = NaN(numDraws, subGroupK);
    
    observationCount = 0;
    attemptCount = 0;
    maximumAttempts = max(100, 50 * numDraws);
    
    h = waitbar(0, sprintf('Sampling %d valid subgroup-frames...', numDraws));
    cleanupWaitbar = onCleanup(@() closeWaitbarIfOpen(h));
    
    while observationCount < numDraws && attemptCount < maximumAttempts
        attemptCount = attemptCount + 1;
    
        % Pick a random valid frame
        t = validFrames(randi(numel(validFrames)));
    
        % IDs that are usable in this frame
        availableIDs = find(usableDaphnia(:, t));
    
        % Safety check (should already be true by construction)
        if numel(availableIDs) < subGroupK
            continue;
        end
    
        % --- NEW: random center + k-nearest neighbors in position space ---
    
        % Positions of all valid Daphnia in this frame
        X_valid = Xmat(availableIDs, t);
        Y_valid = Ymat(availableIDs, t);
    
        % Also allow purely random subgroup selection (without k-NN) for variety:
        if method == "random"
            % purely random k distinct IDs from availableIDs
            rndOrder = randperm(numel(availableIDs), subGroupK);
            chosenIDs = availableIDs(rndOrder);
        else
            % Choose a random center among valid IDs
            centerIdxInAvailable = randi(numel(availableIDs));
            centerX = X_valid(centerIdxInAvailable);
            centerY = Y_valid(centerIdxInAvailable);
        
            % Squared Euclidean distance to all valid Daphnia in frame
            dx = X_valid - centerX;
            dy = Y_valid - centerY;
            dist2 = dx.^2 + dy.^2;
            
            if method == "KNN"
                % Sort by distance and take k nearest (including the center itself)
                [~, sortOrder] = sort(dist2, 'ascend');
                nearestInAvailable = sortOrder(1:subGroupK);
                % Global indices of chosen subgroup
                chosenIDs = availableIDs(nearestInAvailable);
            elseif method == "KFN" 
                % Sort by distance and take k nearest (including the center itself)
                [~, sortOrder] = sort(dist2, 'descend');
                farthestInAvailable = sortOrder(1:subGroupK);
                % Global indices of chosen subgroup
                chosenIDs = availableIDs(farthestInAvailable);
            end
        end
        % ----------------------------------------------------------
    
        % Record positions and velocity data of chosen Daphnia
        curr_X = Xmat(chosenIDs, t);
        curr_Y = Ymat(chosenIDs, t);
        curr_Vx = Vx(chosenIDs, t);
        curr_Vy = Vy(chosenIDs, t);
        curr_speed = V_total(chosenIDs, t);
    
        % Calculate Op and Or
        Op = polarizationOrder(curr_Vx, curr_Vy, curr_speed);
        Or = rotationOrder(curr_X, curr_Y, curr_Vx, curr_Vy, curr_speed);
    
        if ~isfinite(Op) || ~isfinite(Or)
            % e.g. a chosen Daphnia sat exactly at the subgroup's center of mass
            continue;
        end
    
        % Clamp to [0, 1] against tiny floating-point excursions
        Op = min(max(Op, 0), 1);
        Or = min(max(Or, 0), 1);
    
        % Record Op / Or data
        observationCount = observationCount + 1;
        all_Op(observationCount) = Op;
        all_Or(observationCount) = Or;
        sampledFrames(observationCount) = t;
        sampledIDs(observationCount, :) = chosenIDs(:).';
    
        if isgraphics(h)
            waitbar(observationCount / numDraws, h, ...
                sprintf('Valid observations: %d/%d', observationCount, numDraws));
        end
    end
    
    % Trim arrays if fewer than numDraws observations were obtained
    all_Op = all_Op(1:observationCount);
    all_Or = all_Or(1:observationCount);
    sampledFrames = sampledFrames(1:observationCount);
    sampledIDs = sampledIDs(1:observationCount, :);
    
    if observationCount == 0
        error('No geometrically valid Op/Or observations could be calculated.');
    elseif observationCount < numDraws
        warning(['Only %d of %d requested observations were calculated ' ...
            'after %d sampling attempts.'], observationCount, numDraws, attemptCount);
    end
end

function Op = polarizationOrder(Vx, Vy, speed)
%POLARIZATIONORDER Magnitude of the subgroup's average heading unit vector. Op = 1 means every member is moving in the same direction.
    normalized_vx = Vx ./ speed;
    normalized_vy = Vy ./ speed;

    Op = hypot(mean(normalized_vx), mean(normalized_vy));
end

function Or = rotationOrder(X, Y, Vx, Vy, speed)
%ROTATIONORDER Magnitude of the subgroup's average (radial-direction x
%velocity-direction) term, about the subgroup's own center of mass. Or
%is close to 1 when the group is milling/rotating together.
%
% Returns NaN if any member sits exactly at the center of mass, where the
% radial direction is undefined.

    % Calculate center of mass
    comX = mean(X);
    comY = mean(Y);

    % Calculate radius of milling
    r_x = X - comX;
    r_y = Y - comY;
    r_mag = hypot(r_x, r_y);

    % Ignore invalid/insignificant millings
    if any(~isfinite(r_mag)) || any(r_mag <= eps)
        Or = NaN;
        return;
    end

    normalized_rx = r_x ./ r_mag;
    normalized_ry = r_y ./ r_mag;
    normalized_vx = Vx ./ speed;
    normalized_vy = Vy ./ speed;

    crossProduct = normalized_rx .* normalized_vy - normalized_ry .* normalized_vx;
    Or = abs(mean(crossProduct));
end

%% ------------------------------------------------------------------
function plotDensityMap(all_Op, all_Or, subGroupK, observationCount)
%PLOTDENSITYMAP Dot-density style Op-Or plot.

    figure('Color', 'w', 'Position', [100, 100, 520, 420], ...
        'Name', 'Collective Behavior Op-Or Density');
    
    thr = 0.65;
    hold on;
    
    patch([0 0 thr thr], [thr 1 1 thr], [0.70 0.88 0.70], ...
        'FaceAlpha', 0.25, 'EdgeColor', 'none');
    patch([thr 1 1 thr], [0 0 thr thr], [0.86 0.76 0.90], ...
        'FaceAlpha', 0.25, 'EdgeColor', 'none');
    patch([0 thr thr 0], [0 0 thr thr], [0.93 0.83 0.70], ...
        'FaceAlpha', 0.25, 'EdgeColor', 'none');
    patch([thr 1 1 thr], [thr thr 1 1], [0.96 0.88 0.55], ...
        'FaceAlpha', 0.35, 'EdgeColor', 'none');
    
    scatter(all_Op, all_Or, 7, 'k', 'filled', ...
        'MarkerFaceAlpha', 0.15, 'MarkerEdgeAlpha', 0.15);
    
    plot([thr thr], [0 1], 'k-', 'LineWidth', 1);
    plot([0 1], [thr thr], 'k-', 'LineWidth', 1);
    
    sigOp = all_Op >= thr;
    sigOr = all_Or >= thr;
    sigBoth = sigOp & sigOr;
    
    nOpOnly = nnz(sigOp & ~sigOr);
    nOrOnly = nnz(sigOr & ~sigOp);
    nBoth = nnz(sigBoth);
    nNone = nnz(~sigOp & ~sigOr);
    
    pOpOnly = 100 * nOpOnly / observationCount;
    pOrOnly = 100 * nOrOnly / observationCount;
    pBoth = 100 * nBoth / observationCount;
    pNone = 100 * nNone / observationCount;
    
    text(0.10, 0.86, sprintf('%d%%', round(pOrOnly)), ...
        'Color', 'k', 'FontSize', 11, 'FontWeight', 'bold');
    text(0.80, 0.10, sprintf('%d%%', round(pOpOnly)), ...
        'Color', 'k', 'FontSize', 11, 'FontWeight', 'bold');
    text(0.10, 0.10, sprintf('%d%%', round(pNone)), ...
        'Color', 'k', 'FontSize', 11, 'FontWeight', 'bold');
    text(0.80, 0.80, sprintf('%d%%', round(pBoth)), ...
        'Color', 'k', 'FontSize', 11, 'FontWeight', 'bold');
    
    xlabel('O_p Schooling \rightarrow');
    ylabel('O_r Milling \rightarrow');
    title(sprintf('Op-Or Density, K = %d, n = %d', subGroupK, observationCount));
    
    axis square;
    xlim([0 1]);
    ylim([0 1]);
    set(gca, 'Box', 'on', 'LineWidth', 1, 'FontSize', 12);
    
    hold off;
end

%% ------------------------------------------------------------------
function closeWaitbarIfOpen(h)
%CLOSEWAITBARIFOPEN Close a waitbar without throwing a cleanup error.
    if ~isempty(h) && isgraphics(h)
        close(h);
    end
end


[Op, Or] = plotCollectiveBehaviorOpOrDensity('AItracking20minutes', 'C:\Users\vrspr\Videos\Organism_Motion_Data\AItracking20minutes_csv', 30, 2, 3000, [0, 0], 60, 1200, "KFN");
[Op, Or] = plotCollectiveBehaviorOpOrDensity('AItracking20minutes', 'C:\Users\vrspr\Videos\Organism_Motion_Data\AItracking20minutes_csv', 30, 3, 3000, [0, 0], 60, 1200, "KFN");
[Op, Or] = plotCollectiveBehaviorOpOrDensity('AItracking20minutes', 'C:\Users\vrspr\Videos\Organism_Motion_Data\AItracking20minutes_csv', 30, 4, 3000, [0, 0], 60, 1200, 'KFN');
[Op, Or] = plotCollectiveBehaviorOpOrDensity('AItracking20minutes', 'C:\Users\vrspr\Videos\Organism_Motion_Data\AItracking20minutes_csv', 30, 5, 3000, [0, 0], 60, 1200, "KFN");
[Op, Or] = plotCollectiveBehaviorOpOrDensity('AItracking20minutes', 'C:\Users\vrspr\Videos\Organism_Motion_Data\AItracking20minutes_csv', 30, 7, 3000, [0, 0], 60, 1200, 'KFN');
[Op, Or] = plotCollectiveBehaviorOpOrDensity('AItracking20minutes', 'C:\Users\vrspr\Videos\Organism_Motion_Data\AItracking20minutes_csv', 30, 10, 3000, [0, 0], 60, 1200, "KFN");
[Op, Or] = plotCollectiveBehaviorOpOrDensity('AItracking20minutes', 'C:\Users\vrspr\Videos\Organism_Motion_Data\AItracking20minutes_csv', 30, 15, 3000, [0, 0], 60, 1200, "KFN");
[Op, Or] = plotCollectiveBehaviorOpOrDensity('AItracking20minutes', 'C:\Users\vrspr\Videos\Organism_Motion_Data\AItracking20minutes_csv', 30, 30, 3000, [0, 0], 60, 1200, "KFN");


