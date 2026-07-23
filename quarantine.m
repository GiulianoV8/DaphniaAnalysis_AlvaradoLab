function [ref_time, averageCumulativeAngle, cumulativeAnglesMatrix, V_total, Xmat, Ymat, Vx, Vy] = quarantine(base, inputDir, N, quarantine_coords, frame_rate, video_length)
%QUARANTINE Align Daphnia tracks to a common frame grid.
%
% Missing or invalid positions are placed at quarantine_coords. Velocities are calculated only when both endpoint positions are genuine detections.
% Vx, Vy, and V_total are returned as N-by-T matrices so that other analyses can use the corrected data directly.
%
% OUTPUT MATRIX CONVENTION
%   Xmat(i,t), Ymat(i,t)  Position of Daphnia i at frame t.
%   Vx(i,t), Vy(i,t)      Velocity from frame t-1 to frame t.
%   V_total(i,t)          Speed from frame t-1 to frame t.
%
% Column 1 of each velocity matrix is zero because no previous frame exists.
% Invalid, quarantined, and below-threshold velocities are also zero.

    % Check parameter validity
    validateQuarantineInputs(quarantine_coords, N, frame_rate, video_length);

    %% ------------------------------------------------------------------
    function validateQuarantineInputs(quarantine_coords, N, frame_rate, video_length)
    %VALIDATEQUARANTINEINPUTS Basic sanity checks on the function's arguments.
        if numel(quarantine_coords) ~= 2
            error('quarantine_coords must be a two-element vector: [x, y].');
        end
        if N < 1 || N ~= floor(N)
            error('N must be a positive integer.');
        end
        if frame_rate <= 0
            error('frame_rate must be positive.');
        end
        if video_length < 0
            error('video_length cannot be negative.');
        end
    end

    % Fraction of v* used as the movement threshold.
    f = 0.2;

    % Load data
    [X, Y, Time] = loadTrajectoryFiles(base, inputDir, N);
    % Build absolute timeline
    [ref_time, timeLength] = buildReferenceTime(frame_rate, video_length);
    
    % Cumulative angle matrix filled with default NaN
    cumulativeAnglesMatrix = NaN(N, timeLength);

    % Position matrix filled with default quarantine coords
    Xmat = repmat(quarantine_coords(1), N, timeLength);
    Ymat = repmat(quarantine_coords(2), N, timeLength);

    % Velocity matrix filled with default zeros
    Vx = zeros(N, timeLength);
    Vy = zeros(N, timeLength);
    V_total = zeros(N, timeLength);

    h = waitbar(0, 'Aligning positions and calculating velocities...');
    cleanupWaitbar = onCleanup(@() closeWaitbarIfOpen(h));

    % Populate data matrices for each daphnia
    for i = 1:N
        [Xmat(i, :), Ymat(i, :), Vx(i, :), Vy(i, :), V_total(i, :), cumulativeAnglesMatrix(i, :)] = alignOneDaphnia(X{i}, Y{i}, Time{i}, ref_time, quarantine_coords, f, i);

        waitbar(i / N, h, sprintf('Processed Daphnia %d/%d', i, N));
    end
    
    % Calculate cumulative angles
    averageCumulativeAngle = mean(cumulativeAnglesMatrix, 1, 'omitnan');
end

%% ------------------------------------------------------------------
function [X, Y, Time] = loadTrajectoryFiles(base, inputDir, N)
%LOADTRAJECTORYFILES Read each Daphnia's raw X/Y/Time columns from disk.
    X = cell(1, N);
    Y = cell(1, N);
    Time = cell(1, N);

    h = waitbar(0, 'Loading Daphnia trajectory files...');
    cleanupWaitbar = onCleanup(@() closeWaitbarIfOpen(h));

    for i = 0:N-1
        csvFile = fullfile(inputDir, sprintf('%s_daphnia%d.csv', base, i));

        if ~isfile(csvFile)
            error('Could not find trajectory file: %s', csvFile);
        end

        data = readtable(csvFile);
        requiredColumns = {'X', 'Y', 'Time'};

        if ~all(ismember(requiredColumns, data.Properties.VariableNames))
            error('File %s must contain X, Y, and Time columns.', csvFile);
        end

        X{i+1} = data.X;
        Y{i+1} = data.Y;
        Time{i+1} = data.Time;

        waitbar((i + 1) / N, h, sprintf('Loading data... (%d/%d)', i + 1, N));
    end
end

%% ------------------------------------------------------------------
function [ref_time, timeLength] = buildReferenceTime(frame_rate, video_length)
%BUILDREFERENCETIME Construct an absolute timeline for frame reference
    time_interval = 1 / frame_rate;
    ref_time = 0:time_interval:video_length;
    timeLength = numel(ref_time);

    if timeLength < 3
        error('frame_rate and video_length must produce at least three frames.');
    end
end

%% ------------------------------------------------------------------
function [xRow, yRow, vxRow, vyRow, speedRow, angleRow] = alignOneDaphnia(X_raw, Y_raw, Time_raw, ref_time, quarantine_coords, f, daphniaIndex)
%ALIGNONEDAPHNIA Align one Daphnia's raw track onto the common frame grid,
%compute its frame-to-frame velocity, and its cumulative turning angle.
%
% Returns 1-by-timeLength rows. If the track can't be used (mismatched
% column lengths, no finite timestamps, etc.), the defaults below are
% returned unchanged: quarantined position, zero velocity, NaN angle.

    timeLength = numel(ref_time);

    % Position data filled with default quarantine coords
    xRow = repmat(quarantine_coords(1), 1, timeLength);
    yRow = repmat(quarantine_coords(2), 1, timeLength);

    % Velocity/Speed data filled with default quarantine coords
    vxRow = zeros(1, timeLength);
    vyRow = zeros(1, timeLength);
    speedRow = zeros(1, timeLength);
    angleRow = NaN(1, timeLength);

    X_raw = X_raw(:);
    Y_raw = Y_raw(:);
    Time_raw = Time_raw(:);

    % Check for invalid data (container size should match)
    if ~(numel(X_raw) == numel(Y_raw) && numel(X_raw) == numel(Time_raw))
        warning('Daphnia #%d has mismatched X, Y, and Time lengths. Skipping.', daphniaIndex);
        return;
    end

    % --- Fill invalid positions with the quarantine coordinates ---
    invalidPosition = ~isfinite(X_raw) | ~isfinite(Y_raw);

    X_filled = X_raw;
    Y_filled = Y_raw;
    X_filled(invalidPosition) = quarantine_coords(1);
    Y_filled(invalidPosition) = quarantine_coords(2);

    % Check for valid timesteps
    validTime = isfinite(Time_raw);
    sourceIndices = find(validTime);

    if isempty(sourceIndices)
        warning('Daphnia #%d has no finite timestamps. Skipping.', daphniaIndex);
        return;
    end

    % --- Map each recorded timestamp onto its nearest reference frame ---
    matchedIndices = interp1(ref_time, 1:timeLength, Time_raw(sourceIndices), 'nearest', NaN);

    % Mark valid indices
    inRange = isfinite(matchedIndices) & matchedIndices >= 1 & matchedIndices <= timeLength;

    % Only take indices matched to standardized timeline
    sourceIndices = sourceIndices(inRange);

    % Take matched indices and round to whole number
    matchedIndices = round(matchedIndices(inRange));

    if isempty(matchedIndices)
        warning('Daphnia #%d has no timestamps inside the video range. Skipping.', daphniaIndex);
        return;
    end

    % If multiple rows map to one frame, retain the final row.
    [matchedIndices, uniqueLocations] = unique(matchedIndices, 'last');
    sourceIndices = sourceIndices(uniqueLocations);

    positionIsValid = false(1, timeLength);
    
    % Fill data with matched indices
    xRow(matchedIndices) = X_filled(sourceIndices).';
    yRow(matchedIndices) = Y_filled(sourceIndices).';

    % Truth map that tells if a index is a quarantine/invalid index
    positionIsValid(matchedIndices) = ~invalidPosition(sourceIndices);

    % --- Frame-to-frame velocity (only between two genuine detections) ---
    dt = diff(ref_time);
    Vx_step = diff(xRow) ./ dt;
    Vy_step = diff(yRow) ./ dt;

    % Check if velocities are valid calculated between two valid points
    validVelocity = positionIsValid(1:end-1) & positionIsValid(2:end) & isfinite(dt) & dt > 0;

    % Set velocities at quarantine to zero
    Vx_step(~validVelocity) = 0;
    Vy_step(~validVelocity) = 0;

    % Calculate speed
    speed_step = hypot(Vx_step, Vy_step);
    thresholdVelocity = estimateSpeedThreshold(speed_step, validVelocity, f, daphniaIndex);

    % Store interval t-1 -> t in column t. Velocities that fail the validity or v* threshold remain zero.
    vxRow(2:end) = Vx_step .* thresholdVelocity;
    vyRow(2:end) = Vy_step .* thresholdVelocity;
    speedRow(2:end) = speed_step .* thresholdVelocity;

    angleRow = cumulativeTurningAngle(Vx_step, Vy_step, thresholdVelocity, timeLength);
end

%% ------------------------------------------------------------------
function thresholdVelocity = estimateSpeedThreshold(speed_step, validVelocity, f, daphniaIndex)
%ESTIMATESPEEDTHRESHOLD Flag velocity intervals whose speed clears a
%fitted movement threshold v*.
%
% Genuine nonzero speeds are assumed to roughly follow an exponential
% distribution, speed ~ exp(-speed / v*). v* is estimated by histogram-
% binning those speeds and fitting a line to log(probability) vs. speed;
% the fitted decay constant is v*. Intervals with speed >= f * v* are
% kept. If there isn't enough variation in the data to fit v* reliably,
% every genuine nonzero velocity is kept instead.

    % Check for valid speeds (valid index, not NaN or inf, greater than zero)
    speedsForFit = speed_step(validVelocity & isfinite(speed_step) & speed_step > 0);

    % Mark all genuine valid speeds as valid (just filling default list)
    thresholdVelocity = validVelocity & isfinite(speed_step) & speed_step > 0;

    if numel(speedsForFit) < 2 || range(speedsForFit) == 0
        warning(['Daphnia #%d has insufficient speed variation for a ' ...
                 'stable v* fit; all genuine nonzero velocities are retained.'], daphniaIndex);
        return;
    end

    % Build speed histogram
    binSize = 0.5;
    numberOfBins = max(2, ceil(range(speedsForFit) / binSize));
    [counts, edges] = histcounts(speedsForFit, numberOfBins);

    if sum(counts) == 0
        return;
    end

    probabilities = counts / sum(counts);
    binCenters = (edges(1:end-1) + edges(2:end)) / 2;
    % Remove empty bins
    fitMask = probabilities > 0 & isfinite(probabilities);

    if nnz(fitMask) < 2
        return;
    end

    % Fit speed distribution
    fitCoefficients = polyfit(binCenters(fitMask), log(probabilities(fitMask)), 1);
    slope = fitCoefficients(1);

    if ~isfinite(slope) || abs(slope) <= eps
        return;
    end

    % Set threshold
    v_star = 1 / abs(slope);

    % Threshold velocities based on v_star * f
    thresholdVelocity = validVelocity & isfinite(speed_step) & speed_step >= v_star * f;
end

%% ------------------------------------------------------------------
function cumulativeAngle = cumulativeTurningAngle(Vx_step, Vy_step, thresholdVelocity, timeLength)
%CUMULATIVETURNINGANGLE Running turning angle (in revolutions), skipping
%intervals that aren't threshold-valid instead of joining across them.
%
% A turning increment requires two consecutive threshold-valid velocity
% intervals. Increments during a gap are treated as zero, so the running
% sum holds steady across the gap rather than resetting or being
% spuriously incremented by the jump. Columns 1-2 of the returned row are
% always NaN (no turning angle is defined yet); the value stays NaN until
% the first valid turn occurs.

    % Create empty row
    cumulativeAngle = NaN(1, timeLength);

    % Find valid consecutive velocities
    adjacentVelocity = thresholdVelocity(1:end-1) & thresholdVelocity(2:end);

    % Hypotenuse
    denominator = Vx_step(1:end-1).^2 + Vy_step(1:end-1).^2;
    % Check for invalid denominator (e.g., divide by zero)
    adjacentVelocity = adjacentVelocity & isfinite(denominator) & denominator > eps;

    if ~any(adjacentVelocity)
        return;
    end

    % Gaps will be represented by zero 
    delta_phi = zeros(1, timeLength - 2);

    % Calculate change in velocity
    dVx = Vx_step(2:end) - Vx_step(1:end-1);
    dVy = Vy_step(2:end) - Vy_step(1:end-1);

    % Cross product for turn
    numerator = Vx_step(1:end-1) .* dVy - Vy_step(1:end-1) .* dVx;
    delta_phi(adjacentVelocity) = numerator(adjacentVelocity) ./ denominator(adjacentVelocity);

    % Quantify by revolutions
    runningAngle = cumsum(delta_phi / (2 * pi));

    firstValidTurn = find(adjacentVelocity, 1, 'first');
    runningAngle(1:firstValidTurn-1) = NaN;

    cumulativeAngle(3:end) = runningAngle;
end

%% ------------------------------------------------------------------
function closeWaitbarIfOpen(h)
%CLOSEWAITBARIFOPEN Close a waitbar without throwing a cleanup error.
    if ~isempty(h) && isgraphics(h)
        close(h);
    end
end