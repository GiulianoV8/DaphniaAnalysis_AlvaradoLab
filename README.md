# DaphniaAnalysis_AlvaradoLab
Code for analysis of Daphnia collective behavior. Efforts led by Alvarado Lab at Queens College.

Data Visualization: https://docs.google.com/presentation/d/1Bh2P0BvXDtHaI6lCq7zyqyxMzMPgm3qTGtTEco3Orw8/edit?usp=sharing

Inspired by:
Zampetaki, A., Yang, Y., Löwen, H. et al. Dynamical order and many-body correlations in zebrafish show that three is a crowd. Nat Commun 15, 2591 (2024). https://doi.org/10.1038/s41467-024-46426-1

Questions? Contact giulianoyv@gmail.com

# How to use
[TRex](https://doi.org/10.7554/eLife.64000) tracking produces NPZs, converted to CSVs via `NPZtoCSVFIX.py`. These CSVs are used as the data for analysis.

## INPUTS
- base: Filename prefix before "_daphnia#.csv"
- inputDir: Folder containing the trajectory CSV files
- totalN: Total number of tracked Daphnia
- subGroupK: Number of Daphnia used in each Op/Or calculation
- numDraws: Number of random valid frame-subgroup samples
- quarantine_coords: [x, y] coordinates assigned to invalid positions
- frame_rate: Video frame rate in frames per second
- video_length: Video duration in seconds
- selectionMethod: If "random", select Daphnia by random. If "KNN", select by K-nearest neighbors. If "KFN", select by K-Farthest Neighbors.

## OUTPUTS
- all_Op: numDraws-by-1 polarization measurements
- all_Or: numDraws-by-1 rotation measurements
- sampledFrames: Reference-frame index used for each observation
- sampledIDs: Daphnia row indices used for each observation

## Example Call
```
base =  'AItracking20minutes';
inputDir = 'C:\Users\John\Videos\Organism_Motion_Data\AItracking20minutes_csv';
totalN = 30;
subGroupK = 3;
numDraws = 30000;
quarantine_coords: [0,0];
frame_rate: 60;
video_length: 1200;
selectionMethod: 'KNN';

[Op, Or] = plotCollectiveBehaviorOpOr(base,inputDir, totalN, subGroupK, numDraws, quarantine_coords, frame_rate, video_length, selectionMethod);
```
