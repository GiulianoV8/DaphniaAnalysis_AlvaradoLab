import numpy as np
import matplotlib.pyplot as plt
import os
import scipy.spatial
import pandas as pd
import tkinter as tk
from tkinter.filedialog import askdirectory

#def get_directory():
#    root = tk.Tk()
   # root.withdraw()

   # return askdirectory(title='select directory/folder')

base = "Gravid"

# Commenting out this line to test tkinter
input_dir = r"C:\Users\vrspr\OneDrive\Desktop\Daphnia\Daphnia_experiments\experiments\Andrew Summer project\Gravid\videos_pv\data";

output_dir = r"C:\Users\vrspr\Videos\Organism_Motion_Data\\" + base + "_csv"

# Number of daphnia N
N = 31
# Select the range of frame to be analyzed
START = 0
STOP = 35945
# Ensure the output directory exists
os.makedirs(output_dir, exist_ok=True)
 
for i in range(N):
    input_file = os.path.join(input_dir, f"{base}{i}.npz")
    output_file = os.path.join(output_dir, f"{base}_daphnia{i}.csv")

    with np.load(input_file, allow_pickle=True) as npz:
        X = npz["X#wcentroid"]
        Y = -npz["Y#wcentroid"]
        time = npz["time"]
        frame = npz["frame"]
        speed_centroid = npz["SPEED#wcentroid"]
        missing = npz["missing"].astype(bool)

    # Only trim to the experiment window - no missing-frame filtering.
    # time/frame stay complete and gap-free for every animal, which is
    # what fixes ref_time. inf values pass through as-is; MATLAB handles them.
    frame_mask = (frame >= START) & (frame < STOP)

    X_out = X[frame_mask]
    Y_out = Y[frame_mask]
    speed_out = speed_centroid[frame_mask]
    time_out = time[frame_mask]
    frame_out = frame[frame_mask]
    missing_out = missing[frame_mask]

    fps = frame_out[-1] / time_out[-1]

    df = pd.DataFrame({
        'X': X_out,
        'Y': Y_out,
        'Time': time_out,
        'Speed_Centroid': speed_out,
        'Missing': missing_out,
        'fps': fps
    })

    df.to_csv(output_file, index=False)