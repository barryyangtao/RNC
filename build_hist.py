## histogram plotting function
## Barry Y. Li and Tim Duong (2024)

import sys
import numpy as np

def plot_histogram(filename):
    # read data from the file
    try:
        data = np.loadtxt(filename)
    except FileNotFoundError:
        print(f"Error: File '{filename}' not found.")
        return
    
    # calculate histogram
    hist, bin_edges = np.histogram(data, bins=100)

    # find maximum count for scaling
    max_count = np.max(hist)

    # display histogram in terminal (vertical)
    print("Histogram of data:")
    for i in range(len(hist)):
        bin_str = f"{bin_edges[i]:.2f} - {bin_edges[i+1]:.2f}"
        count_str = "#" * int(hist[i] * 50 / max_count)  # Scale the counts to fit in 50 chara.
        print(f"{bin_str.ljust(15)}| {count_str}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python3 histo.py <filename>")
    else:
        filename = sys.argv[1]
        plot_histogram(filename)
