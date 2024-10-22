
# Solution Description: GC Content and GC Skew Analysis

This Python script performs GC content and GC skew analysis on a DNA sequence stored in a FASTA file. It computes these statistics in a windowed manner, sliding across the sequence, and visualizes the results using histograms and plots.

## General results
I peformed computations on the reference genome of *Eschericha coli* Strain: K-12 substrain: MG1655 stored in the NCBI database:  
[NCBI Datasets](https://www.ncbi.nlm.nih.gov/datasets/taxonomy/511145/)  
For final results I have used the following parameters:  
* window_size = 10000
* step = 1000  

The histogram of GC content accross windows is centered around ~0.52.
The GC-skew plot shows characteristic relation with two sharp transitions: from positive to negative and from negative to positive skew. It is well known effect raising from a mechanism of DNA replication in bacteria. Replication of a single strand starts in an origin of replication and is conducted in two directions resulting in leading strand and lagging strand. In general, guanosine is more abundant in leading strand rather than in lagging strand. This is why in leading strand GC-skew will be rather positive and leading strand will be rather negative. Transition points between negative and positive GC-skew are origin of replication and terminator of replication.

## Dependencies
The script uses the following libraries:
- `typing`: For type hints.
- `matplotlib.pyplot`: For plotting and visualizing results.
- `os`: To handle file paths.

## Functions

### `compute_GC_stats(path2file: str, window_size: int = 100, step: int = 50) -> Tuple[List[float], List[float]]`
This function reads a FASTA file containing a single DNA sequence and computes two key statistics:
1. **GC Content**: The percentage of nucleotides in a window that are either Guanine (G) or Cytosine (C).
2. **GC Skew**: Measures the relative abundance of G vs C using the formula:
   \[
   \text{GC Skew} = \frac{G - C}{G + C}
   \]
   
#### Parameters:
- `path2file` (str): Path to the input FASTA file.
- `window_size` (int): The size of the sliding window used to compute statistics (default is 100 nucleotides).
- `step` (int): The number of bases to slide the window after each computation (default is 50 nucleotides).

#### Returns:
- **gc** (`List[float]`): A list of GC content values, one per window.
- **gc_skew** (`List[float]`): A list of GC skew values, one per window.

#### Operation:
- The file is opened and read line by line, ignoring the FASTA header (`>`).
- The sequence is processed in overlapping windows. For each window, GC content and GC skew are calculated.
  
### `visualize_results(gc: List[float], gc_skew: List[float])`
This function generates visualizations for the computed GC content and GC skew statistics.

#### Parameters:
- `gc` (`List[float]`): List of GC content values.
- `gc_skew` (`List[float]`): List of GC skew values.

#### Visualization:
- **GC Content Histogram**: A histogram showing the distribution of GC content values across windows.
- **GC Skew Plot**: A line plot of GC skew values across windows with a baseline (0.0) marked by a dashed red line.
- The figures are displayed using `matplotlib`.

### `main(path2file: str, window_size: int = 10000, step: int = 1000)`
This is the main function that combines the GC computation and visualization steps. It processes a given FASTA file and calls `compute_GC_stats` to obtain the GC content and GC skew, then visualizes the results with `visualize_results`.

#### Parameters:
- `path2file` (str): Path to the FASTA file to be analyzed.
- `window_size` (int): Window size for computation (default is 10000 nucleotides).
- `step` (int): Step size for the sliding window (default is 1000 nucleotides).

### Execution:
When executed as a script, it:
1. Retrieves the current working directory.
2. Specifies a hard-coded file path to a sample FASTA file.
3. Calls the `main()` function to perform the analysis on this file.

## Example of Usage
```python
if __name__ == "__main__":
    CWD = os.getcwd()
    main(os.path.join(CWD, "ncbi_dataset/data/GCA_000005845.2/GCA_000005845.2_ASM584v2_genomic.fna"))
```

The script is designed to be run in a directory containing a specific dataset structure. You may need to update the path to point to your FASTA file location.
