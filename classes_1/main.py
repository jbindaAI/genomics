from typing import List, Tuple
import matplotlib.pyplot as plt
import os


def compute_GC_stats(path2file:str, window_size:int=100, step:int=50) -> Tuple[List[float], List[float]]:
    """
    Read fasta file with ONE sequence, computes statistics: GC content and GC skew.
    Statistics are computed window-wise.
    window_size: specifies length of a sliding window.
    step: specifies sliding step of the window.

    Returns: Tuple of window-wise computed statistics    
    """
    seq = ""
    gc = []
    gc_skew = []
    with open(path2file, "r") as file:
        for line in file:
            if line[0] == ">":
                continue
            seq+=line.strip()
    for start in range(0, len(seq), step):
        if start+window_size < len(seq):
            window = seq[start : start+window_size]
            counts_c = window.count("C")
            counts_g = window.count("G")
            gc.append((counts_c+counts_g)/window_size)
            gc_skew.append((counts_g-counts_c)/(counts_c+counts_g))
        else:
            continue
    return (gc, gc_skew)

def visualize_results(gc:List[float], gc_skew:List[float]):
    fig, axs = plt.subplots(1, 2, figsize=(12, 5))
    axs[0].hist(gc, bins=35)
    axs[0].set_title("Histogram of GC content")
    axs[0].set_xlabel("GC content")
    axs[0].set_ylabel("counts")

    axs[1].plot(gc_skew)
    axs[1].set_title("GC-skew over windows")
    axs[1].set_xlabel("Windows")
    axs[1].set_ylabel("GC-skew")
    axs[1].axhline(y = 0.0, color = 'r', linestyle = '--', label='baseline')
    axs[1].legend(loc="lower right")

    fig.suptitle("Results")
    plt.tight_layout
    plt.savefig("Figure")
    plt.show()

def main(path2file:str, window_size:int=10000, step:int=1000):
    visualize_results(*compute_GC_stats(path2file, window_size, step))

if __name__ == "__main__":
    CWD = os.getcwd()
    main(os.path.join(CWD, "ncbi_dataset/data/GCA_000005845.2/GCA_000005845.2_ASM584v2_genomic.fna"))