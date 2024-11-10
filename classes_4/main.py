import numpy as np
import matplotlib.pyplot as plt
from collections import Counter


NUCLEOTIDES = ['A', 'C', 'G', 'T']

def jc69_transition_probs(alpha, d):
    # Compute probabilities based on JC69 model
    p = 1/4 + 3/4 * np.exp(-4 * alpha * d / 3)
    q = (1 - p) / 3
    return p, q

def simulate_step(sequence, p, q):
    # Generate a new sequence by simulating JC69 substitutions
    new_sequence = []
    for nucleotide in sequence:
        if np.random.rand() < p:
            # No substitution
            new_sequence.append(nucleotide)
        else:
            # Substitute with one of the other nucleotides
            possible_substitutions = [nuc for nuc in NUCLEOTIDES if nuc != nucleotide]
            new_sequence.append(np.random.choice(possible_substitutions))
    return ''.join(new_sequence)

def nucleotide_frequencies(sequence):
    # Calculate the frequency of each nucleotide
    counts = Counter(sequence)
    total = len(sequence)
    return np.array([counts[nuc] / total for nuc in NUCLEOTIDES])

def jc69_simulator(S, alpha, d, K, metric='euclidean'):
    p, q = jc69_transition_probs(alpha, d)
    current_sequence = S
    differences = []
    stationary_distribution = np.array([0.25, 0.25, 0.25, 0.25])
    
    for _ in range(K):
        current_sequence = simulate_step(current_sequence, p, q)
        freqs = nucleotide_frequencies(current_sequence)
        if metric == 'euclidean':
            difference = np.linalg.norm(freqs - stationary_distribution)
        elif metric == 'cityblock':
            difference = np.sum(np.abs(freqs - stationary_distribution))
        differences.append(difference)
    
    return differences


def main():
    # Define initial parameters
    S="ATGGGCACTGCTGGAAAAGTTATTAAGTGCAAAGCAGCTGTGCTTTGGGAGCAGAAGCAACCCTTCTCCATTGAGGAAATAGAAGTTGCCCCACCAAAGACTAAAGAAGTTCGCATTAAGATTTTGGCCACAGGAATCTGTCGCACAGATGACCATGTGATAAAAGGAACAATGGTGTCCAAGTTTCCAGTGATTGTGGGACATGAGGCAACTGGGATTGTAGAGAGCATTGGAGAAGGAGTGACTACAGTGAAACCAGGTGACAAAGTCATCCCTCTCTTTCTGCCACAATGTAGAGAATGCAATGCTTGTCGCAACCCAGATGGCAACCTTTGCATTAGGAGCGATATTACTGGTCGTGGAGTACTGGCTGATGGCACCACCAGATTTACATGCAAGGGCAAACCAGTCCACCACTTCATGAACACCAGTACATTTACCGAGTACACAGTGGTGGATGAATCTTCTGTTGCTAAGATTGATGATGCAGCTCCTCCTGAGAAAGTCTGTTTAATTGGCTGTGGGTTTTCCACTGGATATGGCGCTGCTGTTAAAACTGGCAAGGTCAAACCTGGTTCCACTTGCGTCGTCTTTGGCCTGGGAGGAGTTGGCCTGTCAGTCATCATGGGCTGTAAGTCAGCTGGTGCATCTAGGATCATTGGGATTGACCTCAACAAAGACAAATTTGAGAAGGCCATGGCTGTAGGTGCCACTGAGTGTATCAGTCCCAAGGACTCTACCAAACCCATCAGTGAGGTGCTGTCAGAAATGACAGGCAACAACGTGGGATACACCTTTGAAGTTATTGGGCATCTTGAAACCATGATTGATGCCCTGGCATCCTGCCACATGAACTATGGGACCAGCGTGGTTGTAGGAGTTCCTCCATCAGCCAAGATGCTCACCTATGACCCGATGTTGCTCTTCACTGGACGCACATGGAAGGGATGTGTCTTTGGAGGTTTGAAAAGCAGAGATGATGTCCCAAAACTAGTGACTGAGTTCCTGGCAAAGAAATTTGACCTGGACCAGTTGATAACTCATGTTTTACCATTTAAAAAAATCAGTGAAGGATTTGAGCTGCTCAATTCAGGACAAAGCATTCGAACGGTCCTGACGTTTTGA"
    #S="AAAAAAAAAAAAAAAATTTTGGC"*50
    print("BASELINE:", nucleotide_frequencies(S))

    alpha = 0.1
    d = 0.05
    K = 500

    # Base simulation
    differences = jc69_simulator(S, alpha, d, K)

    # Simulation with rescaled parameters: d --> d/2 and double K
    d_rescaled = d / 2
    K_rescaled = 2 * K
    differences_rescaled = jc69_simulator(S, alpha, d_rescaled, K_rescaled)

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # Plot the original simulation result
    axes[0].plot(range(1, K + 1), differences, label="Original d, K")
    axes[0].set_title("Convergence to Stationary Distribution (Original)")
    axes[0].set_xlabel("Simulation Steps")
    axes[0].set_ylabel("Difference from Stationary Distribution")
    axes[0].legend()

    # Plot the rescaled simulation result
    axes[1].plot(range(1, K_rescaled + 1), differences_rescaled, label="Rescaled d and K", color='orange')
    axes[1].set_title("Convergence to Stationary Distribution (Rescaled)")
    axes[1].set_xlabel("Simulation Steps")
    axes[1].set_ylabel("Difference from Stationary Distribution")
    axes[1].legend()

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()