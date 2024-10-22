import math
import matplotlib.pyplot as plt
import numpy as np

# Ma wyjść pik z długim ogonem. Szczyt około 0.13.

def pJC(a:str, b:str, t:float):
	n, m = _compute_n_m(a, b)
	gt = 1 + 3 * math.e**(-4*t)
	ht = 1 - math.e**(-4*t)
	return (1/4 * gt)**n * (1/4 * ht)**m


def _compute_n_m(seq1, seq2, nucleobases=set('ACGT')):
	# percentage of nucleotide difference in two sequences

	diff_count = 0 # number of nucleotide differences
	valid_nucleotides_count = 0 # number of valid nucleotides (value is float for computing percentage)

	for a, b in zip(seq1, seq2):
		if a in nucleobases and b in nucleobases:
			valid_nucleotides_count += 1
			if a != b: diff_count += 1
	
	return valid_nucleotides_count - diff_count, diff_count


def show_plot(results):
	plt.figure(figsize=(5,5))
	plt.plot(np.linspace(0, 4, num=50), results)
	plt.title("JC plot")
	plt.show()


if __name__ == "__main__":
	seqA = "ACCATAACGA-TGCATCGGA-GACACAAACACGGGGAAACGAGA"
	seqB = "ACCAT--CGC-TCCTTAGGAG---ACAATCTCTGGGAACAGGA-"
	
	res = []
	for t in np.linspace(0, 4, num=50):
		partial = pJC(seqA, seqB, t)
		res.append(partial)
	show_plot(res)
		

		