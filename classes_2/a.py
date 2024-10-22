import math


def pK(a:str, b:str, t:float): #a, b -> sequences (str), t -> time (float)
    ...
    return



#wzor  chyba xD
# g(t) = 1+3e^(-4at)
# h(t) = 1 - e^(-4at)
# f(t) = (1/4 * g(t))^n * (1/4 * h(t))^m
# n = pozycje identyczne, n != 0
# m = pozycje różne, m != 0




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



if __name__ == "__main__":
	res = pJC("ACT-AA","AGGG-A",10)
	print(res)