from mersenne_linalg_untwister import MersenneGF2
import random
from sage.all import vector, Integer, GF
import time

mersenne = MersenneGF2()
temper_matrix  = mersenne.get_temper_matrix(verbose=True)


random.seed(1234509876)
out_nums = [random.getrandbits(32) for _ in range(624)]
after_state = random.getstate()[1][:-1] # index is now 624, still in same twist state

# convert to GF(2) vectors
after_state_vec = [vector(GF(2), Integer(m).digits(base=2,padto=32)) for m in after_state]
out_nums_vec = [vector(GF(2), Integer(m).digits(base=2,padto=32)) for m in out_nums]

# check if temper_matrix * state_int == out_num_int
start_time = time.time()
tempered_state_vec = [temper_matrix * v for v in after_state_vec]
print("Equal:", sum([a + b for a, b in zip(tempered_state_vec, out_nums_vec)]) == 0)  # + in GF(2) is XOR, we check if they are exactly equal
print("Time taken: ", time.time() - start_time)