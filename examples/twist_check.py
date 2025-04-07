from mersenne_linalg_untwister import MersenneGF2
import random
from sage.all import vector, Integer, GF
import time

mersenne = MersenneGF2()
twist_matrix  = mersenne.get_twist_matrix(verbose=True)


random.seed(1234509876)
before_state = random.getstate()[1][:-1] # index is now 624, twist occurs after next call to random.getrandbits
random.getrandbits(1) # this will cause the twist to occur
after_state = random.getstate()[1][:-1] # index is now 1, twist has occurred

# convert to GF(2) vectors
before_state_vec = vector(GF(2), [b for m in before_state for b in Integer(m).digits(base=2,padto=32)])
after_state_vec = vector(GF(2), [b for m in after_state for b in Integer(m).digits(base=2,padto=32)])

# check if twist_matrix * before_state_vec == after_state_vec
start_time = time.time()
print("Equal:", twist_matrix * before_state_vec == after_state_vec)
print("Time taken: ", time.time() - start_time)