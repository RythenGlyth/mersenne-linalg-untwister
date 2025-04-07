from mersenne_linalg_untwister import MersenneSolver
import random

solver = MersenneSolver()

# Create some offset which is different to the initial_idx
random.getrandbits(32)
random.getrandbits(32)
random.getrandbits(32)
random.getrandbits(32)

# Generate some outputs to feed to the solver.
output_nums = [
    random.getrandbits(32) for _ in range(624)
]

# Extract the state of the random generator
actual_state = random.getstate()[1]
# Generate some random numbers to compare against.
actual_next_nums = [
    random.getrandbits(32) for _ in range(624)
]

# Reverse the state using the solver, use an initial_idx of 200 to simulate a wrong index
reversed_state = solver.reverse_state_by_outputs(output_nums, initial_idx=200, check_solution=True, verbose=True)
r = random.Random()
r.setstate((3, tuple(reversed_state), None))
reversed_next_nums = [r.getrandbits(32) for _ in range(624)]

# Compare the actual next numbers with the reversed next numbers
print("diff nums: ", sum([a ^ b for a, b in zip(actual_next_nums, reversed_next_nums)]))
# Compare the actual state with the reversed state
print("diff state: ", sum([a ^ b for a, b in zip(actual_state[:-1], reversed_state[:-1])]))
# Compare the indices of the actual state with the reversed state
print("state_idx: actual: ", actual_state[-1], "reversed: ", reversed_state[-1])