from mersenne_linalg_untwister import MersenneSolver
import random

solver = MersenneSolver()

# Generate some outputs to feed to the solver. Here we only use the LSB of the output numbers, so we need to generate more numbers than usual.
output_nums = [
    (random.getrandbits(32) & 0xFF, 0xFF) for _ in range(624 * 5 + 620) # we need exactly 3740 bytes (rank 19945) to reverse for LSB and initial_idx=0, got this number by trial and error
]
# Generate some random numbers to compare against.
actual_next_nums = [
    random.getrandbits(32) for _ in range(624)
]

reversed_state = solver.reverse_state_by_outputs(output_nums, initial_idx=0, check_solution=True, verbose=True)
# Create a new random generator with the reversed state
r = random.Random()
r.setstate((3, tuple(reversed_state), None))
reversed_next_nums = [r.getrandbits(32) for _ in range(624)]

# Compare the actual next numbers with the reversed next numbers
print("diff nums: ", sum([a ^ b for a, b in zip(actual_next_nums, reversed_next_nums)]))