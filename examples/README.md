# Examples

This directory contains various examples of how you could use `mersenne_linalg_untwister`.

## Example 1: Twist Matrix

See [/temper_check.py](/temper_check.py)

This example shows how to use the `MersenneGF2.get_twist_matrix` function.

## Example 2: Temper Matrix

See [/temper_check.py](/temper_check.py)

This example shows how to use the `MersenneGF2.get_temper_matrix` function.

## Example 3: Basic Usage (Solve inner state by 624 full outputs)

See [/solve624ints.py](/solve624ints.py)

This example shows how to use the `MersenneSolver.reverse_state_by_outputs` function to solve the inner state of a Mersenne Twister generator by using 624 full outputs.

## Example 4: Partial Outputs

See [/solveLSB.py](/solveLSB.py)

This example shows how to use the `MersenneSolver.reverse_state_by_outputs` function to solve the inner state of a Mersenne Twister generator by using partial outputs. It only uses the least significant byte of each output but can reverse the actual state.

## Example 5: Wrong initial_idx

See [/solve624wrongIdx.py](/solve624wrongIdx.py)

This example shows how to use the `MersenneSolver.reverse_state_by_outputs` function to solve the inner state of a Mersenne Twister generator by using 624 full outputs. However, it uses the wrong starting index for the outputs.

This shows that the function works even if the starting index is wrong. It shows that the reversed state is  wrong but yields the correct outputs.