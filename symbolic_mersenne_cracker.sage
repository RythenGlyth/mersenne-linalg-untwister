import numpy as np
from scipy.linalg import block_diag
import time
from typing import List, Tuple, Union
# from sage.all import *
import sage.all as sage
import random
from functools import reduce

GFMatFull = sage.MatrixSpace(sage.GF(2), 624*32, 624*32)
GFMatInt = sage.MatrixSpace(sage.GF(2), 32)

temper_matrix = None
full_temper_matrix = None
twist_matrix = None
inv_twist_matrix = None

class MersenneGF2:
    # twist_matrix @ (Y_i, ..., Y_{i+N-1}) = (Y_{i+N}, ... , Y_{i+2N-1})
    def get_twist_matrix(self, verbose: bool=False):
        global twist_matrix
        if twist_matrix is not None:
            return twist_matrix
        try:
            twist_matrix = sage.load("twist_matrix.sobj")
            return twist_matrix
        except FileNotFoundError:
            pass
        if(verbose): print("calculating twist matrix (this may take a while)")
        start_time = time.time()
        twist_matrix = sage.matrix(sage.GF(2), 32 * 624, 32 * 624, sparse=True)
        for i in range(0,624):
            for b in range(0,32):
                if i < 227:
                    twist_matrix[i*32 + b, (i+397) * 32 + b] += 1
                else:
                    twist_matrix[i*32 + b] += twist_matrix[(i-227) * 32 + b]
                
                if (1 << b) & 0x9908B0DF:
                    if i == 623:
                        twist_matrix[i*32 + b] += twist_matrix[0]
                    else:
                        twist_matrix[i*32 + b, (i+1) * 32] += 1

            for b in range(0,30):
                if i == 623:
                    twist_matrix[i*32 + b] += twist_matrix[b+1]
                else:
                    twist_matrix[i*32 + b, (i+1) * 32 + b + 1] += 1
            twist_matrix[i*32 + 30, i * 32 + 31] += 1
        if(verbose): print("twist matrix calculated in", time.time() - start_time, "seconds")
        sage.save(twist_matrix, "twist_matrix.sobj")
        return twist_matrix
    
    # mat_1 @ x = x ^ (x >> 11)
    def get_temper_matrix_1(self):
        mat_1 = GFMatInt()
        for i in range(32):
            mat_1[i,i] = 1
            if i+11 < 32:
                mat_1[i,i+11] = 1
        return mat_1

    # mat_2 @ x = x ^ (x << 7) & 0x9D2C5680
    def get_temper_matrix_2(self):
        mat_2 = GFMatInt()
        for i in range(32):
            mat_2[i,i] = 1
            if i >= 7 and (1 << i) & 0x9D2C5680:
                mat_2[i,i-7] = 1
        return mat_2

    # mat_3 @ x = x ^ (x << 15) & 0xEFC60000
    def get_temper_matrix_3(self):
        mat_3 = GFMatInt()
        for i in range(32):
            mat_3[i,i] = 1
            if i >= 15 and (1 << i) & 0xEFC60000:
                mat_3[i,i-15] = 1
        return mat_3

    # mat_4 @ x = x ^ (x >> 18)
    def get_temper_matrix_4(self):
        mat_4 = GFMatInt()
        for i in range(32):
            mat_4[i,i] = 1
            if i+18 < 32:
                mat_4[i,i+18] = 1
        return mat_4

    def get_temper_matrix(self, verbose: bool=False):
        global temper_matrix
        if temper_matrix is not None:
            return temper_matrix
        try:
            temper_matrix = sage.load("temper_matrix.sobj")
            return temper_matrix
        except FileNotFoundError:
            pass
        if(verbose): print("calculating temper matrix (this may take a while)")
        start_time = time.time()
        temper_matrix = self.get_temper_matrix_4() * self.get_temper_matrix_3() * self.get_temper_matrix_2() * self.get_temper_matrix_1()
        if(verbose): print("temper matrix calculated in", time.time() - start_time, "seconds")
        sage.save(temper_matrix, "temper_matrix.sobj")
        return temper_matrix
    
    def get_full_temper_matrix(self, verbose: bool=False):
        global full_temper_matrix
        if full_temper_matrix is not None:
            return full_temper_matrix
        if(verbose): print("calculating full temper matrix (this may take a while)")
        start_time = time.time()
        full_temper_matrix = sage.block_diagonal_matrix(
            *([self.get_temper_matrix()] * 624)
        )
        if(verbose): print("full temper matrix calculated in", time.time() - start_time, "seconds")
        return full_temper_matrix
    
class MersenneSolver:
    mersenne = MersenneGF2()

    
    def reverse_state_by_outputs(self, outputs: List[Union[int, Tuple[int, int]]], initial_idx: int=0, check_solution: bool=True, verbose: bool=False) -> list[int]:
        """
        Reverse the state of the Mersenne Twister given outputs, where some
        outputs might be fully known and others only partially known.

        Args:
            outputs: A list where each element corresponds to a consecutive output slot.
                     The length of this list, k, determines the final state index.
                     Each element can be one of two types:
                     - An int: Represents a fully known 32-bit output value (mask 0xFFFFFFFF).
                     - A Tuple[int, int]: Represents (value, mask) for a partially
                       known output. `mask` indicates known bits (1=known, 0=unknown).
            initial_idx: The index of the first output in the list. This is the index
                       of the first output in the Mersenne Twister's internal state. 
                       As it seems it's completely unnecessary to use it, it defaults to 0.
                       If the index is wrong, the function returns a wrong state but the state
                       will produce the same outputs as the actual state. (Needs confirmation)
            check_solution: If True, the function checks if the linear system is constrained enough.
                            If False, it will not check the linear system and will return a possible
                            state which may not be able to predict further outputs.
            verbose: If True, the function will print the progress and time taken for each step.

        Returns:
            A tuple of 625 integers representing the reconstructed internal state
            MT[0]...MT[623] are the state vector. MT[624] is the index of the next value to be returned.
        """
        initial_idx = initial_idx % 624

        twist_matrix = self.mersenne.get_twist_matrix(verbose=verbose)
        temper_matrix = self.mersenne.get_temper_matrix(verbose=verbose)

        ALLOC = len(outputs) * 32 # TODO: make it smaller (currently overhead when mask is not 0xFFFFFFFF for all elements)

        final_matrix_rows = sage.matrix(sage.GF(2), ALLOC, 32 * 624)
        final_matrix_rows_twists = []
        final_matrix_rows_i = 0
        segment_start_row = 0
        right_side = sage.vector(sage.GF(2), ALLOC)
        
        if(verbose): print("Creating linear system...")
        start_time = time.time()

        curr_power = initial_idx // 624
        i = initial_idx % 624
        for out_idx, item in enumerate(outputs):
            if i == 624:
                final_matrix_rows_twists.append((segment_start_row, final_matrix_rows_i, curr_power))
                segment_start_row = final_matrix_rows_i
                i = 0
                curr_power += 1
            
            value: int
            mask: int

            if isinstance(item, int):
                value = item
                mask = 0xFFFFFFFF
            elif isinstance(item, tuple) and len(item) == 2 and all(isinstance(x, int) for x in item):
                value, mask = item
            elif isinstance(item, list) and len(item) == 2 and all(isinstance(x, int) for x in item):
                value, mask = item
            else:
                raise TypeError(
                    f"Element at index {out_idx} in 'outputs' must be an int or a Tuple[int, int], "
                    f"but got {type(item)}"
                )

            mask &= 0xFFFFFFFF # Ensure 32-bit mask

            for b in range(32):
                if (mask >> b) & 1:
                    final_matrix_rows[final_matrix_rows_i, i*32:(i+1)*32] = temper_matrix[b]
                    right_side[final_matrix_rows_i] = (value >> b) & 1
                    final_matrix_rows_i += 1

            i += 1


        final_matrix_rows_twists.append((segment_start_row, final_matrix_rows_i, curr_power)) # append the last twist
        if(verbose): print("Transforming array into GF(2) matrix and twisting...")

        right_side = right_side[:final_matrix_rows_i]
        
        final_matrix = sage.matrix(sage.GF(2), final_matrix_rows_i, 32 * 624, sparse=True)


        # ------- This is somehow slower --------
        # for start_row, end_row, twists in final_matrix_rows_twists:
        #     if end_row <= start_row: # Skip empty segments
        #         continue
        #     final_matrix.set_block(
        #         start_row, 0, final_matrix_rows[start_row:end_row] * (twist_matrix ** twists)
        #     )
        # ----------------------------------------
        # This is faster
        final_matrix = reduce(
            lambda x,y: x.stack(y),
            [
                (final_matrix_rows[
                    start_row
                    :end_row
                ] * (twist_matrix ** twists))
                for (start_row,end_row,twists) in final_matrix_rows_twists
            ]
        )
        # ----------------------------------------



        if(verbose): print("Linear system created in", time.time() - start_time, "seconds")

        # Explaination:
        # Here, we check if Ker(final_matrix) \subset Ker(final_twist_pow) because this means that for every solution to final_matrix*x=right_side, final_twist_pow*x is the same.
        # Given this condition the final output (final_twist_pow*x) is determined.
        # This is equivalent to Ker(final_matrix) = Ker(final_matrix) \cap Ker(final_twist_pow)
        # Which is equivalent to rank(final_matrix) = rank(final_matrix.stack(final_twist_pow)) - which is more efficient to check.
        if(check_solution):
            if(verbose): print("Checking if the linear system is constrained enough...")
            final_matrix_rank = final_matrix.rank()
            stack_rank = final_matrix.stack(twist_matrix ** (curr_power)).rank()
            if(final_matrix_rank != stack_rank):
                raise ValueError(
                    "The linear system is not constrained enough. "
                    "Please provide more outputs or check the inputs."
                )
            if(verbose): print("Linear system is constrained enough")

        # solve the linear system
        if(verbose): print("Solving linear system...")
        start_time = time.time() 

        x_p = final_matrix.solve_right(right_side)

        if(verbose): print("extracted particular solution after", time.time() - start_time, "seconds")

        orig_state_transformed = (twist_matrix ** (curr_power)) * x_p

        if(verbose): print("Solved in", time.time() - start_time, "seconds")
        # convert to list of integers
        orig_state_nums = [
            int(sage.ZZ(orig_state_transformed[i*32:(i+1)*32].list(), base=2)) for i in range(len(orig_state_transformed)//32)
        ]

        return orig_state_nums + [i]  # add the index of the next value to be returned