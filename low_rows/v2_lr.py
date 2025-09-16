from IPython.display import display, Math
import numpy as np
from sympy import symbols, sympify, Matrix, Eq, solve, expand, Abs, pprint, factor
from itertools import permutations, product, chain, repeat
import time
#==============================================================================
# RSK ALGORITHM FUNCTIONS
#==============================================================================

def rsk_insert(tableau, x):
    tableau = tableau.copy()
    rows, cols = tableau.shape
    bumped = x
    for r in range(rows):
        row = tableau[r]
        mask = (row > 0)
        eligible = row[mask]
        idx = np.where(eligible > bumped)[0]
        if idx.size == 0:
            insert_pos = np.sum(mask)
            if insert_pos < cols:
                tableau[r, insert_pos] = bumped
                return tableau, (r, insert_pos)
            else:
                continue
        else:
            i = idx[0]
            bumped, tableau[r, i] = tableau[r, i], bumped
    empty_row = np.zeros(cols, dtype=int)
    empty_row[0] = bumped
    tableau = np.vstack([tableau, empty_row])
    return tableau, (tableau.shape[0] - 1, 0)

def viennot_rsk(biword):
    n = len(biword)
    P = np.zeros((n, n), dtype=int)
    Q = np.zeros((n, n), dtype=int)
    for a, b in biword:
        P, (r, c) = rsk_insert(P, a)
        Q[r, c] = b
    return P, Q

def biword(M: np.array) -> list:
    biword_P = []
    biword_Q = []
    for j in range(M.shape[1]):
        for i in range(M.shape[0]):
            cell_value = int(M[i, j])
            if cell_value != 0:
                biword_P.extend([i + 1] * cell_value)
                biword_Q.extend([j + 1] * cell_value)
    return list(zip(biword_P, biword_Q))

#==============================================================================
# UTILITY FUNCTIONS
#==============================================================================

def print_side_by_side(A, B, sep=' | '):
    for row_a, row_b in zip(A, B):
        row_a_str = ' '.join(str(int(x)) for x in row_a)
        row_b_str = ' '.join(str(int(x)) for x in row_b)
        print(f"{row_a_str}{sep}{row_b_str}")

def masking(arr, k):
    arr[arr <= k] = 0
    return arr

def check_matrix(mat):
    vals = []
    for row in mat:
        nonzero_mask = row != 0
        if not np.any(nonzero_mask):
            continue
        
        nz = row[nonzero_mask]
        if not np.all(nz == nz[0]):
            return False
        vals.append(nz[0])
    
    return len(vals) <= 1 or np.all(np.diff(vals) > 0)

#==============================================================================
# ABC0 GENERATION
#==============================================================================

# Flow of functions:
#   Create A from all_kxk_matrices
#   Create B/C from generateChunks
#   Create all ABC0/AB00/A0C0 Matrices for a given k with gen_giants

def generateChunks(k: int):
    """
    Generates B/C matrices according to the rules
    """
    final = []
    base = np.eye(k)[::-1]
    bases = [base * i for i in range(0, k + 1)]

    fill_positions = [(i, j) for i in range(k) for j in range(k) if i + j < k - 1]
    combos = list(product(range(k + 1), repeat=len(fill_positions)))

    for b in bases:
        for c in combos:
            mat = np.zeros((k, k), dtype=int)
            for idx, (i, j) in enumerate(fill_positions):
                mat[i, j] = c[idx]
            mat = mat + b
            final.append(mat)
    return final

def all_kxk_matrices(k: int):
    """
    Bashes out all kxk combinations for a range of 0 to k (can be changed later)
    """
    vals = range(k + 1)
    return [np.array(combo, dtype=int).reshape((k, k)) for combo in product(vals, repeat=k * k)]

def create_giant(a, b, c, k):
    """
    from kxk matrices A, B, C creates a ABC0 2k x 2k matrix
    """
    giant = np.zeros((2 * k, 2 * k), dtype=a.dtype)
    giant[:k, :k] = a
    giant[:k, k:] = b
    giant[k:, :k] = c
    return giant

def gen_giants(A, B, C, k, max_count=float('inf')):
    """
    Returns dictionary of ABC0/giants created 
    """
    giants = {}
    count = 0
    
    for a in A:
        if not np.any(a):
            continue
        for b in B:
            if not np.any(b):
                continue
            for c in C:
                if not np.any(c):
                    continue
                    
                giants[count] = [a, b, c, create_giant(a, b, c, k)]
                count += 1
                
                if count == max_count:
                    return giants
    return giants

def gen_giants_AB_var(A, B, k, max_count=float('inf')):
    """
    Generates matrices of form AB00
    """
    giants = {}
    count = 0
    c = np.zeros((k, k), dtype=int)
    
    for a in A:
        if not np.any(a):
            continue
        for b in B:
            if not np.any(b):
                continue
            giants[count] = [a, b, c, create_giant(a, b, c, k)]
            count += 1
            
            if count == max_count:
                return giants
    return giants

def gen_giants_AC_var(A, C, k, max_count=float('inf')):
    """
    Generates matrices of form A0C0
    """
    giants = {}
    count = 0
    b = np.zeros((k, k), dtype=int)
    
    for a in A:
        if not np.any(a):
            continue
        for c in C:
            if not np.any(c):
                continue
            giants[count] = [a, b, c, create_giant(a, b, c, k)]
            count += 1
            
            if count == max_count:
                return giants
    return giants

#==============================================================================
# BRACKETING FUNCTIONS
#==============================================================================
def bracketing(vec1, vec2) -> list:
    """
    Given two vectors (row or column), stack them and generate a flat bracket list:
    - ']' for counts from vec1 (first)
    - '[' for counts from vec2 (second)
    """
    if vec1.ndim == 1:
        vec1 = vec1.reshape(1, -1)
    if vec2.ndim == 1:
        vec2 = vec2.reshape(1, -1)

    if vec1.shape[0] == 1 and vec2.shape[0] == 1:
        # Row vectors
        bracket_parts = []
        for j in range(vec1.shape[1]):
            v1_count = int(vec1[0, j])
            v2_count = int(vec2[0, j])
            if v1_count > 0:
                bracket_parts.append(repeat(']', v1_count))
            if v2_count > 0:
                bracket_parts.append(repeat('[', v2_count))

        return list(chain.from_iterable(bracket_parts))
        
    elif vec1.shape[1] == 1 and vec2.shape[1] == 1:
        # Column vectors 
        bracket_parts = []
        for i in range(vec1.shape[0]):
            v1_count = int(vec1[i, 0])
            v2_count = int(vec2[i, 0])
            if v1_count > 0:
                bracket_parts.append(repeat(']', v1_count))
            if v2_count > 0:
                bracket_parts.append(repeat('[', v2_count))

        return list(chain.from_iterable(bracket_parts))
    return []

def clean_brackets(row1, row2):
    """
    Function for checking if there are any unpaired '['
    """
    brackets = bracketing(row1, row2)
    brackets = np.array(brackets)

    while True:
        if len(brackets) == 0 or np.count_nonzero(brackets == '[') == 0:
            return True, brackets
        original = brackets.copy()

        left_mask = np.where(brackets == '[')[0]
        if len(left_mask) > 0:
            first_open_idx = left_mask[0]
            if first_open_idx > 0 and np.all(brackets[:first_open_idx] == ']'):
                brackets = brackets[first_open_idx:]

        right_mask = np.where(brackets == ']')[0]
        if len(right_mask) > 0:
            last_close_idx = right_mask[-1]
            if last_close_idx < len(brackets) - 1 and np.all(brackets[last_close_idx + 1:] == '['):
                brackets = brackets[:last_close_idx + 1]

        i = 0
        while i < len(brackets) - 1:
            if brackets[i] == '[' and brackets[i + 1] == ']':
                brackets = np.delete(brackets, [i, i + 1])
                break
            else:
                i += 1

        if np.array_equal(brackets, original):
            break

    return False, brackets

def relaxed_check(arr, k) -> bool:
    """
    "new" relaxed check with forgiveness on row/col = k
    """
    row_flag = True
    col_flag = True

    for i in chain(range(0, k - 1), range(k, 2 * k - 1)):
        row_flag, _ = clean_brackets(arr[i], arr[i + 1])
        
        if not row_flag:
            return False

    for j in chain(range(0, k - 1), range(k, 2 * k - 1)):
        col_flag, _ = clean_brackets(arr[:, j], arr[:, j + 1])
        
        if not col_flag:
            return False
    return True

def original_check(arr, k) -> bool:
    """
    Standard bracketing check for all consecutive pairs
    """
    row_flag = True
    col_flag = True

    for i in range(2 * k - 1):
        row_flag, _ = clean_brackets(arr[i], arr[i + 1])
        if not row_flag:
            return False

    for j in range(2 * k - 1):
        col_flag, _ = clean_brackets(arr[:, j], arr[:, j + 1])
        if not col_flag:
            return False
    return True

#==============================================================================
# UNIT TESTS
#==============================================================================

def bracket_test():
    vec1 = np.array([[2], [0], [1]])
    vec2 = np.array([[1], [4], [3]])
    
    bracket_word = bracketing(vec1, vec2)
    print(f"Original Brackets: {bracket_word}")

    cleaned = clean_brackets(vec1, vec2)
    print(f"Cleaned Brackets: {cleaned}")

def check_test(arr):
    if original_check(arr, arr.size[0]):
        print("Passed the original check")
        return
    elif relaxed_check(arr, arr.size[0]):
        print("Passed the original check")
        return
    print("Failed both checks")

def check_matrix_check():
    arr = np.array([
        [1, 1, 1, 0],
        [2, 0, 0, 0],
        [3, 3, 0, 0],
        [0, 0, 0, 0]
    ])
    return check_matrix(arr)



#==============================================================================
# MAIN EXECUTION
#==============================================================================

def ABC(k: int):
    """
    B is the same as C anyways its just more readable for me to do this
    """
    B = generateChunks(k)
    A = all_kxk_matrices(k)
    return A, B, B

def filter(k: int):
    # note: information we need to keep are:
    # original matrix

    filtered = []
    relaxed_filtered = []

    A, B, C = ABC(k)
    matrices = gen_giants(A=A, B=B, C=C, k=k)
    # matrices[i] = (a, b, c, ABC0)
    for key in matrices.keys():
        ABC0 = matrices[key][-1]
        if original_check(ABC0, k):
            filtered.append(ABC0)
        elif relaxed_check(ABC0, k):
            relaxed_filtered.append(ABC0)
    return filtered, relaxed_filtered

def main(k: int):
    filtered = {}
    matches = {} # want this to be empty
    A, B, C = ABC(k)
    matrices = gen_giants(A=A, B=B, C=C, k=k, max_count=1000000)
    # matrices[i] = (a, b, c, ABC0)
    for key in matrices.keys():
        ABC0 = matrices[key][-1]
        ABC_biword = biword(ABC0)
        P, Q = viennot_rsk(ABC_biword)
        # Step 1: filter for valid matrices
        if check_matrix(P) and check_matrix(Q):
            if original_check(ABC0, k) or relaxed_check(ABC0, k):
                # Step 2: Create AB00 and A0C0
                ab00_dict = gen_giants_AB_var([matrices[key][0]], [matrices[key][1]], k=k)
                a0c0_dict = gen_giants_AC_var([matrices[key][0]], [matrices[key][2]], k=k)
                ab00 = ab00_dict[0][-1] 
                a0c0 = a0c0_dict[0][-1]
                # Step 3: Extract the biword for the RSK function
                ab00_biword = biword(ab00)
                a0c0_biword = biword(a0c0)
                # Step 4: Perform Standard RSK and extract the necessary tableaus
                ac_p, _ = viennot_rsk(a0c0_biword)
                _, ab_q = viennot_rsk(ab00_biword)
                # Step 5: Mask the necessary tableaus and convert to bytes for fast hashing
                m_acp = masking(ac_p, k)
                m_abq = masking(ab_q, k)
                key = (m_acp.tobytes(), m_abq.tobytes())
                if key in filtered:
                    matches[key] =  [filtered[key], (ABC0, ac_p, ab_q, P, Q)]
                    filtered[key].append((ABC0, ac_p, ab_q, P, Q))
                else:
                    filtered[key] = [(ABC0, ac_p, ab_q, P, Q)]

    return filtered, matches

def test(k: int):
    start_time = time.time()
    filtered, matches = main(k)
    end_time = time.time()
    execution_time = end_time - start_time
    
    print(f"Execution time: {execution_time:.4f} seconds")
    print(f"Found {len(filtered)} unique groups")
    print(f"Found {len(matches)} collision groups")

    for i, (key, values) in enumerate(list(filtered.items())):
        ABC0, ac_p, ab_q, P, Q = values[0]
        print(f"ABC0 matrix:")
        print(ABC0)
        print(f"A0C0 P: ")
        print(ac_p)
        print(f"AB00 Q: ")
        print(ab_q)
        print("(P, Q):")
        print_side_by_side(P, Q)
        print("---")

if __name__ == "__main__":
    test(3)
