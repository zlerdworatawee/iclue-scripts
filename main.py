from low_rows.low_matrices import (
    generateChunks,
    all_kxk_matrices,
    new_examples,
    gen_giants,
    gen_giants_AB_var,
    gen_giants_AC_var,
    masked_examples,
    search_masked_examples,
    testing,
)

def ABC(k, limit=None):
    B = generateChunks(k)
    C = generateChunks(k)
    A = list(all_kxk_matrices(k))

    return A, B, C

def filtered_ABC(k, limit=None):
    A, B, C = ABC(k, limit)
    giants = gen_giants(A, B, C, k)
    new_ex = new_examples(giants, k, limit)

    A = [v[0] for v in new_ex.values()]
    B = [v[1] for v in new_ex.values()]
    C = [v[2] for v in new_ex.values()]
    return A, B, C

def test_MPQ(k, limit=None):
    A, B, C = ABC(k, limit)

    ab_giants = gen_giants_AB_var(A, B, k)
    ac_giants = gen_giants_AC_var(A, C, k)

    masked_examples(ab_giants, ac_giants, k, limit)

def test_new_MPQ(k, limit=None):
    A, B, C = filtered_ABC(k, limit)

    ab_giants = gen_giants_AB_var(A, B, k)
    ac_giants = gen_giants_AC_var(A, C, k)

    masked_examples(ab_giants, ac_giants, k, limit)

def test_old_search(k, limit=None):
    A, B, C = ABC(k, limit)

    ab_giants = gen_giants_AB_var(A, B, k)
    ac_giants = gen_giants_AC_var(A, C, k)

    search_masked_examples(ab_giants, ac_giants, k)

def test_new_search(k, limit=None):
    A, B, C = filtered_ABC(k, limit)

    ab_giants = gen_giants_AB_var(A, B, k)
    ac_giants = gen_giants_AC_var(A, C, k)

    search_masked_examples(ab_giants, ac_giants, k)

def t(k, limit=None):
    A, B, C = filtered_ABC(k, limit)

    ab_giants = gen_giants_AB_var(A, B, k)
    ac_giants = gen_giants_AC_var(A, C, k)

    testing(ab_giants, ac_giants, k)

if __name__ == "__main__":
    
    t(2)
