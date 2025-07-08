from low_rows.low_matrices import (
    generateChunks,
    all_kxk_matrices,
    new_examples,
    gen_giants,
    gen_giants_AB_var,
    gen_giants_AC_var,
    masked_examples,
    search_masked_examples,
    print_side_by_side,
    border_check,
    bracketing, 
    gen_ABC0_giants,
    filter_ABC0_by_bracketing,
    filter_new_ABC0_by_bracketing
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

# def test_old_search(k, limit=None):
#     A, B, C = ABC(k, limit)

#     ab_giants = gen_giants_AB_var(A, B, k)
#     ac_giants = gen_giants_AC_var(A, C, k)

#     return search_masked_examples(ab_giants, ac_giants, k)

# def test_new_search(k, limit=None):
#     A, B, C = filtered_ABC(k, limit)

#     ab_giants = gen_giants_AB_var(A, B, k)
#     ac_giants = gen_giants_AC_var(A, C, k)

    # return search_masked_examples(ab_giants, ac_giants, k, limit)

def find_matching_pairs(list1, list2):
    combined = list1 + list2
    seen = {}
    matches = []

    for a, b, c, p_ac, q_ab in combined:
        key = (a.tobytes(), b.tobytes()) # rmb this makes np arrays hashable
        if key in seen:
            matches.append((seen[key], (c, p_ac, q_ab))) 
        else:
            seen[key] = (c, p_ac, q_ab)

    return matches

def find_matching_pairs(list1):
    combined = list1
    seen = {}
    matches = []

    for a, b, c, p_ac, q_ab in combined:
        key = (a.tobytes(), b.tobytes()) # rmb this makes np arrays hashable
        if key in seen:
            matches.append((seen[key], (c, p_ac, q_ab))) 
        else:
            seen[key] = (c, p_ac, q_ab)

    return matches

def test_old_search(A, B, C, k, limit=100000000):
    abc0_giants = gen_ABC0_giants(A, B, C, k, max_count=limit)
    ab_pairs, ac_pairs = filter_ABC0_by_bracketing(abc0_giants, k)

    ab_giants = gen_giants_AB_var([a for a, _ in ab_pairs], [b for _, b in ab_pairs], k)
    ac_giants = gen_giants_AC_var([a for a, _ in ac_pairs], [c for _, c in ac_pairs], k)

    result = search_masked_examples(ab_giants, ac_giants, k, limit)
    return result


def test_new_search(A, B, C, k, limit=10000000):
    abc0_giants = gen_ABC0_giants(A, B, C, k, max_count=limit)
    ab_pairs, ac_pairs = filter_new_ABC0_by_bracketing(abc0_giants, k)

    ab_giants = gen_giants_AB_var([a for a, _ in ab_pairs], [b for _, b in ab_pairs], k)
    ac_giants = gen_giants_AC_var([a for a, _ in ac_pairs], [c for _, c in ac_pairs], k)

    result = search_masked_examples(ab_giants, ac_giants, k, limit)
    return result

if __name__ == "__main__":
    A, B, C = ABC(2)
    old_list = test_old_search(A, B, C, 2)
    A, B, C = filtered_ABC(2)
    new_list = test_new_search(A, B, C, 2)
    new_list.extend(old_list)
    
    matches = find_matching_pairs(new_list)

    for one, two in matches:
        if not (one[0] == two[0]).all():
            print("Match found (1, 2):")
            c1, p1, q1 = one
            c2, p2, q2 = two
            # b1 = bracketing(c1, 2)
            # b2 = bracketing(c2, 2)

            print_side_by_side(c1, c2)
            print("1: P, Q")
            print_side_by_side(p1, q1)
            print("2: P, Q")
            print_side_by_side(p2, q2)
            # print("brackets 1")
            # print(b1)
            # print("brackets 2")
            # print(b2)
            print()
    print(sum)


