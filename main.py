from low_rows.low_matrices import generateChunks, all_kxk_matrices, new_examples, gen_giants


if __name__ == "__main__":
    B = generateChunks(3)
    C = generateChunks(3)
    A = all_kxk_matrices(3)

    giants = gen_giants(A, B, C, 3)
    new_ex = new_examples(giants, 3, 15)
