# iclue-scripts
General repo for all work done for ICLUE from SP25-present

## Project Directory

```
iclue-scripts/
├── main.py                          # Main execution script for testing multiple k values
├── README.md                        # Project documentation
├── requirements.txt                 # Python dependencies
├── anti-diagonal-code/              # Anti-diagonal matrix analysis
│   ├── anti_diag_final_var.ipynb   # Jupyter notebook for anti-diagonal computations
│   └── README.md                    # Documentation for anti-diagonal methods
├── low_rows/                        # Core matrix operations and RSK algorithm
│   ├── low_matrices.py             # Original matrix generation and analysis functions
│   ├── sum_low_matrices.py         # Optimized matrix summation operations
│   ├── v2_lr.py                    # Refactored and optimized version with RSK implementation
│   ├── data/                       # Generated search results and matrix data
│   │   ├── 2x2_search.txt          # Search results for 2x2 matrices
│   │   ├── 2x2leq10_search.txt     # Results with constraints ≤ 10
│   │   └── 2x2leq20search.txt      # Results with constraints ≤ 20
│   └── highest_weight_data/        # Highest weight vector analysis
│       ├── 2x2_new_definition.txt  # 2x2 results with new definition
│       ├── 2x2_new_MPQ.txt         # 2x2 MPQ tableaux results
│       ├── 2x2_old_MPQ.txt         # 2x2 legacy MPQ results
│       ├── 3x3_new_definition.txt  # 3x3 results with new definition
│       ├── 3x3_new_MPQ.txt         # 3x3 MPQ tableaux results
│       └── 3x3_old_MPQ.txt         # 3x3 legacy MPQ results
└── rsk-linear-operator-code/       # RSK algorithm and linear operator implementations
    ├── Algo7_1_draft.ipynb        # Draft implementation of Algorithm 7.1
    ├── algo7_1_draft.py           # Python version of Algorithm 7.1 draft
    ├── algo7_1.py                 # Final Algorithm 7.1 implementation
    ├── README.md                  # Documentation for RSK implementations
    └── rsk_linear_op_mnd_var.py   # RSK linear operator with MND variations
```

## Quick Setup

To run the code:
1. Clone the repository from Git
2. Create a new Python environment (recommended)
3. Install all dependencies:
   ```bash
   pip install -r requirements.txt
   ```
4. Run files as needed

## Project Overview

This repository contains mathematical algorithms and analysis tools for matrix operations, particularly focusing on:

- **RSK Algorithm Implementation**: Robinson-Schensted-Knuth correspondence for tableaux
- **Matrix Bracketing Operations**: Validation and transformation of matrix structures  
- **Anti-diagonal Analysis**: Specialized computations for anti-diagonal matrix patterns
- **Highest Weight Vectors**: Analysis of weight vector properties in various matrix spaces
- **Performance Optimization**: Optimized implementations using numpy and itertools for large-scale computations
