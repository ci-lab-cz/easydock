import pytest
from rdkit import Chem

from easydock.preparation_for_docking import pdbqt2molblock

OFFENDING_PDBQT = """MODE 1
REMARK VINA RESULT:      -9.5      0.000      0.000
REMARK SMILES Cn1nc(C2CCCN2C(=O)OC(C)(C)C)cc1C(=O)NC1CCc2ccccc2NC1=O
REMARK SMILES IDX 1 1 2 2 3 3 4 4 17 5 18 6 19 7 20 8 21 9 22 11 23 12 32 14
REMARK SMILES IDX 31 15 33 16 30 18 25 19 29 20 26 21 28 22 27 23 24 24 5 26
REMARK SMILES IDX 6 27 9 28 7 29 8 30 10 31 11 32 12 33 13 34 14 35 15 36
REMARK SMILES IDX 16 37
REMARK H PARENT 21 10 31 17
REMARK Flexibility Score: 80.00
ROOT
ATOM      1  C   UNL     1      43.287  48.724 -16.282  1.00  0.00     0.174 C 
ATOM      2  N   UNL     1      43.855  47.595 -15.582  1.00  0.00    -0.263 N 
ATOM      3  N   UNL     1      43.204  46.765 -14.751  1.00  0.00    -0.175 NA
ATOM      4  C   UNL     1      44.043  45.832 -14.283  1.00  0.00     0.086 A 
ATOM      5  C   UNL     1      45.261  46.119 -14.865  1.00  0.00     0.061 A 
ATOM      6  C   UNL     1      45.143  47.228 -15.683  1.00  0.00     0.126 A 
ENDROOT
BRANCH   6   7
ATOM      7  C   UNL     1      46.185  47.868 -16.483  1.00  0.00     0.270 C 
ATOM      8  O   UNL     1      46.279  49.114 -16.583  1.00  0.00    -0.268 OA
ATOM      9  N   UNL     1      47.134  47.078 -17.180  1.00  0.00    -0.339 N 
ATOM     10  H   UNL     1      47.113  46.015 -17.139  1.00  0.00     0.164 HD
BRANCH   9  11
ATOM     11  C   UNL     1      48.162  47.679 -17.969  1.00  0.00     0.172 C 
BRANCH  11  12
ATOM     12  C   UNL     1      48.988  48.557 -17.058  1.00  0.00     0.040 CG0
ATOM     13  G   UNL     1      49.807  49.575 -17.802  1.00  0.00     0.000 G0
ENDBRANCH  11  12
BRANCH  11  14
ATOM     14  C   UNL     1      48.970  46.709 -18.728  1.00  0.00     0.246 C 
ATOM     15  N   UNL     1      48.738  46.269 -20.061  1.00  0.00    -0.324 N 
ATOM     16  O   UNL     1      49.963  46.238 -18.103  1.00  0.00    -0.272 OA
ATOM     17  H   UNL     1      49.367  45.473 -20.360  1.00  0.00     0.170 HD
BRANCH  15  18
ATOM     18  C   UNL     1      47.791  46.746 -21.031  1.00  0.00     0.044 A 
ATOM     19  C   UNL     1      47.878  47.837 -21.873  1.00  0.00    -0.024 A 
ATOM     20  C   UNL     1      46.668  45.947 -21.058  1.00  0.00     0.026 A 
ATOM     21  C   UNL     1      46.880  48.177 -22.762  1.00  0.00     0.006 A 
ATOM     22  C   UNL     1      45.660  46.262 -21.931  1.00  0.00     0.002 A 
ATOM     23  C   UNL     1      45.755  47.357 -22.772  1.00  0.00     0.000 A 
BRANCH  19  24
ATOM     24  C   UNL     1      49.049  48.763 -21.923  1.00  0.00     0.040 CG0
ATOM     25  G   UNL     1      48.593  50.063 -21.322  1.00  0.00     0.000 G0
ENDBRANCH  19  24
ENDBRANCH  15  18
ENDBRANCH  11  14
ENDBRANCH   9  11
ENDBRANCH   6   7
BRANCH   4  26
ATOM     26  C   UNL     1      43.811  44.720 -13.359  1.00  0.00     0.139 C 
ATOM     27  C   UNL     1      44.782  44.726 -12.195  1.00  0.00     0.034 C 
ATOM     28  N   UNL     1      43.988  43.406 -13.893  1.00  0.00    -0.300 N 
ATOM     29  C   UNL     1      44.747  43.301 -11.719  1.00  0.00     0.024 C 
ATOM     30  C   UNL     1      43.998  42.510 -12.761  1.00  0.00     0.124 C 
BRANCH  28  31
ATOM     31  C   UNL     1      44.126  43.012 -15.240  1.00  0.00     0.410 C 
ATOM     32  O   UNL     1      43.229  43.229 -16.105  1.00  0.00    -0.224 OA
BRANCH  31  33
ATOM     33  O   UNL     1      45.255  42.356 -15.723  1.00  0.00    -0.444 OA
BRANCH  33  34
ATOM     34  C   UNL     1      45.369  40.927 -15.656  1.00  0.00     0.106 C 
ATOM     35  C   UNL     1      44.048  40.263 -15.443  1.00  0.00     0.056 C 
ATOM     36  C   UNL     1      45.959  40.505 -17.008  1.00  0.00     0.056 C 
ATOM     37  C   UNL     1      46.390  40.540 -14.606  1.00  0.00     0.056 C 
ENDBRANCH  33  34
ENDBRANCH  31  33
ENDBRANCH  28  31
ENDBRANCH   4  26
TORSDOF 6
ENDMDL
MODE 2
REMARK VINA RESULT:      -9.4      2.398      9.352
REMARK SMILES Cn1nc(C2CCCN2C(=O)OC(C)(C)C)cc1C(=O)NC1CCc2ccccc2NC1=O
REMARK SMILES IDX 1 1 2 2 3 3 4 4 17 5 18 6 19 7 20 8 21 9 22 11 23 12 32 14
REMARK SMILES IDX 31 15 33 16 30 18 25 19 29 20 26 21 28 22 27 23 24 24 5 26
REMARK SMILES IDX 6 27 9 28 7 29 8 30 10 31 11 32 12 33 13 34 14 35 15 36
REMARK SMILES IDX 16 37
REMARK H PARENT 21 10 31 17
REMARK Flexibility Score: 80.00
ROOT
ATOM      1  C   UNL     1      42.929  48.250 -15.760  1.00  0.00     0.174 C 
ATOM      2  N   UNL     1      44.355  48.024 -15.795  1.00  0.00    -0.263 N 
ATOM      3  N   UNL     1      45.253  48.717 -16.514  1.00  0.00    -0.175 NA
ATOM      4  C   UNL     1      46.484  48.234 -16.302  1.00  0.00     0.086 A 
ATOM      5  C   UNL     1      46.314  47.196 -15.407  1.00  0.00     0.061 A 
ATOM      6  C   UNL     1      44.975  47.064 -15.087  1.00  0.00     0.126 A 
ENDROOT
BRANCH   6   7
ATOM      7  C   UNL     1      44.359  46.096 -14.183  1.00  0.00     0.270 C 
ATOM      8  O   UNL     1      43.236  46.289 -13.662  1.00  0.00    -0.268 OA
ATOM      9  N   UNL     1      45.034  44.890 -13.865  1.00  0.00    -0.339 N 
ATOM     10  H   UNL     1      45.991  44.662 -14.271  1.00  0.00     0.164 HD
BRANCH   9  11
ATOM     11  C   UNL     1      44.455  43.927 -12.983  1.00  0.00     0.172 C 
BRANCH  11  12
ATOM     12  C   UNL     1      43.987  42.759 -13.820  1.00  0.00     0.040 CG0
ATOM     13  G   UNL     1      44.789  42.576 -15.078  1.00  0.00     0.000 G0
ENDBRANCH  11  12
BRANCH  11  14
ATOM     14  C   UNL     1      45.340  43.529 -11.875  1.00  0.00     0.246 C 
ATOM     15  N   UNL     1      44.989  43.429 -10.500  1.00  0.00    -0.324 N 
ATOM     16  O   UNL     1      46.528  43.256 -12.215  1.00  0.00    -0.272 OA
ATOM     17  H   UNL     1      45.815  43.251  -9.864  1.00  0.00     0.170 HD
BRANCH  15  18
ATOM     18  C   UNL     1      43.699  43.533  -9.875  1.00  0.00     0.044 A 
ATOM     19  C   UNL     1      42.629  42.662  -9.931  1.00  0.00    -0.024 A 
ATOM     20  C   UNL     1      43.601  44.693  -9.138  1.00  0.00     0.026 A 
ATOM     21  C   UNL     1      41.440  42.894  -9.272  1.00  0.00     0.006 A 
ATOM     22  C   UNL     1      42.430  44.950  -8.473  1.00  0.00     0.002 A 
ATOM     23  C   UNL     1      41.362  44.072  -8.534  1.00  0.00     0.000 A 
BRANCH  19  24
ATOM     24  C   UNL     1      42.634  41.380 -10.699  1.00  0.00     0.040 CG0
ATOM     25  G   UNL     1      41.741  40.430  -9.951  1.00  0.00     0.000 G0
ENDBRANCH  19  24
ENDBRANCH  15  18
ENDBRANCH  11  14
ENDBRANCH   9  11
ENDBRANCH   6   7
BRANCH   4  26
ATOM     26  C   UNL     1      47.777  48.648 -16.851  1.00  0.00     0.139 C 
ATOM     27  C   UNL     1      47.936  50.155 -16.878  1.00  0.00     0.034 C 
ATOM     28  N   UNL     1      48.027  48.294 -18.214  1.00  0.00    -0.300 N 
ATOM     29  C   UNL     1      48.999  50.364 -17.919  1.00  0.00     0.024 C 
ATOM     30  C   UNL     1      49.197  49.040 -18.614  1.00  0.00     0.124 C 
BRANCH  28  31
ATOM     31  C   UNL     1      47.313  47.403 -19.040  1.00  0.00     0.410 C 
ATOM     32  O   UNL     1      47.196  46.173 -18.769  1.00  0.00    -0.224 OA
BRANCH  31  33
ATOM     33  O   UNL     1      46.682  47.803 -20.215  1.00  0.00    -0.444 OA
BRANCH  33  34
ATOM     34  C   UNL     1      47.294  47.594 -21.496  1.00  0.00     0.106 C 
ATOM     35  C   UNL     1      46.343  47.012 -22.489  1.00  0.00     0.056 C 
ATOM     36  C   UNL     1      47.752  48.984 -21.957  1.00  0.00     0.056 C 
ATOM     37  C   UNL     1      48.549  46.762 -21.336  1.00  0.00     0.056 C 
ENDBRANCH  33  34
ENDBRANCH  31  33
ENDBRANCH  28  31
ENDBRANCH   4  26
TORSDOF 6
ENDMDL
MODE 3
REMARK VINA RESULT:      -8.8      1.976      3.446
REMARK SMILES Cn1nc(C2CCCN2C(=O)OC(C)(C)C)cc1C(=O)NC1CCc2ccccc2NC1=O
REMARK SMILES IDX 1 1 2 2 3 3 4 4 17 5 18 6 19 7 20 8 21 9 22 11 23 12 32 14
REMARK SMILES IDX 31 15 33 16 30 18 25 19 29 20 26 21 28 22 27 23 24 24 5 26
REMARK SMILES IDX 6 27 9 28 7 29 8 30 10 31 11 32 12 33 13 34 14 35 15 36
REMARK SMILES IDX 16 37
REMARK H PARENT 21 10 31 17
REMARK Flexibility Score: 80.00
ROOT
ATOM      1  C   UNL     1      44.438  47.703 -15.770  1.00  0.00     0.174 C 
ATOM      2  N   UNL     1      45.090  46.461 -15.425  1.00  0.00    -0.263 N 
ATOM      3  N   UNL     1      44.484  45.301 -15.126  1.00  0.00    -0.175 NA
ATOM      4  C   UNL     1      45.399  44.359 -14.858  1.00  0.00     0.086 A 
ATOM      5  C   UNL     1      46.618  44.990 -15.006  1.00  0.00     0.061 A 
ATOM      6  C   UNL     1      46.424  46.313 -15.363  1.00  0.00     0.126 A 
ENDROOT
BRANCH   6   7
ATOM      7  C   UNL     1      47.440  47.332 -15.617  1.00  0.00     0.270 C 
ATOM      8  O   UNL     1      47.825  48.127 -14.728  1.00  0.00    -0.268 OA
ATOM      9  N   UNL     1      48.025  47.448 -16.903  1.00  0.00    -0.339 N 
ATOM     10  H   UNL     1      47.745  46.807 -17.704  1.00  0.00     0.164 HD
BRANCH   9  11
ATOM     11  C   UNL     1      49.021  48.436 -17.178  1.00  0.00     0.172 C 
BRANCH  11  12
ATOM     12  C   UNL     1      50.057  47.807 -18.079  1.00  0.00     0.040 CG0
ATOM     13  G   UNL     1      49.753  46.374 -18.419  1.00  0.00     0.000 G0
ENDBRANCH  11  12
BRANCH  11  14
ATOM     14  C   UNL     1      48.474  49.695 -17.713  1.00  0.00     0.246 C 
ATOM     15  N   UNL     1      48.039  49.934 -19.046  1.00  0.00    -0.324 N 
ATOM     16  O   UNL     1      48.396  50.640 -16.876  1.00  0.00    -0.272 OA
ATOM     17  H   UNL     1      47.805  50.946 -19.241  1.00  0.00     0.170 HD
BRANCH  15  18
ATOM     18  C   UNL     1      47.880  49.010 -20.136  1.00  0.00     0.044 A 
ATOM     19  C   UNL     1      48.518  48.984 -21.359  1.00  0.00    -0.024 A 
ATOM     20  C   UNL     1      46.938  48.050 -19.834  1.00  0.00     0.026 A 
ATOM     21  C   UNL     1      48.261  48.026 -22.318  1.00  0.00     0.006 A 
ATOM     22  C   UNL     1      46.662  47.086 -20.768  1.00  0.00     0.002 A 
ATOM     23  C   UNL     1      47.307  47.067 -21.992  1.00  0.00     0.000 A 
BRANCH  19  24
ATOM     24  C   UNL     1      49.553  49.976 -21.783  1.00  0.00     0.040 CG0
ATOM     25  G   UNL     1      48.811  51.215 -22.202  1.00  0.00     0.000 G0
ENDBRANCH  19  24
ENDBRANCH  15  18
ENDBRANCH  11  14
ENDBRANCH   9  11
ENDBRANCH   6   7
BRANCH   4  26
ATOM     26  C   UNL     1      45.237  42.953 -14.484  1.00  0.00     0.139 C 
ATOM     27  C   UNL     1      44.576  42.138 -15.578  1.00  0.00     0.034 C 
ATOM     28  N   UNL     1      44.392  42.691 -13.361  1.00  0.00    -0.300 N 
ATOM     29  C   UNL     1      44.026  40.955 -14.833  1.00  0.00     0.024 C 
ATOM     30  C   UNL     1      44.151  41.267 -13.363  1.00  0.00     0.124 C 
BRANCH  28  31
ATOM     31  C   UNL     1      43.876  43.599 -12.413  1.00  0.00     0.410 C 
ATOM     32  O   UNL     1      42.762  44.179 -12.567  1.00  0.00    -0.224 OA
BRANCH  31  33
ATOM     33  O   UNL     1      44.548  43.927 -11.239  1.00  0.00    -0.444 OA
BRANCH  33  34
ATOM     34  C   UNL     1      43.954  43.694  -9.954  1.00  0.00     0.106 C 
ATOM     35  C   UNL     1      43.071  42.488  -9.949  1.00  0.00     0.056 C 
ATOM     36  C   UNL     1      45.134  43.486  -8.995  1.00  0.00     0.056 C 
ATOM     37  C   UNL     1      43.236  44.941  -9.484  1.00  0.00     0.056 C 
ENDBRANCH  33  34
ENDBRANCH  31  33
ENDBRANCH  28  31
ENDBRANCH   4  26
TORSDOF 6
ENDMDL
MODE 4
REMARK VINA RESULT:      -8.2      1.821      2.815
REMARK SMILES Cn1nc(C2CCCN2C(=O)OC(C)(C)C)cc1C(=O)NC1CCc2ccccc2NC1=O
REMARK SMILES IDX 1 1 2 2 3 3 4 4 17 5 18 6 19 7 20 8 21 9 22 11 23 12 32 14
REMARK SMILES IDX 31 15 33 16 30 18 25 19 29 20 26 21 28 22 27 23 24 24 5 26
REMARK SMILES IDX 6 27 9 28 7 29 8 30 10 31 11 32 12 33 13 34 14 35 15 36
REMARK SMILES IDX 16 37
REMARK H PARENT 21 10 31 17
REMARK Flexibility Score: 80.00
ROOT
ATOM      1  C   UNL     1      42.990  47.766 -15.952  1.00  0.00     0.174 C 
ATOM      2  N   UNL     1      43.534  46.869 -14.958  1.00  0.00    -0.263 N 
ATOM      3  N   UNL     1      42.835  46.164 -14.055  1.00  0.00    -0.175 NA
ATOM      4  C   UNL     1      43.665  45.438 -13.295  1.00  0.00     0.086 A 
ATOM      5  C   UNL     1      44.929  45.723 -13.770  1.00  0.00     0.061 A 
ATOM      6  C   UNL     1      44.848  46.623 -14.817  1.00  0.00     0.126 A 
ENDROOT
BRANCH   6   7
ATOM      7  C   UNL     1      45.943  47.190 -15.601  1.00  0.00     0.270 C 
ATOM      8  O   UNL     1      46.525  48.248 -15.267  1.00  0.00    -0.268 OA
ATOM      9  N   UNL     1      46.382  46.541 -16.783  1.00  0.00    -0.339 N 
ATOM     10  H   UNL     1      45.932  45.639 -17.125  1.00  0.00     0.164 HD
BRANCH   9  11
ATOM     11  C   UNL     1      47.450  47.074 -17.568  1.00  0.00     0.172 C 
BRANCH  11  12
ATOM     12  C   UNL     1      48.701  46.295 -17.232  1.00  0.00     0.040 CG0
ATOM     13  G   UNL     1      49.766  46.405 -18.286  1.00  0.00     0.000 G0
ENDBRANCH  11  12
BRANCH  11  14
ATOM     14  C   UNL     1      47.162  47.114 -19.012  1.00  0.00     0.246 C 
ATOM     15  N   UNL     1      47.461  48.177 -19.909  1.00  0.00    -0.324 N 
ATOM     16  O   UNL     1      46.590  46.083 -19.470  1.00  0.00    -0.272 OA
ATOM     17  H   UNL     1      47.049  48.045 -20.874  1.00  0.00     0.170 HD
BRANCH  15  18
ATOM     18  C   UNL     1      48.227  49.372 -19.688  1.00  0.00     0.044 A 
ATOM     19  C   UNL     1      49.032  50.069 -20.567  1.00  0.00    -0.024 A 
ATOM     20  C   UNL     1      48.081  49.819 -18.392  1.00  0.00     0.026 A 
ATOM     21  C   UNL     1      49.713  51.214 -20.211  1.00  0.00     0.006 A 
ATOM     22  C   UNL     1      48.747  50.954 -18.010  1.00  0.00     0.002 A 
ATOM     23  C   UNL     1      49.552  51.646 -18.897  1.00  0.00     0.000 A 
BRANCH  19  24
ATOM     24  C   UNL     1      49.255  49.671 -21.990  1.00  0.00     0.040 CG0
ATOM     25  G   UNL     1      49.968  50.820 -22.648  1.00  0.00     0.000 G0
ENDBRANCH  19  24
ENDBRANCH  15  18
ENDBRANCH  11  14
ENDBRANCH   9  11
ENDBRANCH   6   7
BRANCH   4  26
ATOM     26  C   UNL     1      43.385  44.521 -12.188  1.00  0.00     0.139 C 
ATOM     27  C   UNL     1      43.815  45.087 -10.849  1.00  0.00     0.034 C 
ATOM     28  N   UNL     1      44.070  43.266 -12.218  1.00  0.00    -0.300 N 
ATOM     29  C   UNL     1      43.957  43.860  -9.994  1.00  0.00     0.024 C 
ATOM     30  C   UNL     1      43.875  42.671 -10.918  1.00  0.00     0.124 C 
BRANCH  28  31
ATOM     31  C   UNL     1      44.796  42.682 -13.276  1.00  0.00     0.410 C 
ATOM     32  O   UNL     1      45.993  43.002 -13.534  1.00  0.00    -0.224 OA
BRANCH  31  33
ATOM     33  O   UNL     1      44.258  41.712 -14.117  1.00  0.00    -0.444 OA
BRANCH  33  34
ATOM     34  C   UNL     1      44.590  41.670 -15.512  1.00  0.00     0.106 C 
ATOM     35  C   UNL     1      45.968  41.143 -15.750  1.00  0.00     0.056 C 
ATOM     36  C   UNL     1      43.552  40.736 -16.148  1.00  0.00     0.056 C 
ATOM     37  C   UNL     1      44.374  43.034 -16.133  1.00  0.00     0.056 C 
ENDBRANCH  33  34
ENDBRANCH  31  33
ENDBRANCH  28  31
ENDBRANCH   4  26
TORSDOF 6
ENDMDL
"""


def test_pdbqt_meeko_macrocycle():
    molecule = Chem.MolFromSmiles(
        "Cn1nc(C2CCCN2C(=O)OC(C)(C)C)cc1C(=O)NC1CCc2ccccc2NC1=O"
    )
    result = pdbqt2molblock(OFFENDING_PDBQT, molecule, "Test-molecule")

    assert result is not None
