
Thanks to Uwe, test failure on M1mac actually due to eigen decomposition sign identifiability. A change so that only the non-decomposed matrix is tested for equality should make it pass.

