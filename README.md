# MLatomization
Development of a Python code to predict atomization energies for molecular systems

dsgdb7ae2.xyz contains the coordinates of about 7000 molecules and the corresponding atomization energies 
the coordinates are the only quantity used to build the features
900 molecules are used to train (chosen randomly), 100 molecules are used as hold out test set, and the remaining are used to predict. TODO: Cross validation could be more reliable

test_molecule_laplacian.py: contains the best working script (mean absolute error of ~8 kcal/mol) that uses the kernel ridge regression algorithm with a Laplacian kernel

test_molecule.py uses kernel ridge regression algorithm with a Gaussian kernel

grid_molecule.py search on a grid for the best possible parameters

test_molecule_scikit.py uses scikit-learn libraries to test different algorithms. SVM gives the same results as kernel ridge rgression while neural networks give worst results 
