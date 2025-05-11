import numpy as np

a = [1, 2, -7, -8, 12]
roots = np.roots(a)
for i in range(len(roots)):
    print(roots[i])