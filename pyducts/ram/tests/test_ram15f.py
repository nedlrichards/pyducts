import numpy as np
import matplotlib.pyplot as plt

plt.ion()

# load transmission loss from tl.line
test_file = "tl.line"

tl = np.loadtxt(test_file)

fig, ax = plt.subplots()
ax.plot(tl[:, 0], tl[:, 1])
ax.set_ylim(100, 50, 'o')
