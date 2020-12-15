import numpy as np
import matplotlib.pyplot as plt
import shutil
from subprocess import check_output

plt.ion()

ram_exe = '/home/nedrichards/bin/ram.exe'

# run example from RAM documentation
shutil.copy('readme.in', 'ram.in')
rout = check_output(ram_exe)

# load transmission loss from tl.line
test_file = "tl.line"

tl = np.loadtxt(test_file)

fig, ax = plt.subplots()
ax.plot(tl[:, 0], tl[:, 1])
ax.set_ylim(100, 50, 'o')
