import numpy as np
import shutil
from subprocess import check_output
from ram_wrapper import read_line, read_grid

ram_exe = '/home/nedrichards/bin/ram.exe'
ram_typed_exe = '/home/nedrichards/bin/ram_typed.exe'

# run example from RAM documentation
shutil.copy('./tests/readme.in', 'ram.in')
rout = check_output(ram_exe)

# load transmission loss from tl.line
r, tl = read_line("tl.line")
rg, zg, tlg = read_grid("tl.grid", num_bytes=4)

# compare with typed ram output
rout = check_output(ram_typed_exe)
r_t, tl_t = read_line("tl.line")
rg_t, zg_t, tlg_t = read_grid("tl.grid", num_bytes=8)

# typing changes tl by a bit, put absolute threshold on accuracy
eps = 1e-2
assert(eps > np.max(np.abs(tl-tl_t)))
assert(eps > np.max(np.abs(tlg-tlg_t)))


