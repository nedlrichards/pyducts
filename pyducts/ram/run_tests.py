import numpy as np
import shutil
from subprocess import check_output
from ram_wrapper import read_line, read_grid

ram_exe = '/home/nedrichards/bin/ram.exe'
ram_typed_exe = '/home/nedrichards/bin/ram_typed.exe'

# typing changes tl by a bit, put absolute threshold on accuracy
eps = 0.1

# run example from RAM documentation
shutil.copy('./tests/readme.in', 'ram.in')
rout = check_output(ram_exe)

# load transmission loss from tl.line
r, tl = read_line("tl.line", is_cmpx=False)
rg, zg, tlg = read_grid("tl.grid", num_bytes=4, is_cmpx=False)

# compare with typed ram output
rout = check_output(ram_typed_exe)
r_t, tl_t = read_line("tl.line")
rg_t, zg_t, tlg_t = read_grid("tl.grid", num_bytes=8)

# confirm that all error occurs at very low tl
if np.any(np.abs(tl-tl_t) > eps):
    assert(np.min(tl_t[np.abs(tl-tl_t) > eps]) > 80)
if np.any(np.abs(tlg-tlg_t) > eps):
    assert(np.min(tlg_t[np.abs(tlg-tlg_t) > eps]) > 80)

# run example from RAM documentation
shutil.copy('./tests/spice.in', 'ram.in')
rout = check_output(ram_exe)

# load transmission loss from tl.line
r, tl = read_line("tl.line", is_cmpx=False)
rg, zg, tlg = read_grid("tl.grid", num_bytes=4, is_cmpx=False)

# compare with typed ram output
rout = check_output(ram_typed_exe)
r_t, tl_t = read_line("tl.line")
rg_t, zg_t, tlg_t = read_grid("tl.grid", num_bytes=8, is_cmpx=True)

# confirm that all error occurs at very low tl
if np.any(np.abs(tl-tl_t) > eps):
    assert(np.min(tl_t[np.abs(tl-tl_t) > eps]) > 80)
if np.any(np.abs(tlg-tlg_t) > eps):
    assert(np.min(tlg_t[np.abs(tlg-tlg_t) > eps]) > 80)


