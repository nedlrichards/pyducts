import numpy as np
import shutil
from subprocess import check_output
from ram_wrapper import read_line, read_grid

ram_exe = './bin/ram'
ram_typed_exe = './bin/ram_v1'

# typing changes tl by a bit, put absolute threshold on accuracy
eps = 0.1

# run example from RAM documentation
shutil.copy('./tests/readme.in', 'ram.in')
rout = check_output(ram_exe)

# load transmission loss from tl.line
r, tl = read_line(is_standard=True)
rg, zg, tlg = read_grid(is_standard=True)

# compare with typed ram output
rout = check_output(ram_typed_exe)
r_t, p_t = read_line(is_standard=False)
tl_t = -20 * np.log10(np.abs(p_t) + np.spacing(1))
rg_t, zg_t, p_t = read_grid(is_standard=False)
tlg_t = -20 * np.log10(np.abs(p_t) + np.spacing(1))

# confirm that all error occurs at very low tl
if np.any(np.abs(tl-tl_t) > eps):
    assert(np.min(tl_t[np.abs(tl-tl_t) > eps]) > 80)
if np.any(np.abs(tlg-tlg_t) > eps):
    assert(np.min(tlg_t[np.abs(tlg-tlg_t) > eps]) > 80)

# run spice example
shutil.copy('./tests/spice.in', 'ram.in')
rout = check_output(ram_exe)

# load transmission loss from tl.line
r, tl = read_line(is_standard=True)
rg, zg, tlg = read_grid(is_standard=True)

# compare with typed ram output
rout = check_output(ram_typed_exe)
r_t, p_t = read_line(is_standard=False)
tl_t = -20 * np.log10(np.abs(p_t) + np.spacing(1))
rg_t, zg_t, p_t = read_grid(is_standard=False)
tlg_t = -20 * np.log10(np.abs(p_t) + np.spacing(1))

# confirm that all error occurs at very low tl
if np.any(np.abs(tl-tl_t) > eps):
    assert(np.min(tl_t[np.abs(tl-tl_t) > eps]) > 80)
if np.any(np.abs(tlg-tlg_t) > eps):
    assert(np.min(tlg_t[np.abs(tlg-tlg_t) > eps]) > 80)
