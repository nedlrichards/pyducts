import numpy as np
from os import path
import subprocess
import struct

ram_path = path.abspath(path.expanduser('~/bin'))

def run_ram():
    """run kraken on an env file"""
    run_cmd = path.join(ram_path, 'ram_typed.exe')
    # TODO: see if status gives anything
    status = subprocess.run(run_cmd, check=True)
    return

def write_ram_in(fc, z_src, r_max, z, ssp, bottom_HS=[1600., 1.8], z_rcr=None):
    """basic capability ram.in writting"""

    # front matter
    title = "auto_gen"

    # ssp
    z = np.abs(z)
    sigma = 0
    zb = z[-1]

    # steps and decimation
    dr = 200.
    ndr = 1
    dz = 2.
    ndz = 1
    zmplt = zb

    # check if TL line depth is specified
    if z_rcr is None:
        z_rcr = z_src

    # bottom is assumed flat
    rb = 0
    z_hs = zb + 300

    # add a sponge layer
    z_sponge = z_hs + 300
    attn_sponge = [0.5, 10]

    # ram parameters are kept generic
    c0 = 1500.
    num_p = 6
    num_stab = 1
    stab_range_max = 0.

    with open("ram.in", 'w') as writer:
        writer.write(title + '\n')
        writer.write(f'{fc:.2f} {z_src:.2f} {z_rcr:.2f}\n')
        writer.write(f'{r_max * 1e3:.2f} {dr:.2f} {ndr}\n')
        writer.write(f'{z_sponge:.2f} {dz:.2f} {ndz} {zmplt:.2f}\n')
        writer.write(f'{c0:.2f} {num_p} {num_stab} {stab_range_max:.2f}\n')

        # bathy depth
        writer.write(f"0.0 {zb:.2f}\n" )
        writer.write("-1 -1\n" )

        # ssp
        np.savetxt(writer, np.array([z, ssp]).T, fmt='%.2f', newline='\n', delimiter=' ')
        writer.write("-1 -1\n" )

        # bathy sound speed
        writer.write(f"0.0 {bottom_HS[0]:.2f}\n" )
        writer.write("-1 -1\n" )

        # bathy density
        writer.write(f"0.0 {bottom_HS[1]:.2f}\n" )
        writer.write("-1 -1\n" )

        # bathy attn
        writer.write(f"{z_hs:.2f} {attn_sponge[0]:.2f}\n" )
        writer.write(f"{z_sponge:.2f} {attn_sponge[1]:.2f}\n" )
        writer.write("-1 -1\n" )

def read_grid():
    """read tl.grid"""
    with open("tl.grid", "rb") as fid:
        # read header
        # fortran adds a field at every write
        temp = struct.unpack("i", fid.read(4))
        freq = struct.unpack("d", fid.read(8))[0]
        z_src = struct.unpack("d", fid.read(8))[0]
        z_rcr = struct.unpack("d", fid.read(8))[0]
        r_max = struct.unpack("d", fid.read(8))[0]
        dr = struct.unpack("d", fid.read(8))[0]
        ndr = struct.unpack("q", fid.read(8))[0]
        z_max = struct.unpack("d", fid.read(8))[0]
        dz = struct.unpack("d", fid.read(8))[0]
        ndz = struct.unpack("q", fid.read(8))[0]
        zmplt = struct.unpack("d", fid.read(8))[0]
        c0 = struct.unpack("d", fid.read(8))[0]
        num_p = struct.unpack("q", fid.read(8))[0]
        num_s = struct.unpack("q", fid.read(8))[0]
        range_s_max = struct.unpack("d", fid.read(8))[0]
        num_zplot = struct.unpack("q", fid.read(8))[0]
        # fortran adds a field at every write
        temp = struct.unpack("i", fid.read(4))

        num_r_read = int((r_max / dr) / ndr)

        zplot = np.arange(num_zplot) * zmplt / (num_zplot-1)
        rplot = (np.arange(num_r_read) + 1) * r_max / num_r_read

        cmpx_TL = np.zeros((num_zplot, num_r_read)) \
                + 1j * np.zeros((num_zplot, num_r_read))

        for i in np.arange(num_r_read):
            temp = struct.unpack("i", fid.read(4))
            TL = np.fromfile(fid, count=2 * num_zplot)
            cmpx_TL[:, i] = TL[::2] + 1j * TL[1::2]
            temp = struct.unpack("i", fid.read(4))
    return zplot, rplot, cmpx_TL
