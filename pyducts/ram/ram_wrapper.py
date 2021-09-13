import numpy as np
from os import path
import subprocess
import struct

ram_path = path.abspath(path.expanduser('~/bin'))

def run_ram():
    """run kraken on an env file"""
    run_cmd = path.join(ram_path, 'ram_typed.exe')
    #run_cmd = path.join(ram_path, 'ram.exe')
    # TODO: see if status gives anything
    status = subprocess.run(run_cmd, check=True)

def read_grid(is_standard=False):
    """read tl.grid"""
    with open("tl.grid", "rb") as fid:

        if is_standard:
            num_bytes = 4
            int_type = "i"
            float_type = "f"
        else:
            num_bytes = 8
            int_type = "q"
            float_type = "d"

        # read header
        # fortran adds a field at every write
        temp = struct.unpack("i", fid.read(4))

        freq = struct.unpack(float_type, fid.read(num_bytes))[0]
        z_src = struct.unpack(float_type, fid.read(num_bytes))[0]
        z_rcr = struct.unpack(float_type, fid.read(num_bytes))[0]
        r_max = struct.unpack(float_type, fid.read(num_bytes))[0]
        dr = struct.unpack(float_type, fid.read(num_bytes))[0]
        ndr = struct.unpack(int_type, fid.read(num_bytes))[0]
        z_max = struct.unpack(float_type, fid.read(num_bytes))[0]
        dz = struct.unpack(float_type, fid.read(num_bytes))[0]
        ndz = struct.unpack(int_type, fid.read(num_bytes))[0]
        zmax_plot = struct.unpack(float_type, fid.read(num_bytes))[0]
        c0 = struct.unpack(float_type, fid.read(num_bytes))[0]
        num_p = struct.unpack(int_type, fid.read(num_bytes))[0]
        num_s = struct.unpack(int_type, fid.read(num_bytes))[0]
        range_s_max = struct.unpack(float_type, fid.read(num_bytes))[0]
        num_zplot = struct.unpack(int_type, fid.read(num_bytes))[0]

        # fortran adds a field at every write
        temp = struct.unpack("i", fid.read(4))

        num_r_read = int((r_max / dr) / ndr)

        zplot = np.arange(num_zplot) * ndz * dz
        rplot = (np.arange(num_r_read) + 1) * ndr * dr

        p_out = []

        for i in np.arange(num_r_read):
            temp = struct.unpack("i", fid.read(4))
            if is_standard:
                TL = np.fromfile(fid, count=num_zplot, dtype=np.float32)
                p_out.append(TL)
            else:
                p = np.fromfile(fid, count=2 * num_zplot, dtype=np.float64)
                p_out.append(p[::2] + 1j * p[1::2])
            temp = struct.unpack("i", fid.read(4))
        p_out = np.array(p_out)

    return zplot, rplot, p_out

def read_line(is_standard=False):
    """read tl.line output by ram"""
    line = np.loadtxt("tl.line")
    if is_standard:
        p_line = line[:, 1]
    else:
        p_line = line[:, 1] + 1j * line[:, 2]
    r_line = line[:, 0]
    return r_line, p_line

class RamIn:
    """basic capability ram.in writting"""

    def __init__(self, fc, z_src, r_max, z_bottom, bottom_HS=[1600., 1800., 0.5],
                 z_rcr=None, dr=200., zmax_plot=None, ndz=None):
        """setup front matter in ram.in"""
        self.fc = fc
        self.z_src = z_src
        self.r_max = r_max
        self.z_bottom = z_bottom
        self.bottom_HS = bottom_HS
        self.dr = dr

        # ram parameters are kept generic
        self.c0 = 1500.
        self.num_p = 6
        self.num_stab = 1
        self.stab_range_max = 0.

        # ssp
        self.sigma = 0

        # steps and decimation
        self.ndr = 1
        dz_decimation = 15

        num_points = int(dz_decimation * np.ceil(self.z_bottom * fc / self.c0))
        self.dz = self.z_bottom / (num_points - 1)

        if ndz is None:
            # aim for about lambda / 2 spacing on save grid
            self.ndz = dz_decimation // 2
        else:
            self.ndz = ndz

        # ram seems to clip the bottom of plot grid
        if zmax_plot is None:
            self.zmax_plot = self.z_bottom + 2 * self.dz
        else:
            self.zmax_plot = zmax_plot + 2 * self.dz

        # check if TL line depth is specified
        if z_rcr is None:
            self.z_rcr = z_src
        else:
            self.z_rcr = z_rcr

        # bottom is assumed flat
        rb = 0
        self.z_hs = self.z_bottom + 300

        # add a sponge layer
        self.z_sponge = self.z_hs + 300
        self.attn_sponge = [0.5, 10]


    def write_frontmatter(self):
        """write frontmatter common to all profiles"""

        with open("ram.in", 'w') as writer:
            writer.write("auto_gen\n")
            writer.write(f'{self.fc:.4f} {self.z_src:.4f} {self.z_rcr:.4f}\n')
            writer.write(f'{self.r_max:.4f} {self.dr:.4f} {self.ndr}\n')
            writer.write(f'{self.z_sponge:.4f} {self.dz:.8f} {self.ndz} {self.zmax_plot:.4f}\n')
            writer.write(f'{self.c0:.4f} {self.num_p} {self.num_stab} {self.stab_range_max:.4f}\n')

            # bathy depth
            writer.write(f"0.0 {self.z_bottom:.4f}\n" )
            writer.write("-1 -1\n" )


    def write_profile(self, rprof, z, ssp):
        """Add profile to ram file"""
        with open("ram.in", 'a') as writer:
            # need to write profile range except when r == 0
            if rprof > 1e-3:
                writer.write(f"{rprof:.4f}\n" )

            # ssp
            np.savetxt(writer, np.array([z, ssp]).T, fmt='%.4f', newline='\n', delimiter=' ')
            writer.write("-1 -1\n" )

            # bathy sound speed
            writer.write(f"0.0 {self.bottom_HS[0]:.4f}\n" )
            writer.write("-1 -1\n" )

            # bathy density
            writer.write(f"0.0 {self.bottom_HS[1]*1e-3:.4f}\n" )
            writer.write("-1 -1\n" )

            # bathy attn
            if len(self.bottom_HS) > 2:
                writer.write(f"0.0 {self.bottom_HS[2]:.4f}\n" )
                writer.write(f"{self.z_hs - 0.1:.4f} {self.bottom_HS[2]:.4f}\n")

            writer.write(f"{self.z_hs:.4f} {self.attn_sponge[0]:.4f}\n" )
            writer.write(f"{self.z_sponge:.4f} {self.attn_sponge[1]:.4f}\n" )
            writer.write("-1 -1\n" )
