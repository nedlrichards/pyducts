import numpy as np
from math import pi, e

eta = 1 / (40 * pi * np.log10(e)) # attenuation unit conversion

class RamIn:
    """Read in from ram.in"""

    def __init__(self, file_in="ram.in"):
        """Read inputs from ram.in"""

        def _bathy(rmax, fid):
            """Read bathymetric data line by line"""
            rb = 0.
            while rb >= 0.:
                oput = fid.readline()
                rb, zb = np.array(oput.split()[:2], dtype=np.float64)
                if rb >= 0.:
                    last_zb = zb
                    yield rb, zb

            # extend constant bathymetry past end of computation range
            yield 2 * rmax, last_zb

        # read file
        with open(file_in, "r") as fid:
            title = fid.readline()

            # read tun specifications from file
            specs = []
            specs += fid.readline().split()[:3]
            specs += fid.readline().split()[:3]
            specs += fid.readline().split()[:4]
            specs += fid.readline().split()[:4]

            self.bathy = np.array(list(_bathy(self.rmax, fid)))

            all_profiles = list(self.read_profiles(fid))

            cw         = np.array([p[0] for p in all_profiles])
            cb         = np.array([p[1] for p in all_profiles])
            self.rho_b = np.array([p[2] for p in all_profiles])
            attn       = np.array([p[3] for p in all_profiles])

            # 0 range is added to ranges of profiles
            rp         = np.array([0.] + [p[4] for p in all_profiles])

            # assign specifications to variables
            specs = specs[::-1]  # flip to make pop more intuative
            self.freq = float(specs.pop())
            self.z_src = float(specs.pop())
            self.z_rcr = float(specs.pop())
            self.rmax = float(specs.pop())
            self.dr = float(specs.pop())
            self.range_decimation = int(specs.pop())
            self.zmax = float(specs.pop())
            self.dz = float(specs.pop())
            self.z_decimation = int(specs.pop())
            self.zmax_plot = float(specs.pop())
            self.c0 = float(specs.pop())
            self.num_pade = int(specs.pop())
            self.num_stability = int(specs.pop())
            self.range_stability = float(specs.pop())

            self.omega = 2.0 * pi * self.freq
            self.k0 = self.omega / self.c0

            # reciver specifications
            receiver_zbins = 1.0 + self.z_rcr / self.dz
            self.receiver_zi = int(receiver_zbins)
            self.receiver_zmod = receiver_zbins - self.receiver_zi

            num_z = int(self.zmax / self.dz - 0.5)
            self.zaxis = (np.arange(num_z) + 1) * self.dz
            self.z_save = self.zaxis[self.zaxis < self.zmax_plot]
            self.z_save = self.z_save[::self.z_decimation]

            # find closest index of the bottom
            z_bottom = self.bathy[0, 1]
            bottom_index = int(1.0 + z_bottom / self.dz)
            bottom_index = max(2, bottom_index)
            self.bottom_index = min(num_z, bottom_index)

            if self.range_stability < self.dr:
                self.range_stability = 2.0 * self.rmax

            # computed quantities from interpolated profiles
            self.ksqw = (self.omega / cw) ** 2 - self.k0 ** 2
            self.ksqb = ((self.omega / cb) * (1.0 + 1j * eta * attn)) ** 2 \
                      - self.k0 ** 2
            self.alpw = np.sqrt(cw / self.c0)
            self.alpb = np.sqrt(self.rho_b * cb / self.c0)

    def read_profiles(self, fid):
        """Read profile blocks of ram.in in sequence"""

        def _profile(zmax, fid):
            """Read single property profile in line by line"""
            zprof = 0.

            while zprof >= 0.:
                oput = fid.readline().split()[:2]
                zprof, value = np.array(oput, dtype=np.float64)
                if zprof >= 0.:
                    last_value = value
                    yield zprof, value

            yield zmax, last_value

        # read profiles at each defined range
        rp = 0.
        while rp < self.rmax:
            cw_points = np.array(list(_profile(self.zmax, fid)))
            cb_points = np.array(list(_profile(self.zmax, fid)))
            rho_b_points = np.array(list(_profile(self.zmax, fid)))
            attn_points = np.array(list(_profile(self.zmax, fid)))

            # end of file will occur here
            next_l = fid.readline()
            if next_l:
                rp = float(next_l.split()[0])
            else:
                # break condition
                rp = 2 * self.rmax

            # interpolate profiles to ram grid
            cw = np.interp(self.zaxis, cw_points[:, 0], cw_points[:, 1])
            cb = np.interp(self.zaxis, cb_points[:, 0], cb_points[:, 1])
            rho_b = np.interp(self.zaxis, rho_b_points[:, 0], rho_b_points[:, 1])
            attn = np.interp(self.zaxis, attn_points[:, 0], attn_points[:, 1])

            yield cw, cb, rho_b, attn, rp
