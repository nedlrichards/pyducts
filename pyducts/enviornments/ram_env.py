import os
import numpy as np
from scipy.interpolate import interp1d

class RAMEnv:
    """Basic enviornment for RAM run"""
    def __init__(self):
        """setup common parameters for all ranges"""
        self.z_src = None
        self.r_axis = None
        self.z_axis = None
        self.ndr = None
        self.ndz = None
        self.facous = None
        self.ram_info = {"c0":None, "num_pade":None, "stability":None}
        self.bathy = None
        self.profile_range = None
        self.profiles = None
        self.eta = 1.0 / (40.0 * np.pi * np.log10(np.e))

    def setup(self, z_src=None, r_axis=None, z_axis=None, facous=None,
              ram_info=None, bathy=None, profile_range=None, profiles=None):
        """selectively update information, overwrite all input information"""
        if z_src is not None: self.z_src = z_src
        if r_axis is not None: self.r_axis = r_axis
        if z_axis is not None: self.z_axis = z_axis
        if facous is not None: self.facous = facous
        if ram_info is not None: self.ram_info = ram_info
        if bathy is not None: self.bathy = bathy
        if profile_range is not None: self.profile_range = profile_range
        if profiles is not None: self.profiles = profiles

    def write_in(self, file_name='ram.in', title=None, z_rcr=-50,
                 write_profiles=True):
        """write current enviornment to an input file"""
        # find bathymetry index related to z_axis
        bz_i = list(self.profiles.shape).index(z_axis.size)
        bp_i = 1 if bz_i is 0 else 0
        num_p = profiles.shape[bp_i]

        # check if bathymetry is repeated, or same number of changes as profiles
        if len(bathy.shape) == 1:
            is_same_bathy = True
        else:
            is_same_bathy = False

        with open(file_name, 'w') as f:
            if title is not None:
                f.write(title + '\n')
            else:
                f.write('autogen' + '\n')

            f.write("{:.2f}  {:.2f}  {:.2f}".format(self.facous, self.z_src, z_rcr))
            f.write('\n')
            dr = (self.r_axis[-1] - self.r_axis[0]) / (self.r_axis.size - 1)
            f.write("{:.2f}  {:.2f}  {:d}".format(self.r_axis[-1], dr, self.ndr))
            f.write('\n')
            dz = (self.z_axis[-1] - self.z_axis[0]) / (self.z_axis.size - 1)
            f.write("{:.2f}  {:.2f}  {:d}  {:d}".format(self.z_axis[-1], dz, self.ndz, self.ndz))
            f.write('\n')
            f.write("{:.2f}  {:.2f}  {:d}  {:.2f}".format(self.ram_info["c0"],
                                                          self.ram_info["num_pade"],
                                                          self.ram_info["stability"][0],
                                                          self.ram_info["stability"][1]))
            # TODO: remove hardcode of bathy layer
            f.write("{:.2f}  {:.2f}".format(0.0, 5199))
            f.write("{:d}  {:d}".format(-1, -1))
            # write profiles to file
            #for i in num_p:
            i = 0
            self._append_profile(i, bz_i, is_same_bathy, f)
            f.write("{:d}  {:d}".format(-1, -1))
            self._append_bathy(i, bz_i, is_same_bathy, f)


    def _append_profile(self, ind, bz_i, is_same_bathy, fid):
        """Append a profile to file"""
        bp_i = (Ellipsis, ind) if bz_i is 1 else (ind, Ellipsis)
        printer = np.array([self.profiles[bp_i], self.z_axis]).T
        np.tofile(fid, printer)

    def _append_bathy(self, ind, bz_i, is_same_bathy, fid):
        # TODO: remove hardcode
        f.write("{:.2f}  {:.2f}".format(0.0, 1544.9))
        fid.write("{:d}  {:d}".format(-1, -1))
        f.write("{:.2f}  {:.2f}".format(0.0, 1.0))
        fid.write("{:d}  {:d}".format(-1, -1))
        f.write("{:.2f}  {:.2f}".format(0.0, 0.05))
        fid.write("{:d}  {:d}".format(-1, -1))

    def read_in(self, file_name, read_profiles=True):
        """read RAM imput from file_name"""
        with open(file_name, 'r') as f:
            title = f.readline()  # title
            [facous, z_src, z_rcr] = np.array(f.readline().split()[:3]).astype(float)
            self.facous = facous
            self.z_src = z_src
            [rmax, dr, ndr] = np.array(f.readline().split()[:3]).astype(float)
            self.ndr = int(ndr)
            [zmax, dz, ndz, zmplt] = np.array(f.readline().split()[:4]).astype(float)
            self.ndz = int(ndz)
            [c0, num_pade, num_stab, rs] = np.array(f.readline().split()[:4]).astype(float)
            self.ram_info["c0"] = c0
            self.ram_info["num_pade"] = int(num_pade)
            self.ram_info["stability"] = (int(num_stab), rs)

            self.z_axis = np.arange(np.ceil(zmax / dz) + 1) * dz
            self.r_axis = (np.arange(np.ceil(rmax / dr)) + 1) * dr

            if not read_profiles:
                return

            # read bottom position
            bathy = [np.array(f.readline().split()).astype(float)]
            while bathy[-1][0] >= 0:
                bathy.append(np.array(f.readline().split()).astype(float))
            bathy[-1] = np.array([2 * rmax, bathy[-2][-1]])
            bathy = np.array(bathy)
            self.bathy = bathy

            # read profiles one at a time
            profiles = []
            rp = 0.0
            while isinstance(rp, float):
                cw = [np.array(f.readline().split()).astype(float)]
                while cw[-1][0] >= 0:
                    cw.append(np.array(f.readline().split()).astype(float))
                cw = np.array(cw[:-1])

                cb = [np.array(f.readline().split()).astype(float)]
                while cb[-1][0] >= 0:
                    cb.append(np.array(f.readline().split()).astype(float))
                cb = np.array(cb[:-1])

                rhob = [np.array(f.readline().split()).astype(float)]
                while rhob[-1][0] >= 0:
                    rhob.append(np.array(f.readline().split()).astype(float))
                rhob = np.array(rhob[:-1])

                attn = [np.array(f.readline().split()).astype(float)]
                while attn[-1][0] >= 0:
                    attn.append(np.array(f.readline().split()).astype(float))
                attn = np.array(attn[:-1])

                profiles.append({"rp":rp,
                                "cw":cw,
                                "cb":cb,
                                "rhob":rhob,
                                "attn":attn})

                rp = f.readline().split()
                if rp:
                    rp = float(rp[0])
                else:
                    break
            self.profiles = profiles

    def populate_quantities(self, profile_num):
        """interpolate input parameters to grid used by ram"""
        omega = 2 * np.pi * self.facous
        k0 = omega / self.ram_info["c0"]

        cw = self.profiles[profile_num]['cw']
        cb = self.profiles[profile_num]['cb']
        rhob = self.profiles[profile_num]['rhob']
        attn = self.profiles[profile_num]['attn']

        # populate verical profiles
        if cw.shape[0] == 1:
            cw_up = np.full_like(self.zaxis, cw[0, 1])
        else:
            ier = interp1d(cw[:, 0], cw[:, 1], kind=1, fill_value="extrapolate")
            cw_up = ier(self.zaxis)

        if cb.shape[0] == 1:
            cb_up = np.full_like(self.zaxis, cb[0, 1])
        else:
            ier = interp1d(cb[:, 0], cb[:, 1], kind=1, fill_value="extrapolate")
            cb_up = ier(self.zaxis)

        if rhob.shape[0] == 1:
            rhob_up = np.full_like(self.zaxis, rhob[0, 1])
        else:
            ier = interp1d(rhob[:, 0], rhob[:, 1], kind=1, fill_value="extrapolate")
            rhob_up = ier(self.zaxis)

        if attn[0,0] >= 0:
            attn_up = np.zeros_like(self.zaxis)
            # assume zero attenuation at z=0

        if attn.shape[0] == 1:
            attn_up = np.full_like(self.zaxis, attn[0, 1])

        else:
            ier = interp1d(attn[:, 0], attn[:, 1], kind=1, fill_value="extrapolate")
            zi = self.zaxis > attn[0,0]
            attn_up[zi] = ier(self.zaxis[zi])

        # compute values used in ram calculations
        ksqw = (omega / cw_up) ** 2 - k0 ** 2
        ksqb = ((omega / cb_up) *(1.0 + 1.0j * self.eta * attn_up)) ** 2 - k0 ** 2
        alpw = np.sqrt(cw_up / self.ram_info["c0"])
        alpb = np.sqrt(rhob_up * cb_up / self.ram_info["c0"])

        return cw_up, cb_up, rhob_up, attn_up, ksqw, ksqb, alpw, alpb
