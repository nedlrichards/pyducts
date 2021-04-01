import os
import numpy as np
from scipy.interpolate import interp1d

class Enviornment:
    """Basic enviornment for RAM run"""
    def __init__(self):
        """setup common parameters for all ranges"""
        self.z_src = None
        self.r_axis = None
        self.z_axis = None
        self.facous = None
        self.ram_info = {"c0":None, "num_pade":None, "stability":None}
        self.bathy = None
        self.profiles = None
        self.eta = 1.0 / (40.0 * np.pi * np.log10(np.e))

    def read_in(self, file_name):
        """read RAM imput from file_name"""
        with open(file_name, 'r') as f:
            title = f.readline()  # title
            [facous, z_src, z_rcr] = np.array(f.readline().split()).astype(float)
            self.facous = facous
            self.z_src = z_src
            [rmax, dr, ndr] = np.array(f.readline().split()).astype(float)
            ndr = int(ndr)
            [zmax, dz, ndz, zmplt] = np.array(f.readline().split()).astype(float)
            ndz = int(ndz)
            [c0, num_pade, num_stab, rs] = np.array(f.readline().split()).astype(float)
            self.ram_info["c0"] = c0
            self.ram_info["num_pade"] = int(num_pade)
            self.ram_info["stability"] = (int(num_stab), rs)

            self.zaxis = np.arange(np.floor(zmax / dz) + 1) * dz
            self.raxis = (np.arange(np.floor(rmax / dr)) + 1) * dr

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
        eta = 1.0 / (40.0 * np.pi * np.log10(np.e))

        ksqw = (omega / cw_up) ** 2 - k0 ** 2
        ksqb = ((omega / cb_up) *(1.0 + 1.0j * self.eta * attn_up)) ** 2 - k0 ** 2
        alpw = np.sqrt(cw_up / self.ram_info["c0"])
        alpb = np.sqrt(rhob_up * cb_up / self.ram_info["c0"])

        return cw_up, cb_up, rhob_up, attn_up, ksqw, ksqb, alpw, alpb


class Layer:
    """A medium with smoothly varing properties"""
    def __init__(self, rzw_l):
        """initialize layer with continous material properties"""
        # 3xN dimensional array, range, z origin, and width of layer
        self.rzw_l = rzw_l
        self.p_range = []
        self.properties = []

    def add_profile(self, range_profile, cp_l, cs_l=None, rho_l=1e3,
                    ap_l=0., as_l=None):
        """Add properties to layer, specified as a vertical profile"""
        # Following properties may be scalar, or 2xN specified as:
        # depth from layer origin (m)
        # property value
        #
        # cp_l: compressional sound speed in layer (m/s)
        # rho_l: density in layer (kg / m^3)
        # ap_l: compressional attenuation
        # cs_l: shear sound speed in layer (m/s)
        # as_l: shear attenuation

        self.p_range.append(range_profile)
        self.properties.append([cp_l, rho_l, ap_l, cs_l, as_l])

    def populate_grid(self, acoustic_property, raxis, zaxis, grided_data=None):
        """interpolate an acoustic property to a (r, z) grid"""


class Profile:
    """Define properties at a single range value"""
    def __init__(self, r_profile):
        """define location of profile"""
        self.r_profile = r_profile
        self.layers = []

    def add_layer(self, width_layer, z_layer, cp_layer,
                  cs_layer=None, rho_layer=None, attn_layer=None):
        """define a layer"""
        new_l = Layer(layer_width, z_layer, cp_layer, cs_layer,
                      rho_layer, attn_layer)
        self.layers.append(new_l)

    def add_halspace(self, width_hs, z_hs, cp_hs,
                     cs_hs=None, rho_hs=None, attn_hs=None):
        """defined a half-space, or absorbing layer at grid edge"""
        self.width_hs = width_hs  # height of absorbing layer
        self.z_hs = z_hs  # hs depth corrdinate
        self.cp_hs = cp_hs
        self.cs_hs = cs_hs
        self.rho_hs = rho_hs
        self.attn_hs = attn_hs


