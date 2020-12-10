import numpy as np
import os

class Enviornment:
    """Basic enviornment for RAM run"""
    def __init__(self):
        """setup common parameters for all ranges"""
        self.z_src = None
        self.r_axis = None
        self.z_axis = None
        self.bottom_position = None

    def read_ram_in(self, file_name):
        """Read ram.in"""

        with open(file_name, 'r') as f:
            _ = f.readline()
            [facous, z_src, z_rcr] = f.readline().split()
            self.z_src = z_src
            [rmax, dr, ndr] = f.readline().split()
            self.r_axis = np.arange(np.ceil(rmax / dr)) * dr
            [zmax, dz, ndz, zmplt] = f.readline().split()
            self.z_axis = np.arange(np.ceil(zmax / dz)) * dz
            [c0, npr, ns, rs] = f.readline().split()

            # read bottom position
            bathy = [np.array(f.readline().split()).astype(float)]
            while bathy[-1][0] >= 0:
                bathy.append(np.array(f.readline().split()).astype(float))
            self.bottom_position = np.array(bathy[:-1])

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
                rp = f.readline().split()
            return cw, cb, rhob, attn

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

class Layer:
    """A medium with smoothly varing properties"""
    def __init__(self, width_l, z_l, cp_l, cs_l, rho_l, attn_l):
        """initialize layer with material properties"""
        self.width_l = width_l
        self.z_l = z_l
        self.cp_l = cp_l
        self.cs_l = cs_l
        self.rho_l = rho_l
        self.atten_l = atten_l
