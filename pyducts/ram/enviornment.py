import numpy as np
import os

class Enviornment:
    """Basic enviornment for RAM run"""
    def __init__(self):
        """setup common parameters for all ranges"""
        self.z_src = None
        self.r_axis = None
        self.z_axis = None
        self.layers = []

    def read_ram_in(self, file_name):
        """Read ram.in"""

        with open(file_name, 'r') as f:
            _ = f.readline()
            [facous, z_src, z_rcr] = np.array(f.readline().split()).astype(float)
            self.z_src = z_src
            [rmax, dr, ndr] = np.array(f.readline().split()).astype(float)
            self.r_axis = np.arange(np.ceil(rmax / dr)) * dr
            [zmax, dz, ndz, zmplt] = np.array(f.readline().split()).astype(float)
            self.z_axis = np.arange(np.ceil(zmax / dz)) * dz
            [c0, npr, ns, rs] = np.array(f.readline().split()).astype(float)

            # read bottom position
            bathy = [np.array(f.readline().split()).astype(float)]
            while bathy[-1][0] >= 0:
                bathy.append(np.array(f.readline().split()).astype(float))
            bottom_position = np.array(bathy[:-1])

            # ram assumes two layers, one water and one bottom
            rzw = np.concatenate([bottom_position[:, 0][:, None],
                                  np.zeros((2, 1)),
                                  bottom_position[:, 1][:, None]],
                                  axis=1)
            water_layer = Layer(rzw)

            # have bottom layer cover entire grid, rely on layer order to
            # insert water layer where it is defined
            rzw = np.array([0., 0., zmax])
            bottom_layer = Layer(rzw)

            self.layers.append(water_layer)
            self.layers.append(bottom_layer)

            # read profiles one at a time
            profiles = []
            rp = 0.0
            while rp < rmax:
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
                # attn is assumed zero if not specified at z=0
                if attn[0, 0] > 0:
                    attn = np.concatenate([np.zeros((1, 2)), attn])

                self.layers[0].add_profile(rp, cw)
                self.layers[1].add_profile(rp, cb, rho_l=rhob, ap_l=attn)

                rp = f.readline().split()
                if rp:
                    rp = float(rp[0])
                else:
                    break


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


