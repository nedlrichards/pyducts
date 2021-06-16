import numpy as np
import json

class Medium:
    """A depth dependent profile of ocean acoustic properties"""
    def __init__(self, z, cp, rho=1000., cs=None, attn_p=None,
                 attn_s=None, order=0):
        """define the depths of interfaces"""
        self.properties = ['z', 'cp', 'rho']
        self.order = order
        z = np.array(z)
        cp = np.array(cp)
        rho = np.array(rho)
        ssp = [z]

        if cp.size == 1 and z.size > 1:
            cp = np.full_like(z, cp.item())
        ssp += [cp]

        if rho.size == 1 and z.size > 1:
            rho = np.full_like(z, rho.item())
        ssp += [rho]

        if cs is not None:
            self.properties += ['cs']
            cs = np.array(cs)
            if cs.size == 1 and z.size > 1:
                cs = np.full_like(z, cs.item())
            ssp += [cs]

        if attn_p is not None:
            self.properties += ['attn_p']
            attn_p = np.array(attn_p)
            if attn_p.size == 1 and z.size > 1:
                attn_p = np.full_like(z, attn_p.item())
            ssp += [attn_p]

        if attn_s is not None:
            self.properties += ['attn_s']
            attn_s = np.array(attn_s)
            if attn_s.size == 1 and z.size > 1:
                attn_s = np.full_like(z, attn_s.item())
            ssp += [attn_s]

        self.z = z
        self.cp = cp
        self.rho = rho
        self.cs = cs
        self.attn_p = attn_p
        self.attn_s = attn_s
        self.ssp = np.array(ssp)

    def __repr__(self):
        rstr = 'medium\n\n' \
             + 'property list:\n{.properties}\n\n'.format(self) \
             + 'ssp:\n{.ssp.T}'.format(self)
        return rstr

    def dumps(self):
        # make stream as small as possible
        stream = json.dumps(water,
                            cls=MediumEncoder,
                            separators=(',', ':'))
        return stream


class MediumEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, Medium):
            # take care of optional arguments first
            z = obj.z if obj.z is None else obj.z.tolist()
            cp = obj.cp if obj.cp is None else obj.cp.tolist()
            rho = obj.rho if obj.rho is None else obj.rho.tolist()
            cs = obj.cs if obj.cs is None else obj.cs.tolist()
            attn_p = obj.attn_p if obj.attn_p is None else obj.attn_p.tolist()
            attn_s = obj.attn_s if obj.attn_s is None else obj.attn_s.tolist()
            # dictionary representation of medium
            i_dct = {"__medium__": True,
                     "properties": obj.properties,
                     "z": obj.z.tolist(),
                     "cp": obj.cp.tolist(),
                     "rho": rho,
                     "cs": cs,
                     "attn_p": attn_p,
                     "attn_s": attn_s
                     }
            return i_dct
        else:
            return json.JSONEncoder.default(self, obj)



if __name__ == '__main__':

    def decode_interface(dct):
        if "__medium__" in dct:
            return Medium(dct["z"],
                          dct["cp"],
                          rho=dct['rho'],
                          cs=dct["cs"],
                          attn_p=dct["attn_p"],
                          attn_s=dct["attn_s"],
                          )
        return dct

    water = Medium([100., 200.],  [1500., 1515.], cs=100.)
    sediment = Medium(200., 1800., rho=1800., cs=100., attn_s=0.1)
    data = water.dumps()
    with open("test.txt", 'w') as fid:
        json.dump([water, sediment], fid, cls=MediumEncoder)

    iface = json.loads(data, object_hook=decode_interface)
