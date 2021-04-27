import numpy as np
import xarray as xr

class PropagationEnviornment:
    """General specification of a 2-D propagation enviornment"""

    def __init__(self, name='autogen'):
        """setup an empty framework"""
        self.name = name
        self.source = None
        self.p_acs = None

    def src(self, source_depth):
        """specify the source location"""
        source_depth = np.array(source_depth, ndmin=1)

        # use 1 as a value for source existance
        source_value = np.ones_like(source_depth)

        self.source = xr.DataArray(source_value,
                                   coords=[source_depth],
                                   dims=['depth'])

    #TODO: the pressure observations are still tied to the profile, which may not be a bad thing but should be either made explicit or changed.
    def pressure(self, obs_depth, obs_range):
        """specify the observation locations"""
        obs_depth = np.array(obs_depth, ndmin=1)
        obs_range = np.array(obs_range, ndmin=1)

        # use 0 as an empty placeholder
        p_val = np.zeros((obs_depth.size, obs_range.size))

        self.source = xr.DataArray(p_val,
                                   coords=[obs_depth, obs_range],
                                   dims=['depth', 'range'])


if __name__ == "__main__":
    test = PropagationEnviornment('hello')
    test.src(50.)
