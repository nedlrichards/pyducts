import numpy as np
import xarray as xr

class PropagationEnviornment:
    """General specification of a 2-D propagation enviornment"""

    def __init__(self, name='autogen'):
        """setup an empty framework"""
        self.name = name
        self.source = xr.DataArray()
        self.p_acs = xr.DataArray()
        self.profile = xr.DataArray()

    def src(self, source_depth, frequency):
        """specify the source location"""
        source_depth = np.array(source_depth, ndmin=1)
        frequency = np.array(frequency, ndmin=1)
        self.source = xr.DataArray(1+0j,
                                   dims=['source_depth', 'frequency'],
                                   coords={'source_depth': source_depth,
                                           'frequency':frequency},
                                   attrs={'units':'(m), (Hz)'})

    def pressure(self, dz, num_z, dr, num_r):
        """specify the observation locations"""
        self.p_acs = xr.DataArray(np.nan + 1j * np.nan,
                                  coords={'dz':dz,
                                          'num_z':num_z,
                                          'dr':dr,
                                          'num_r':num_r})

    def profiles(self, zs, rs, cp, interpolation=1):
        """Populate a profile matrix"""
        zs = np.array(zs, ndmin=1)
        rs = np.array(rs, ndmin=1)
        cp = np.array(cp, ndmin=1)

        self.profile = xr.DataArray(cp,
                                    dims=['zs', 'rs'],
                                    coords={'zs':zs,
                                            'rs':rs},
                                    attrs={'units':'(m)'})

    def to_netCDF(self):
        """Put all data into a dataset and convert to netCDF file"""
        all_data = xr.Dataset(data_vars={'source':self.source,
                                         'p_acs':self.p_acs,
                                         'profiles':self.profile})
        all_data.to_netcdf(self.name + '.nc', format="NETCDF4")


if __name__ == "__main__":
    test = PropagationEnviornment('hello')
    test.src(50., 1e3)
