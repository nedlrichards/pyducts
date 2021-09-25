from numpy.distutils.misc_util import Configuration
from numpy.distutils.core import setup

def configuration(parent_package='',top_path=None):
    config = Configuration('ram_pe', parent_package, top_path)
    sources = ['src/pade_coeffs.f90',
               'src/pade_coeffs_wrapper.f90',
               'src/ram_v4.f90',
               'src/ram_v4_wrappper.f90']
    config.add_library('ram_src', sources=sources)
    config.add_extension('ram_pe', sources=['src/epade.c', 'src/main_ram_v4.c'], libraries=['ram_src'])
    return config

if __name__ == '__main__':
    setup(configuration=configuration)
