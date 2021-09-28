from numpy.distutils.misc_util import Configuration
from numpy.distutils.core import setup

def configuration(parent_package='',top_path=None):
    config = Configuration('epade', parent_package, top_path)
    sources = ['src/constants.f90',
               'src/cmplx_roots_sg.f90',
               'src/pade_coeffs.f90',
               'src/pade_coeffs_wrapper.f90',
               'src/ram_io.f90',
               'src/ram_v4.f90',
               'src/ram_v4_wrapper.f90']
    config.add_library('ram_src', sources=sources)
    config.add_extension('epade',
                         #sources=['src/epade.c',],
                         sources=['src/ram_main_v4.c',],
                         #sources=['src/epade.c',
                                  #'src/ram_main_v4.c'],
                         libraries=['ram_src'])
    return config

if __name__ == '__main__':
    setup(configuration=configuration)
