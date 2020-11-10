from numpy.distutils.core import Extension

ext1 = Extension(name='kraken',
                 sources = ['./misc/subtabulate.f90',
                            './misc/FatalError.f90',
                            './misc/beampattern.f90',
                            './misc/MathConstants.f90',
                            './misc/RefCoef.f90',
                            './misc/SourceReceiverPositions.f90',
                            './misc/pchipMod.f90',
                            './misc/AttenMod.f90',
                            './misc/sspMod.f90',
                            './misc/RWSHDFile.f90',
                            './misc/interpolation.f90',
                            './misc/MergeVectorsMod.f90',
                            './misc/munk.f90',
                            './misc/ReadEnvironmentMod.f90',
                            './misc/SortMod.f90',
                            './misc/splinec.f90',
                            './misc/subtabulate.f90',
                            './misc/calculateweights.f90',
                            './misc/norms.f90',
                            './misc/cross_products.f90',
                            './misc/monotonicMod.f90',
                            './misc/PolyMod.f90',
                            './misc/PekRoot.f90',
                            './misc/RootFinderSecantMod.f90'
                     ])

ext2 = Extension(name='misc',
                 sources = ['./misc/subtabulate.f90',
                            './misc/FatalError.f90',
                            './misc/beampattern.f90',
                            './misc/MathConstants.f90',
                            './misc/PolyMod.f90',
                            './misc/norms.f90',
                            './misc/munk.f90',
                            './misc/RootFinderSecantMod.f90',
                            './misc/AttenMod.f90',
                            #'./misc/sspMod.f90',
                            #'./misc/ReadEnvironmentMod.f90',
                            #'./misc/calculateweights.f90',
                            #'./misc/RefCoef.f90',
                            #'./misc/SourceReceiverPositions.f90',
                            #'./misc/pchipMod.f90',
                            #'./misc/RWSHDFile.f90',
                            #'./misc/interpolation.f90',
                            #'./misc/MergeVectorsMod.f90',
                            #'./misc/SortMod.f90',
                            #'./misc/splinec.f90',
                            #'./misc/subtabulate.f90',
                            #'./misc/cross_products.f90',
                            #'./misc/monotonicMod.f90',
                            #'./misc/PekRoot.f90',
                     ])

                     #'./Kraken/kraken.f90'])

ext3 = Extension(name='attn',
                 sources = ['./misc/MathConstants.f90',
                            './misc/FatalError.f90',
                            './misc/AttenMod.f90'])


if __name__ == "__main__":
    from numpy.distutils.core import setup

    setup(
        name="pyducts",
        version="0.1",
        ext_modules = [ext3]
    )


