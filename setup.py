from setuptools import setup
#from distutils.extension import Extension
#from Cython.Distutils import build_ext

setup(
    name="WMM",
    version="0.1.0dev",
    #cmdclass = {'build_ext': build_ext},
    #ext_modules = [Extension("wmm.grid",
    #    sources=["wmm/grid.pyx", "GeomagnetismLibrary.c"],
    #    include_dirs=["."],
    #)],
    packages=["wmm"],
    package_data={'wmm': ['wmm/WMM.COF', "wmm/wmm_grid.exe"]},
    include_package_data=True
)
