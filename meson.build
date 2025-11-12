# SPDX-FileCopyrightText: 2023 Sunip K. Mukherjee
#
# SPDX-License-Identifier: Apache-2.0

project('glowpython2', 'c',
  version : '1.0',
  license: 'BSD-3',
  meson_version: '>=0.64.0',
  default_options : ['warning_level=2'],
)

add_languages('fortran')

# Find Python interpreter and dependency
py_mod = import('python')
py = py_mod.find_installation(pure: false)
py_dep = py.dependency()

# Define installation directory
install_dir = join_paths(py.get_install_dir(), 'glowpython2')

# Get numpy include directories
incdir_numpy = run_command(py,
  [
    '-c',
    'import numpy; print(numpy.get_include())'
  ],
  check: true
).stdout().strip()

# Get numpy.f2py include directory
incdir_f2py = run_command(py,
    ['-c', 'import numpy.f2py; print(numpy.f2py.get_include())'],
    check : true
).stdout().strip()

# Combine include directories
inc_np = include_directories(incdir_numpy, incdir_f2py)

# Install Python source files
py.install_sources(
    'src/glowpython2/__init__.py',
    'src/glowpython2/atmo_iri90.py',
    'src/glowpython2/atmo_msis00.py',
    'src/glowpython2/igrf.py',
    'src/glowpython2/pogo68.py',
    'src/glowpython2/base.py',
    'src/glowpython2/utils.py',
    'src/glowpython2/version.py',
    'src/glowpython2/plots.py',
    'src/glowpython2/examples.py',
    subdir: 'glowpython2',
    pure: false,
)

# Install data directories
install_subdir(
    'src/GLOW/data',
    install_dir: install_dir,
    strip_directory: false,
)

install_subdir(
    'src/ATMO/data',
    install_dir: install_dir,
    strip_directory: false,
)

install_subdir(
    'src/IGRF/data',
    install_dir: install_dir,
    strip_directory: false,
)

# GLOW Fortran sources
glow_fort_srcs = [
  'src/GLOW/cglow.f90', 
  'src/GLOW/glow.f90', 
  'src/GLOW/bands.f90', 
  'src/GLOW/conduct.f90',
  'src/GLOW/pyconduct.f90',
  'src/GLOW/egrid.f90', 
  'src/GLOW/ephoto_init.f90', 
  'src/GLOW/ephoto.f90', 
  'src/GLOW/etrans.f90', 
  'src/GLOW/exsect.f90',
  'src/GLOW/gchem.f90',
  'src/GLOW/maxt.f90', 
  'src/GLOW/qback.f90', 
  'src/GLOW/rcolum.f90',
  'src/GLOW/solzen.f90', 
  'src/GLOW/ssflux_init.f90', 
  'src/GLOW/ssflux.f90',
]

# GLOW Fortran F2PY wrapper
glowpython_source = custom_target('glowfortmodule.c',
  input : glow_fort_srcs,
  output : ['glowfortmodule.c', 'glowfort-f2pywrappers.f', 'glowfort-f2pywrappers2.f90'],
  command : [py, '-m', 'numpy.f2py', '@INPUT@', '-m', 'glowfort', '--lower']
)

# MSIS-00 and IRI-90 Fortran sources
atmo_fort_srcs = [
  'src/ATMO/iri90.f',
  'src/ATMO/nrlmsise00.f',
  'src/ATMO/dfp.f90',
  'src/ATMO/geomag.f90',
  'src/ATMO/snoem.f90',
]

# ATMO static library
atmo_static = static_library('atmo',
  atmo_fort_srcs,
  fortran_args : ['-cpp', '-march=native', '-ffast-math', '-std=legacy'],
  install : false,
)

# ATMO F2PY wrapper
atmoshim_tgt = custom_target('glowatmomodule.c',
  input : ['src/ATMO/atmo.f90'],
  output : ['glowatmomodule.c', 'glowatmo-f2pywrappers.f'],
  command : [py, '-m', 'numpy.f2py', '@INPUT@', '-m', 'glowatmo', '--lower']
)

# FIELDM Fortran sources
fieldm_fort_srcs = [
  'src/FIELDM/fieldm.f',
  'src/FIELDM/geomag.f90',
]

# FIELDM static library
fieldm_static = static_library('fieldm',
  fieldm_fort_srcs,
  fortran_args : ['-cpp', '-march=native', '-ffast-math', '-std=legacy'],
  install : false,
)

# FIELDM F2PY wrapper
fieldm_tgt = custom_target('glowfieldmmodule.c',
  input : ['src/FIELDM/magfield.f90'],
  output : ['glowfieldmmodule.c', 'glowfieldm-f2pywrappers.f'],
  command : [py, '-m', 'numpy.f2py', '@INPUT@', '-m', 'glowfieldm', '--lower']
)

# IGRF Fortran static library
igrf_static = static_library('igrf',
  ['src/IGRF/igrf.f', 'src/IGRF/dfp.f90'],
  fortran_args : ['-cpp', '-march=native', '-ffast-math', '-std=legacy'],
  install : false,
)

# IGRF F2PY wrapper
igrf_tgt = custom_target('glowigrfmodule.c',
  input : ['src/IGRF/igrf14_shim.f90'],
  output : ['glowigrfmodule.c', 'glowigrf-f2pywrappers2.f90'],
  command : [py, '-m', 'numpy.f2py', '@INPUT@', '-m', 'glowigrf', '--lower']
)

# GLOW Python extension module
py.extension_module('glowfort',
  glow_fort_srcs + [glowpython_source],
  incdir_f2py / 'fortranobject.c',
  include_directories: inc_np,
  dependencies : py_dep,
  install : true,
  install_dir: install_dir
)

# FIELDM Python extension module
py.extension_module('glowfieldm',
  ['src/FIELDM/magfield.f90', fieldm_tgt],
  incdir_f2py / 'fortranobject.c',
  include_directories: inc_np,
  dependencies : py_dep,
  link_with: [fieldm_static],
  install : true,
  install_dir: install_dir
)

# ATMO Python extension module
py.extension_module('glowatmo',
  ['src/ATMO/atmo.f90', atmoshim_tgt],
  incdir_f2py / 'fortranobject.c',
  include_directories: inc_np,
  dependencies : py_dep,
  link_with: [atmo_static],
  install : true,
  install_dir: install_dir
)

# IGRF Python extension module
py.extension_module('glowigrf',
  ['src/IGRF/igrf14_shim.f90', igrf_tgt],
  incdir_f2py / 'fortranobject.c',
  include_directories: inc_np,
  dependencies : py_dep,
  link_with: [igrf_static],
  install : true,
  install_dir: install_dir
)