
libip_src = [
  'explodegif.c',
  'fits.c',
  'fsmatch.c',
  'gaussfit.c',
  'lstsqr.c',
  'median.c',
  'sqr.c',
  'stars.c',
  'stats.c',
  'walk.c',
  'wcs.c',
]

libip_inc = include_directories('../libastro')
libip = static_library('ip', libip_src,
                       include_directories: libip_inc)
