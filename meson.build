project('XEphem', 'c',
  version : '4.1.0',
  default_options : [
    'warning_level=3',
    'b_lto=true',
    'optimization=2'
    ])

inc_dirs = include_directories('libastro', 'libip', 'liblilxml')

src=[
  'GUI/xephem/aavso.c',
  'GUI/xephem/annotmenu.c',
  'GUI/xephem/broadcast.c',
  'GUI/xephem/calmenu.c',
  'GUI/xephem/closemenu.c',
  'GUI/xephem/compiler.c',
  'GUI/xephem/coordsmenu.c',
  'GUI/xephem/datamenu.c',
  'GUI/xephem/db.c',
  'GUI/xephem/dbmenu.c',
  'GUI/xephem/earthmap.c',
  'GUI/xephem/earthmenu.c',
  'GUI/xephem/fallbacks.c',
  'GUI/xephem/favmenu.c',
  'GUI/xephem/formats.c',
  'GUI/xephem/fsmenu.c',
  'GUI/xephem/gallerymenu.c',
  'GUI/xephem/glance.c',
  'GUI/xephem/gsc.c',
  'GUI/xephem/gscnet.c',
  'GUI/xephem/helpmenu.c',
  'GUI/xephem/homeio.c',
  'GUI/xephem/hznmenu.c',
  'GUI/xephem/imregmenu.c',
  'GUI/xephem/indimenu.c',
  'GUI/xephem/jpeg2pm.c',
  'GUI/xephem/jupmenu.c',
  'GUI/xephem/listmenu.c',
  'GUI/xephem/mainmenu.c',
  'GUI/xephem/marsmenu.c',
  'GUI/xephem/marsmmenu.c',
  'GUI/xephem/moonmenu.c',
  'GUI/xephem/moviemenu.c',
  'GUI/xephem/msgmenu.c',
  'GUI/xephem/netmenu.c',
  'GUI/xephem/objmenu.c',
  'GUI/xephem/obslog.c',
  'GUI/xephem/patchlevel.c',
  'GUI/xephem/plot_aux.c',
  'GUI/xephem/plotmenu.c',
  'GUI/xephem/preferences.c',
  'GUI/xephem/progress.c',
  'GUI/xephem/ps.c',
  'GUI/xephem/query.c',
  'GUI/xephem/rotated.c',
  'GUI/xephem/satmenu.c',
  'GUI/xephem/saveres.c',
  'GUI/xephem/scope.c',
  'GUI/xephem/sites.c',
  'GUI/xephem/skybinary.c',
  'GUI/xephem/skyeyep.c',
  'GUI/xephem/skyfifos.c',
  'GUI/xephem/skyfiltmenu.c',
  'GUI/xephem/skyfits.c',
  'GUI/xephem/skyhist.c',
  'GUI/xephem/skyip.c',
  'GUI/xephem/skylist.c',
  'GUI/xephem/skytoolbar.c',
  'GUI/xephem/skyviewmenu.c',
  'GUI/xephem/solsysmenu.c',
  'GUI/xephem/splash.c',
  'GUI/xephem/srchmenu.c',
  'GUI/xephem/sunmenu.c',
  'GUI/xephem/time.c',
  'GUI/xephem/tips.c',
  'GUI/xephem/trailmenu.c',
  'GUI/xephem/ucac.c',
  'GUI/xephem/uranusmenu.c',
  'GUI/xephem/usno.c',
  'GUI/xephem/versionmenu.c',
  'GUI/xephem/webdbmenu.c',
  'GUI/xephem/xe2.c',
  'GUI/xephem/xe3.c',
  'GUI/xephem/xephem.c',
  'GUI/xephem/xmisc.c',
]

ccompiler = meson.get_compiler('c')

empty_dep = dependency('', required:false)

motif_dep = ccompiler.find_library('Xm')
xt_dep = ccompiler.find_library('Xt')
x_dep = dependency('X11')

openssl_dep = dependency('openssl')
m_dep = ccompiler.find_library('m', required: false)
if m_dep.found()
  add_project_link_arguments('-lm', language : 'c')

if get_option('with-sys-zlib')
  zlib_dep = dependency('zlib')
  libz = []
  zlib_inc = []
else
  zlib_dep = empty_dep
  subdir('libz')
  zlib_inc = include_directories('libz')
endif

if get_option('with-sys-jpeg')
  jpeg_dep = dependency('libjpeg')
  libjpeg = []
  jpeg_inc = []
else
  jpeg_dep = empty_dep
  subdir('libjpegd')
  jpeg_inc = include_directories('libjpegd')
endif

if get_option('with-sys-png')
  png_dep = dependency('libpng')
  libpng = []
  png_inc = []
else
  png_dep = empty_dep
  subdir('libpng')
  png_inc = include_directories('libpng')
endif

xephem_data_dir = get_option('datadir') / 'xephem'

subdir('libastro')
subdir('libip')
subdir('liblilxml')
subdir('GUI/xephem')

executable('xephem', src,
           dependencies: [motif_dep, xt_dep, x_dep, zlib_dep, \
               jpeg_dep, png_dep, openssl_dep, m_dep],
           include_directories: [inc_dirs, png_inc, zlib_inc, jpeg_inc],
           link_with: [libastro, libip, liblilxml, libpng, libz, libjpeg],
           install : true)
endif

run_target('cppcheck', command:['scripts/cppcheck.sh', '--project=' +
           join_paths(meson.build_root(), 'compile_commands.json')])

conf_data = configuration_data()
conf_data.set('DATADIR', get_option('prefix') / xephem_data_dir)
configure_file(input: 'XEphem.in',
               output : 'XEphem',
               configuration: conf_data,
               install_dir : get_option('sysconfdir') / 'X11/app-defaults/')
