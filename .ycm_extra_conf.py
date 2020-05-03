
TDSEHOME = __file__
from os.path import join, exists
assert exists(TDSEHOME)
TDSEINC = join(TDSEHOME, "include")

from os import getenv
MATRIXHOME = getenv('MATRIXHOME')
assert exists(MATRIXHOME)
MATRIX_INC = join(MATRIXHOME, "include")

PARAMHOME = getenv('PARAMHOME')
assert exists(PARAMHOME)
PARAM_INC = join(PARAMHOME, "include")

def Settings( **kwargs ):
  return {
    'flags': [ 
        '-x', 'c++', '-Wall', '-Wextra', '-Werror',
        '-I', TDSEINC, 
        '-I', MATRIX_INC,
        '-I', PARAM_INC,
#        '-DCOULOMB',  # add this flag when testing Coulomb potential test code
        ],
  }

