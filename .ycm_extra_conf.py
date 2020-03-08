
#TDSEHOME = getenv('TDSEHOME')
TDSEHOME = __file__
from os.path import join, exists
#TDSEHOME = "/home/ahn/Dropbox/tdse/numerical/package/tdse"
assert exists(TDSEHOME)
TDSEINC = join(TDSEHOME, "include")
#assert exists(TDSEINC)

from os import getenv
MATRIXHOME = getenv('MATRIXHOME')
assert exists(MATRIXHOME)
MATRIX_INC = join(MATRIXHOME, "include")

def Settings( **kwargs ):
  return {
    'flags': [ 
        '-x', 'c++', '-Wall', '-Wextra', '-Werror',
        '-I', TDSEINC, 
        '-I', MATRIX_INC ],
  }

