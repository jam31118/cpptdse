from os.path import join, exists
TDSEHOME = "/home/ahn/Dropbox/tdse/numerical/package/tdse"
assert exists(TDSEHOME)
TDSEINC = join(TDSEHOME, "include")
assert exists(TDSEINC)

def Settings( **kwargs ):
  return {
    'flags': [ '-x', 'c++', '-Wall', '-Wextra', '-Werror', '-I', TDSEINC ],
  }

