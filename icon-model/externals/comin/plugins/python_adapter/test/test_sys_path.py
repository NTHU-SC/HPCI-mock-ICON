from dummy_module import dummy
import sys

dummy()

# check for duplicate entries in sys.path
assert len(sys.path) == len(set(sys.path)), "Duplicate entry in sys.path detected!"
