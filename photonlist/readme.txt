boxer.so (exact name may vary by platform) needs to be created
from fortran boxer.f using command:

f2py -c -m boxer boxer.f

This command should work on any platform, I think, though the
output file may have a different file extension. In any case,
should all then be importable in python using
'import boxer'.
