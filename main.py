#!/usr/bin/env python3
from sys import argv


from drawer import drawObj

filename = argv[1] if len(argv) > 1 else "obj/duck.obj"

print(filename)
drawObj(filename)
