#!/usr/bin/env python
__author__    = 'Hoby Rakotoarivelo'
__license__   = 'GPL v3'
__version__   = '1.0'
__copyright__ = 'Copyright 2016, The trinity project'

import os

def replace_ext(path, ext):
    (base,ext) = os.path.splitext(path)
    return base+".bb"

def handle_nodes(stream,mesh):
    #
    nb_nodes = int(stream.readline().split()[0])
    #
    for i in range(nb_nodes):
        coords = [float(c) for c in stream.readline().split()]
        label  = int(coords.pop())
        mesh["nodes"][i] = tuple(coords)
    return stream
#
def handle_elems(stream,mesh):
    nb_elems = int(stream.readline().split()[0])

    for i in range(nb_elems):
        vertices = [int(k) for k in stream.readline().split()]
        mesh["elems"][i] = tuple(vertices)
    return stream
#
def handle_default(stream,mesh):
    return stream
#
def handle_end(stream,mesh):
    return None

#
_keywords = {
    'Dimension'     : handle_default,
    'Vertices'      : handle_nodes,
    'Edges'         : handle_default,
    'Triangles'     : handle_elems,
    'Quadrangles'   : handle_default,
    'Quadrilaterals': handle_default,
    'Tetrahedra'    : handle_default,
    'End'           : handle_end
}

#
def parse_medit(path):

    global _keywords;
    #
    mesh = {'dim':2,'nodes':{},'darts':{},'elems':{}}
    # parse input mesh
    with open(path,'r') as f:
        while f:
            line = f.readline().split()
            if line and line[0] in _keywords.keys():
                f = _keywords[line[0]](f,mesh)

    return mesh

