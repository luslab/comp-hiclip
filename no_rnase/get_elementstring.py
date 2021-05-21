#!/usr/bin/env python

# import sys
import forgi.graph.bulge_graph as fgb
import re
# import sys

# db = '..((((((.....(((((((((.,,..(((........))).,}...))))))))).......))))))..{{{.....}}}..................'
db = input()
# db = sys.argv[1]
db = re.sub(',|{|}', '.', db) # Remove low probability matches

bg = fgb.BulgeGraph.from_dotbracket(db)


# bg = fgb.BulgeGraph.from_dotbracket('..((((..).)))..')

# bg = fgb.BulgeGraph.from_dotbracket('..((((((.....(((((((((.,,..(((........))).,}...))))))))).......))))))..{{{.....}}}..................')
# bg = fgb.BulgeGraph.from_dotbracket('..((((((.....(((((((((.....(((........)))......))))))))).......))))))..(((.....)))..................')
print(bg.to_bg_string())

# es = bg.to_element_string(with_numbers=False)

# es.find('h')

# bg.get_domains()

# f'{db}\n{es}'
# print(f'{db}\n{es}')
# print(db + '\n' + es)

# for s in bg.stem_iterator():
#     print(bg.stem_length(s))

# bg.get_connected_residues(s1, s2, bulge=None)
# bg.traverse_graph()
# bg.junctions()