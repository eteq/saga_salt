Procedure:

1. identify sci/arc pairs and add them to XXXX_prepare.py
2. run ``python XXX_prepare.py`` from inside pysalt docker container (requires access to arcs)
3. from pysalt docker, go into each to reduce, and do ``pyraf < XXXX_pysalt.txt``
4. Look at xcm...fits file, identify upper/lower bounds of galaxy trace.  Modify ``XXX_post_pysalt.py`` appropriately
5. Do ``python XXX_post_pysalt.py``
6. Inspect resulting png
7. Profit
