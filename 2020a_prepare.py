import saga_salt

base = 'SALT_2020-1-SCI-034/product/'

sciarcs = [
(base+'mbxgpP202006240064.fits', base+'mbxgpP202006240065.fits'), # 208-pgc66934-201237793
(base+'mbxgpP202006190134.fits', base+'mbxgpP202006190135.fits'), # 002-pgc66934-203016260
(base+'mbxgpP202006180107.fits', base+'mbxgpP202006180108.fits'), # 091-pgc66318-188047679
(base+'mbxgpP202006170059.fits', base+'mbxgpP202006170060.fits'), # 157-pgc67782-219143443
(base+'mbxgpP202005290034.fits', base+'mbxgpP202005290035.fits'), # 134-pgc66318-188048484
#(base+'mbxgpP202005210063.fits', base+'.fits'), # 091-pgc66318-188047679   ... focus did not go back properly after a move and it drifted. We didn't notice. Reject.
(base+'mbxgpP202005210064.fits', base+'mbxgpP202005210065.fits'), # 040-pgc67782-247125749
#(base+'mbxgpP202005190077.fits', base+'mbxgpP202005190078.fits'), # 017-pgc66318-192117557
]

saga_salt.prepare_reduction(sciarcs, clobber=True, basepath='2020a', pathtoarcs='/home/pysalt/linelists/')