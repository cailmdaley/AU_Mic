from glob import glob

files = glob("data_files/*.ms")
tclean(vis=files,
imagename='aumic_composite_dirty',
imsize=512,
cell='0.03arcsec',
niter=1000,
usemask='auto-thresh')
