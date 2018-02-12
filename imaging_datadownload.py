import numpy as np
from astropy.io import fits
import urllib
from urllib.request import urlretrieve


color = sys.argv[1]
color_array = ['u', 'g', 'r', 'i', 'z']

pointer = color_array.index(color)

data_set = sys.argv[2]


if data_set == 'DR7':

    Q_catalog = fits.open('/Volumes/work_dir/catalog/QSOs_DR7_by_Schneider.fits')
    Mg_catalog = fits.open('/Volumes/work_dir/catalog/Trimmed_SDSS_DR7_107.fits')

    Q_data = Q_catalog[1].data
    #print Q_catalog[1].columns
    Mg_data = Mg_catalog[1].data



    #run rerun col run col filld
    run = Q_data['RUN']
    rerun = Q_data['RERUN']
    col = Q_data['CAMCOL']
    field = Q_data['FIELD']



    for i in np.arange(run.size):
        #run-rerun-col-run-col-field
        download_link_image = 'http://das.sdss.org/imaging/%d/%d/corr/%d/fpC-%06d-'+pointer+'%d-%04d.fit.gz'% \
                         (int(run[i]), int(rerun[i]), int(col[i]), int(run[i]), int(col[i]), int(filed[i]))

        donwload_link_table = 'http://das.sdss.org/imaging/%d/%d/objcs/%d/fpObjc-%06d-'+pointer+'%d-%04d.fit.gz' % \
                          (int(run[i]), int(rerun[i]), int(col[i]), int(run[i]), int(col[i]), int(filed[i]))



        out_image_name = 'fpC-%06d-'+pointer+'%d-%04d-%d.fit' % ( int(run[i]), int(col[i]), int(filed[i]), int(rerun[i]))
        out_table_name = 'fpObjc-%06d-'+pointer+'%d-%04d-%d.fit' % ( int(run[i]), int(col[i]), int(filed[i]), int(rerun[i]))




        urlretrieve(download_link_image, out_image_name)
        urlretrieve(download_link_table, out_table_name)



if data_set == 'DR12':

    Q_catalog = fits.open('/Volumes/work_dir/catalog/DR12Q.fits')
    Mg_catalog = fits.open('/Volumes/work_dir/catalog/Trimmed_BOSS_DR12_107.fits')

    Q_data = Q_catalog[1].data
    #print Q_catalog[1].columns
    Mg_data = Mg_catalog[1].data


    #run rerun col run col filld
    run = Q_data['RUN_NUMBER']
    rerun = Q_data['RERUN_NUMBER']
    col = Q_data['COL_NUMBER']
    field = Q_data['FIELD_NUMBER']


    data_path = '/image_quasar/DR12/'


    for i in np.arange(run.size):

        download_link_image = 'https://dr12.sdss.org/sas/dr12/boss/photoObj/frames/301/%d/%d/frame-g-%06d-%d-%04d.fits.bz2' % \
                              (int(run[i]), int(col[i]), int(run[i]), int(col[i]), int(filed[i]))

        download_link_table = 'https://data.sdss.org/sas/dr12/boss/photo/redux/301/%d/objcs/%d/fpObjc-%06d-%d-%04d.fit' % \
                              (int(run[i]), int(col[i]), int(run[i]), int(col[i]), int(filed[i]))



        out_image_name = 'fpC-%06d-'+pointer+'%d-%04d-%d.fit' % ( int(run[i]), int(col[i]), int(filed[i]), int(rerun[i]))
        out_table_name = 'fpObjc-%06d-'+pointer+'%d-%04d-%d.fit' % ( int(run[i]), int(col[i]), int(filed[i]), int(rerun[i]))


        urlretrieve(download_link_image, out_image_name)
        urlretrieve(download_link_table, out_table_name)
