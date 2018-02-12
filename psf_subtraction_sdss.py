import numpy as np
import sys
from image_subtraction_lib import *

color = sys.argv[1]
color_array = ['u', 'g', 'r', 'i', 'z']

pointer = color_array.index(color)



Q_catalog = fits.open('/Volumes/work_dir/catalog/DR12Q.fits')
Mg_catalog = fits.open('/Volumes/work_dir/catalog/Trimmed_BOSS_DR12_107.fits')

Q_data = Q_catalog[1].data
Mg_data = Mg_catalog[1].data


#run rerun col run col filld
run = Q_data['RUN_NUMBER']
rerun = Q_data['RERUN_NUMBER']
col = Q_data['COL_NUMBER']
field = Q_data['FIELD_NUMBER']



image_path = '/image_quasar/DR12/'

for i in np.arange(run.size):
    file_name  = 'frame-%06d-g%d-%04d-%d.fit' % ( int(run[i]), int(col[i]), int(filed[i]), int(rerun[i]))
    image_frame = fits.open(image_path+file_name)
    image = image_frame[0].data.astype(float)

    #fpobjc table here to extract the pixel positions:
    table = Table.read()
    q_id = objid_extract(int(Mg_data['OBJ_ID']), False)['id']

    quasar = obj_table[q_id - 1]



	quasar_xc = quasar['COLC'][pointer]
	quasar_yc = quasar['ROWC'][pointer]

    chunk_size = 50

# here: check the quasar zoom-in image quality and perform the masking algorithm.

    final_image =  mask_algorithms(image, quasar_xc, quasar_yc, obj_table, chunk_size, pointer)

    if final_image.size > 1.0:
        input_shape = star_array.shape
        star_array = []

#start to build the training images sample here:
        for j in range(len(obj_table)):
            sx = obj_table['COLC'][j][pointer]
    		sy = obj_table['ROWC'][j][pointer]
    		flags1 = detflags(obj_table['OBJC_FLAGS'][j])
    		flags2 = detflags(obj_table['OBJC_FLAGS2'][j])


            # select the stars to get the eigen-images
            criterion1 = obj_table['OBJC_TYPE'][j] == 6 and flags1[12] == False and flags1[17] == False and \
                          flags1[18] == False and flags2[27] == False

            criterion2 = distance(sx, sy, quasar_xc, quasar_yc) > 5 and obj_table['PSFMAG'][j][pointer] < 18

            criterion3 = obj_table['M_RR_CC'][j][pointer] > 0 and abs(obj_table['M_RR_CC'][j][pointer] - \
                            obj_table['M_RR_CC'][obj_id - 1][pointer]) < 0.1 * obj_table['M_RR_CC'][obj_id - 1][pointer]


            if criterion1 and criterion2 and criterion3:

                final_ref =  mask_algorithms(image, quasar_xc, quasar_yc, obj_table, chunk_size, pointer)

                final_ref /= obj_table['PSFFLUX'][j][pointer]
                final_ref *= obj_table['PSFFLUX'][obj_id - 1][pointer]


                star_array.append(final_ref.reshape(-1))

        star_array = np.array(star_array)
        print(np.shape(largearr))

#start to do the psf subtraction here:
#maybe 80% of the components?

        numcomp = star_array.shape[0]
        mean_vector = np.mean(star_array, axis =0)
        star_array -= mean_vector
        ipca = IncrementalPCA(n_components=numcomp)
        ipca.fit(star_array)
        ipca_comp = ipca.components_

        ipca_comp = ipca_comp.T
        final_image_1d = final_image.reshape(-1) - mean_vector
        coeff = np.dot(final_image_1d, ipca_comp)
        final_fit = np.dot(ipca_comp, coeff)
        final_fit += mean_vector
        final_fit = final_fit.reshape((input_shape))
        mean, median, stddev = sigma_clipped_stats(final_fit)
        final_fit -= median

        residue = final_image - final_fit
        mean, median, stddev = sigma_clipped_stats(residue, sigma=3.0, iters=5)
        residue -= median
