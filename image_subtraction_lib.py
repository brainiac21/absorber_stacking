import numpy as np
'''
this is the lib
'''


def distance(x, y, x1, y1):
	return math.sqrt((x - x1)**2 + (y - y1)**2)



def inbounds(x, y):
	return x < 2048 and x >= 0 and y < 1489 and y >= 0




def objid_extract(obj_id, full=False):
    '''
    use to extract id information from objid.
    '''

	masks={'sky_version':0x7800000000000000,
		   'rerun':0x07FF000000000000,
		   'run':0x0000FFFF00000000,
		   'camcol':0x00000000E0000000,
		   'first_field':0x0000000010000000,
		   'field':0x000000000FFF0000,
		   'id':0x000000000000FFFF}

	run=(obj_id & masks['run']) >> 32
	rerun=(obj_id & masks['rerun']) >> 48
	camcol=(obj_id & masks['camcol']) >> 29
	field=(obj_id & masks['field']) >> 16
	id=(obj_id & masks['id']) >> 0
	sky_version=(obj_id & masks['sky_version']) >> 59
	first_field=(obj_id & masks['first_field']) >> 28

	return {'run':run,
			'rerun':rerun,
			'camcol':camcol,
			'field':field,
			'id':id,
			'first_field':first_field,
			'sky_version':sky_version}


def floodfill(data, x, y, mean, threshold):
    '''
    recursive function to perform the floodfill algorithm to mask out objects with
    irregular shape.

    data: M * N image.
    x: starting pixel position in x axis
    y: starting pixel position in y axis
    mean: replaced value
    threshold: threshold value to perform floodfill.
    visited: maybe not necessary ?

    '''
	if x >= 0 and x < data.shape[0] and y >= 0 and y < data.shape[1] and data[y, x] >= threshold:
		data[y][x] = mean

	else:
		return data

	floodfill(data, x - 1, y, mean, threshold) #left
	floodfill(data, x + 1, y, mean, threshold) #right
	floodfill(data, x, y - 1, mean, threshold) #up
	floodfill(data, x, y + 1, mean, threshold) #down

	return data




def checkInner(data, sources, qX, qY, mean, stddev, pointer, gthresh):
    '''
    function to mask out possible sources between 5 and 200 pixels
    '''

	for i in np.arange(source.shape[1]):
        xs = int(sources['COLC'][i, pointer])
        ys = int(sources['ROWC'][i, pointer])


		if distance(xs,ys, qX, qY) > 5 and distance(xs,ys, qX, qY) < 200:


			data = floodfill(data, sx, sy, mean, mean + stddev)

	return data




def perimeter(data, mean, stddev):
    '''
    checks for cut-off bright sources on the boundaries of the image
    i.e. centroid not on cutout
    '''
    x_pixels = data.shape[0]
    y_pxiels = data.shape[1]

    for i in np.arange(y_pxiels):
		if data[x_pixels -1 , i] > mean + 3 * stddev:
			data = floodfill(data, i, x_pixels -1, mean, mean + stddev)
		if data[i, x_pixels - 1] > mean + 3 * stddev:
			data = floodfill(data, x_pixels - 1, i, mean, mean + stddev)
		if data[0, i] > mean + 3 * stddev:
			data = floodfill(data, i, 0, mean, mean + stddev)
		if data[i, 0] > mean + 3 * stddev:
			data = floodfill(data, 0, i, mean, mean + stddev)

	return data


def checkOutter(data, mean, std):
	'''
    Method that checks if a source > 1 sigma is outside of a certain radius
    and if so, masks it by putting mean value
    '''
	for i in np.arange(data.shape[1]):
		for j in np.arange(data.shape[0]):
			if data[i, j] > mean + 3 * std and distance(i, j, 50, 50) > 50:
			   data[i, j] = mean

	return data




def calc_background(data, x, y, radius1, radius2):
    '''
    Calculate the background by taking the sigma clipped median value
    between the radii specified
    '''
	bkg_array = []
	for i in np.arange(data.shape[0]):
		for j in np.arange(data.shape[1]):
			if abs(x - i) < radius2 and abs(y - j) < radius2 and inbounds(i, j) and (i - x)**2 + \
                         (j - y)**2 >= radius1**2 and (i - x)**2 + (j - y)**2 <= radius2**2:
				bkg_array.append(data[j, i])

	true_bkg = sigma_clipped_stats(bkg_array, sigma=3.0, iters=5)[0]

	return true_bkg





def detflags(total):
    '''
    Determine all flags associated with an object
    '''
	flags = np.zeros(32, dtype=bool)

	if total < 0:
		total += 2<<31
		flags[31] = True

	for i in range(1, 32):
		if total > (2 << (31 - i)):
			total -= 2 << (31 - i)
			flags[31 - i] = True

	return flags


def mask_algorithms(image, xc, yc, obj_table, chunk_size, pointer):

   if inbounds(xc + chunk_size + 6, yc + chunk_size + 6) and inbounds(xc - chunk_size - 5, yc - chunk_size - 5):
        upper_y = int(yc + chunk_size + 6)
        lower_y = int(yc - chunk_size - 5)
        upper_x = int(xc + chunk_size + 6)
        lower_x = int(xc - chunk_size - 5)


        image_cut = image[lower_y : upper_y, lower_x : upper_x]
        mean1, median1, std1 = sigma_clipped_stats(image_cut, sigma=3.0, iters=5)

        # Mask all sources that are deemed unlikely to be an absorber host galaxy
        image = checkInner(image, obj_table, xc, yc, mean1, std1, pointer, 100)

        # estimate the background here:
        mean1 = calc_background(image, xc, yc, 200, 250):



        # Interpolate the image to the centroid of the QSO
        image_cut = image[lower_y : upper_y, lower_x : upper_x]

        image_cut -= mean1

        spline = interpolate.interp2d(np.arange(lower_x, upper_x),
                                      np.arange(lower_y, upper_y),
                                      image_cut)
        xrang = np.arange(lower_x +5, upper_x - 5)
        yrang = np.arange(lower_y +5, upper_y - 5)

        if len(xrang) > 2 * chunk_size + 1:
            xrang = xrang[:-1]
        if len(yrang) > 2 * chunk_size + 1:
            yrang = yrang[:-1]

        shifted = spline(xrang, yrang)
        shifted = perimeter(shifted, mean1, std1)
        mean1, median1, std1 = sigma_clipped_stats(shifted, sigma=3.0, iters=5)
        shifted -= median1
    else:
        print 'this is a bad target.'
        shifted = 0.0

    return shifted
