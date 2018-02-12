import numpy as np


def distance(x, y, x1, y1):
	return math.sqrt((x - x1)**2 + (y - y1)**2)



def inbounds(x, y):
	return x < 2048 and x >= 0 and y < 1489 and y >= 0


def objid_extract(obj_id, full=False):
    #use to extract id information from objid.
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



def floodfill(data, x, y, mean, threshold, visited):
	if x >= 0 and x < len(data[0]) and y >= 0 and y < len(data) and data[y][x] >= threshold and visited[y][x] == False:
		data[y][x] = mean
		visited[y][x] = True
	else:
		return data

	floodfill(data, x - 1, y, mean, threshold, visited)
	floodfill(data, x + 1, y, mean, threshold, visited)
	floodfill(data, x, y - 1, mean, threshold, visited)
	floodfill(data, x, y + 1, mean, threshold, visited)

	return data



def checkInner(data, sources, qX, qY, mean, stddev, pointer, gthresh):


	for i in np.arange(source.shape[1]):



		if distance(sources['COLC'][i][pointer], sources['ROWC'][i][pointer], qX, qY) > 5 \
             and distance(sources['COLC'][i][pointer], sources['ROWC'][i][pointer], qX, qY) < 200:
			sx = int(sources['COLC'][i][pointer])
			sy = int(sources['ROWC'][i][pointer])




			visited = np.zeros((len(data), len(data[0])), dtype=bool)


			data = floodfill(data, sx, sy, mean, mean + stddev, visited)

	return data



# Specifically checks for cut-off bright sources on the boundaries of the image i.e. centroid not on cutout

def perimeter(data, mean, stddev):
	visited = np.zeros((len(data), len(data)), dtype=bool)
	for i in range(len(data)):
		if data[len(data) - 1][i] > mean + 3 * stddev:
			data = floodfill(data, i, len(data) - 1, mean, mean + stddev, visited)
		if data[i][len(data) - 1] > mean + 3 * stddev:
			data = floodfill(data, len(data) - 1, i, mean, mean + stddev, visited)
		if data[0][i] > mean + 3 * stddev:
			data = floodfill(data, i, 0, mean, mean + stddev, visited)
		if data[i][0] > mean + 3 * stddev:
			data = floodfill(data, 0, i, mean, mean + stddev, visited)

	return data



# Method that checks if a source > 1 sigma is outside of a certain radius, and if so, masks it by putting mean value

def checkOutter(data, mean, std):
	count = 0
	for i in range(len(data)):
		for j in range(len(data[0])):
			if data[i][j] > mean + 3 * std and distance(i, j, 50, 50) > 50:
				data[i][j] = mean
	return data


def radiusmask(data, xc, yc, isophotal_radius, mean, visited):
	if xc < 0 or xc >= len(data) or yc < 0 or yc >= len(data[0]):
		return data

	for r in range(len(data)):
		for c in range(len(data)):
			if distance(r, c, xc, yc) < 1.2 * isophotal_radius / 2 + 3:
				data[c][r] = mean

	return data

# Stand in for Atlas files, determines which pixels are associated with each source by checking if greater than threshold (2 sigma)
# Returns a boolean array

def connected(data, x, y, threshold, visited, radius):
	if x >= 0 and x < len(data) and y >= 0 and y < len(data) and data[y][x] >= threshold and \
        visited[y][x] == False and distance(x, y, 50, 50) <= radius:
		visited[y][x] = True
	else:
		return visited

	connected(data, x - 1, y, threshold, visited, radius)
	connected(data, x + 1, y, threshold, visited, radius)
	connected(data, x, y - 1, threshold, visited, radius)
	connected(data, x, y + 1, threshold, visited, radius)

	return visited


# Calculate the background by taking the sigma clipped median value between the radii specified

def calc_background(data, x, y, radius1, radius2):
	bkg_array = []
	for i in range(len(data[0])):
		for j in range(len(data)):
			if abs(x - i) < radius2 and abs(y - j) < radius2 and inbounds(i, j) and (i - x)**2 + \
                         (j - y)**2 >= radius1**2 and (i - x)**2 + (j - y)**2 <= radius2**2:
				bkg_array.append(data[j, i])

	true_bkg = sigma_clipped_stats(bkg_array, sigma=3.0, iters=5)[0]
	print(true_bkg)
	return true_bkg


# Method that checks if a particular residue is noise or not by drawing a ring of similar distance as the "bright" point around the centroid

def checkNoise(x, y, qX, qY, data):
	halo = []
	for i in range(42):
		for j in range(42):
			if abs(distance(i, j, qX, qY) - distance(x, y, qX, qY)) <= 2 and distance(i, j, x, y) >= 2:
				halo.append(data[j, i])
	mean, median, std = sigma_clipped_stats(halo, sigma=3.0, iters=5)

	if data[y, x] > mean + 3 * std:
		return True
	return False



# Method that calculates the total photon count within 3 sigma of the quasar centroid

def photonCount(xc, yc, sigma, data):
	count = 0
	for i in range(len(data)):
		for j in range(len(data)):
			if distance(i, j, xc, yc) <= sigma:
				count += data[i][j]
	return count





# Normalizes a PSF by recursively checking all connecting bright points greater than a specified and treshold (2 sigma) and scaling them

def normalize(data, x, y, p1, p2, threshold, visited):
	if x >= 0 and x < len(data) and y >= 0 and y < len(data) and data[y][x] >= threshold and visited[y][x] == False:
		data[y][x] /= p1
		data[y][x] *= p2
		visited[y][x] = True
	else:
		return data

	normalize(data, x - 1, y, p1, p2, threshold, visited)
	normalize(data, x + 1, y, p1, p2, threshold, visited)
	normalize(data, x, y - 1, p1, p2, threshold, visited)
	normalize(data, x, y + 1, p1, p2, threshold, visited)

	return data



# Determine all flags associated with an object

def detflags(total):
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



        # Interpolate the image to the centroid of the QSO
        image_cut = image[lower_y : upper_y, lower_x : upper_x]

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
