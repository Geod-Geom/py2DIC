# -*- coding: utf-8 -*-
'''
 py2DIC
 2D Digital Image Correlation software
 developed by Geodesy and Geomatics Division   
 Sapienza University of Rome
 
 The information in this file is
 Copyright(c) 2017, 
 Andrea Nascetti    <andrea.nascetti@uniroma1.it>,  
 Valeria Belloni    <valeria.belloni@uniroma1.it>,
 Roberta Ravanelli  <roberta.ravanelli@uniroma1.it>,
 Martina Di Rita    <martina.dirita@uniroma1.it> 
 and is subject to the terms and conditions of the
 GNU Lesser General Public License Version 2.1
 The license text is available from
 http://www.gnu.org/licenses/lgpl.html
 
 More information in the following scientific papers:
 Ravanelli R., Nascetti A., Di Rita M., Belloni V., Mattei D., Nisticò N., and Crespi M.: A new Digital Image Correlation software for displacements field measurement in structural applications, The International Archives of the Photogrammetry, Remote Sensing and Spatial Information Sciences, XLII-4/W2, 139-145,
 https://doi.org/10.5194/isprs-archives-XLII-4-W2-139-2017, 2017.
 Belloni V., Ravanelli, R., Nascetti, A., Di Rita, M., Mattei, D., and Crespi, M.: DIGITAL IMAGE CORRELATION FROM COMMERCIAL TO FOS SOFTWARE: A MATURE TECHNIQUE FOR FULL-FIELD DISPLACEMENT MEASUREMENTS,The International Archives of the Photogrammetry, Remote Sensing and Spatial Information Sciences, XLII-2, 91-95, 
 https://doi.org/10.5194/isprs-archives-XLII-2-91-2018, 2018. 
Belloni V., Ravanelli, R., Nascetti, A., Di Rita, M., Mattei, D., and Crespi, M.: py2DIC: A New Free and Open Source Software for Displacement and Strain Measurements in the Field of Experimental Mechanics, Sensors 2019, 19, 3832. https://doi.org/10.3390/s19183832, 2019
'''

import cv2
import numpy as np
from matplotlib import pyplot as plt
import glob
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import matplotlib.colorbar as mcolorbar
from scipy import signal as sg
import os
from PyQt5.QtCore import pyqtRemoveInputHook
from joblib import Parallel, delayed
from joblib import load, dump
import timeit
import multiprocessing
from scipy.ndimage import gaussian_filter
from astropy.convolution import convolve
from astropy.convolution import Gaussian2DKernel
#import matplotlib.ticker as ticker


def unwrap_self(arg, **kwarg):
    return matching_class.modify_array2(*arg, **kwarg)

    
def template_match(img_master, img_slave, method = 'cv2.TM_CCOEFF_NORMED', mlx = 1, mly = 1, show=True):    
                
                # Apply image oversampling 
                img_master = cv2.resize(img_master,None,fx=mlx, fy=mly, interpolation = cv2.INTER_CUBIC)
                img_slave  = cv2.resize(img_slave,None,fx=mlx, fy=mly, interpolation = cv2.INTER_CUBIC)

                res = cv2.matchTemplate(img_slave,img_master,eval(method))

                w, h = img_master.shape[::-1]    
                min_val, max_val, min_loc, max_loc = cv2.minMaxLoc(res)
                
                # Control if the method is TM_SQDIFF or TM_SQDIFF_NORMED, take minimum value
                if method in [cv2.TM_SQDIFF, cv2.TM_SQDIFF_NORMED]:
                    top_left = min_loc
                else:
                    top_left = max_loc 
                bottom_right = (top_left[0] + w, top_left[1] + h)

                # Retrieve center coordinates
                px = (top_left[0]+bottom_right[0])/(2.0*mlx)
                py = (top_left[1]+bottom_right[1])/(2.0*mly)

                # Scale images for visualization
                img_master_scaled = cv2.convertScaleAbs(img_master, alpha=(255.0/500))
                img_slave_scaled  = cv2.convertScaleAbs(img_slave, alpha=(255.0/500))

                cv2.rectangle(img_slave_scaled,top_left, bottom_right, 255, 2*mlx) 

                if show == True:
                    plt.figure(figsize=(20,10))
                    plt.subplot(131),plt.imshow(res,cmap = 'gray')
                    plt.title('Matching Result'), plt.xticks([]), plt.yticks([])
                    plt.subplot(132),plt.imshow(img_master_scaled,cmap = 'gray')
                    plt.title('Detected Point'), plt.xticks([]), plt.yticks([])
                    plt.subplot(133),plt.imshow(img_slave_scaled, cmap = 'gray')
                    plt.suptitle(method)
                    plt.show()
                
                return px, py, max_val

def convolve2d(slab,kernel,max_missing=0.5,verbose=True):

    '''2D convolution with missings ignored

    <slab>: 2d array. Input array to convolve. Can have numpy.nan or masked values.
    <kernel>: 2d array, convolution kernel, must have sizes as odd numbers.
    <max_missing>: float in (0,1), max percentage of missing in each convolution
                   window is tolerated before a missing is placed in the result.

    Return <result>: 2d array, convolution result. Missings are represented as
                     numpy.nans if they are in <slab>, or masked if they are masked
                     in <slab>.
    '''

    from scipy.ndimage import convolve as sciconvolve

    assert np.ndim(slab)==2, "<slab> needs to be 2D."
    assert np.ndim(kernel)==2, "<kernel> needs to be 2D."
    assert kernel.shape[0]%2==1 and kernel.shape[1]%2==1, "<kernel> shape needs to be an odd number."
    assert max_missing > 0 and max_missing < 1, "<max_missing> needs to be a float in (0,1)."

    # Get mask for missings
    if not hasattr(slab,'mask') and np.any(np.isnan(slab))==False:
        has_missing=False
        slab2=slab.copy()

    elif not hasattr(slab,'mask') and np.any(np.isnan(slab)):
        has_missing=True
        slabmask=np.where(np.isnan(slab),1,0)
        slab2=slab.copy()
        missing_as='nan'

    elif (slab.mask.size==1 and slab.mask==False) or np.any(slab.mask)==False:
        has_missing=False
        slab2=slab.copy()

    elif not (slab.mask.size==1 and slab.mask==False) and np.any(slab.mask):
        has_missing=True
        slabmask=np.where(slab.mask,1,0)
        slab2=np.where(slabmask==1,np.nan,slab)
        missing_as='mask'

    else:
        has_missing=False
        slab2=slab.copy()

    # No missing 
    if not has_missing:
        result=sciconvolve(slab2,kernel,mode='constant',cval=0.)
    else:
        H,W=slab.shape
        hh=int((kernel.shape[0]-1)/2)  # half height
        hw=int((kernel.shape[1]-1)/2)  # half width
        min_valid=(1-max_missing)*kernel.shape[0]*kernel.shape[1]

        # Do not forget to flip the kernel
        kernel_flip=kernel[::-1,::-1]

        result=sciconvolve(slab2,kernel,mode='constant',cval=0.)
        slab2=np.where(slabmask==1,0,slab2)

        # Get nan holes 
        miss_idx=zip(*np.where(slabmask==1))

        if missing_as=='mask':
            mask=np.zeros([H,W])

        for yii,xii in miss_idx:

            # Recompute at each new nan in result 
            hole_ys=range(max(0,yii-hh),min(H,yii+hh+1))
            hole_xs=range(max(0,xii-hw),min(W,xii+hw+1))

            for hi in hole_ys:
                for hj in hole_xs:
                    hi1=max(0,hi-hh)
                    hi2=min(H,hi+hh+1)
                    hj1=max(0,hj-hw)
                    hj2=min(W,hj+hw+1)

                    slab_window=slab2[hi1:hi2,hj1:hj2]
                    mask_window=slabmask[hi1:hi2,hj1:hj2]
                    kernel_ij=kernel_flip[max(0,hh-hi):min(hh*2+1,hh+H-hi),
                                     max(0,hw-hj):min(hw*2+1,hw+W-hj)]
                    kernel_ij=np.where(mask_window==1,0,kernel_ij)

                    # Fill with missing if not enough valid data 
                    ksum=np.sum(kernel_ij)
                    if ksum<min_valid:
                        if missing_as=='nan':
                            result[hi,hj]=np.nan
                        elif missing_as=='mask':
                            result[hi,hj]=0.
                            mask[hi,hj]=True
                    else:
                        result[hi,hj]=np.sum(slab_window*kernel_ij)

        if missing_as=='mask':
            result=np.ma.array(result)
            result.mask=mask

    return result

def gaussian_smooting_with_Nan(img, sigma=7, k_size=15):
    kernel = Gaussian2DKernel(sigma)
    img_smoothed = gaussian_filter(img, sigma=sigma, mode='nearest')
    img_smoothed_2 = convolve(sg.medfilt2d(img, k_size), kernel)

    idx_NaN = np.isnan(img)
    indx_NaN_filtered = np.isnan(img_smoothed)
    mask = idx_NaN != indx_NaN_filtered

    img_smoothed[mask] = img_smoothed_2[mask]
    img_smoothed[idx_NaN] = np.NaN
    img_smoothed = sg.medfilt2d(img_smoothed, 5)

    return img_smoothed

def create_circular_mask(h, w, center=None, radius=None):

    if center is None: # use the middle of the image
        center = [int(w/2), int(h/2)+20]
    if radius is None: # use the smallest distance between the center and image walls
        radius = min(center[0], center[1], w-center[0], h-center[1])

    Y, X = np.ogrid[:h, :w]
    dist_from_center = np.sqrt((X - center[0])**2 + (Y-center[1])**2)

    mask = dist_from_center <= radius
    return mask
 
def myfmt(x, pos):
    return '{0:.3f}'.format(x)

def myfmt2(x, pos):
    return '{0:2.1e}'.format(x)

def gif(path1,filename):

    imgs = []
    img_names = glob.glob(path1 + "/*.png")
    img_names.sort()
    for img_name in img_names:
           img = cv2.imread(img_name) 
           imgs.append(cv2.cvtColor(img, cv2.COLOR_BGR2RGB))

    writeGif(path1+filename, imgs, duration=1)

class matching_class:
    def __init__(self, A, B,temp_dim, b, c, d, Dim_y, Dim_x, h, w, H, V, step_y, step_x,number_of_processors, ml):
        self.A = A
        self.B = B
        self.temp_dim=temp_dim
        self.b = b
        self.c = c
        self.d = d
        self.Dim_y = Dim_y
        self.Dim_x = Dim_x
        self.h = h
        self.w = w
        self.H = H 
        self.V = V 
        self.step_y = step_y
        self.step_x = step_x
        self.results = np.zeros((self.H*self.V * self.w * self.h, 6))
        self.number_of_processors = number_of_processors
        self.ml = ml
        
    def run(self, num):
        #output = []#backend="threading"o multiprocess?
        # sommo sulf.results perche' potrei avrere i risultati del precedente livello temporale
        self.results[num] +=  Parallel(n_jobs= self.number_of_processors, backend="threading")\
            (delayed(unwrap_self)(i) for i in zip([self]*len(num), num))

        
    def modify_array2(self, k):

                    j = int(k / (self.H*self.V *self.w))
                    m = int(self.step_y * (k/ (self.w * self.H) - j * self.V))
                    n = int(self.step_x * (k % self.H))
                    i = int((k % (self.w * self.H))/ self.H)
                    # i, j sono coordinate griglia

                    #print (str(k))#, str(self.H*self.V*self.w*self.h)#, m, n, self.w, self.h, self.H, self.V
                    # occhio ai delta: se d e b non sono moltiplicati per 2, allora avra una sovrapposizione
                    # parziale delle search area (i template comunque non si sovrappongono mai, e restano sempre centrati sulla corrispondente search area)
                    # se non voglio che le search area si sovrappongano, devo moltiplicare b e d per due (anche nel ciclo for)
                    Delta_X = i*(self.c+2*self.d+self.temp_dim)# i*(c+d+temp_dim) aggiungendo il *2, le search area non si accavallano piu'
                    Delta_Y = j*(self.c+2*self.b+self.temp_dim)# j*(c+b+temp_dim)
                    TP_temp_x = Delta_X + self.c + self.d + self.temp_dim/2.0 + n # x start (pixel)
                    TP_temp_y = Delta_Y + self.c + self.b + self.temp_dim/2.0 + m # y start (pixel)

                    start_x_template_slice = self.c + self.d + Delta_X + n
                    stop_x_template_slice  = self.c + self.d + Delta_X + self.temp_dim + n 
                    start_y_template_slice = self.c + self.b + Delta_Y + m
                    stop_y_template_slice  = self.c + self.b + Delta_Y + self.temp_dim + m 
                    
                    shape_template = np.shape( self.A [  start_y_template_slice : stop_y_template_slice , start_x_template_slice : stop_x_template_slice ])
                   
                    # Check dimensions of template
                    assert np.allclose(shape_template[0], self.temp_dim)
                    assert np.allclose(shape_template[1], self.temp_dim)
                    assert np.allclose(stop_y_template_slice - start_y_template_slice, self.temp_dim)
                    assert np.allclose(stop_x_template_slice - start_x_template_slice, self.temp_dim)
                    assert np.allclose(TP_temp_x, (start_x_template_slice+ stop_x_template_slice)/2.0)
                    assert np.allclose(TP_temp_y, (start_y_template_slice+ stop_y_template_slice)/2.0)

                    start_x_search_slice = self.c + Delta_X + n 
                    stop_x_search_slice  = self.c + Delta_X + 2 * self.d + self.temp_dim + n 
                    start_y_search_slice = self.c + Delta_Y + m 
                    stop_y_search_slice  = self.c + Delta_Y + 2 * self.b + self.temp_dim + m 
                    
                    shape_search = np.shape( self.B [  start_y_search_slice : stop_y_search_slice , start_x_search_slice : stop_x_search_slice ])

                    # Check dimensions of search area
                    assert np.allclose((shape_search[0] - self.temp_dim) /2.0, self.b)
                    assert np.allclose((shape_search[1] - self.temp_dim) /2.0, self.d)
                    assert np.allclose(TP_temp_x, (start_x_search_slice+ stop_x_search_slice)/2.0)
                    assert np.allclose(TP_temp_y, (start_y_search_slice+ stop_y_search_slice)/2.0)

                    # in teoria dovrei fare un controllo anche sugli start e stop del template e del search, ma tanto stanno cmq piu a sinstra e piu' sopra dello stop del template
                    # if (stop_y_search_slice>=Dim_y) or (stop_x_search_slice>=Dim_x):
                    if (self.c + Delta_Y + 2*self.b +self.temp_dim + m >=self.Dim_y) or (self.c + Delta_X + 2*self.d + self.temp_dim + n >=self.Dim_x):
                        # mi devo fermare altrimenti vado fuori dai bordi esterni dell'immaginepass
                        '''output[k,0] = float('NaN')
                        output[k,1] = float('NaN')
                        output[k,2] = float('NaN')
                        output[k,3] = float('NaN')'''
                        return np.array(['NaN','NaN','NaN','NaN','NaN','NaN'])
                        
                    else:

                        temp  = self.A[start_y_template_slice :   stop_y_template_slice, start_x_template_slice : stop_x_template_slice]
                        search_area = self.B[start_y_search_slice : stop_y_search_slice,   start_x_search_slice : stop_x_search_slice]

                        indx, indy, maxcc = template_match(temp.astype('uint8'), search_area.astype('uint8'), mlx = self.ml, mly = self.ml, show = False)

                        TP_search_x = Delta_X + self.c + indx + n       # end point x [pixel]
                        TP_search_y = Delta_Y + self.c + indy + m       # end point y [pixel]

                        return np.array([TP_temp_x,TP_temp_y,TP_search_x-TP_temp_x,TP_search_y-TP_temp_y,0,0])


def DIC(images_absolute_path, dim_pixel, start_index, levels, image_time_sampling , temp_dim, b, d, recty1, recty2, rectx1, rectx2, H, V, defor, ml, c=0):
    start = timeit.default_timer()
    num_cores = multiprocessing.cpu_count()
    print ('used CPU cores:', num_cores)
    plt.close("all")
    plt.switch_backend("Qt5Agg")

    # Find the last folder, which we assume being the name of the test
    test_name = os.path.basename(os.path.normpath(images_absolute_path))

    initial_start_index = start_index
    print ('Test:', test_name)
    plt.rc('text', usetex=False)

    msg = test_name + "\n"
    msg = msg + '-----------\n'

    img_names = glob.glob( images_absolute_path +"/*")
    img_names.sort()

    # Read the first image to obtain information about the camera resolution 
    img1 = cv2.imread(img_names[0], 0)
    crop_img1 = img1[:,:]
    Dim_x = np.shape(crop_img1)[1] # width 
    Dim_y = np.shape(crop_img1)[0] # height 

    # If statement to define width and height in ubuntu
    windows = True
    if windows == False:
        #print ("UBUNTU")
        Dim_x = np.shape(crop_img1)[0] # width 
        Dim_y = np.shape(crop_img1)[1]

    ########### TEMPLATE PARAMETERS ################
    #H = 2*d + temp_dim # numero di griglie che si ripetono sulla singola riga della griglia grande
    #V = 2*b + temp_dim # numero di griglie che si ripetono sulla singola colonna della griglia grande
    ltot = H*V
    step_x = 1
    step_y = 1
    
    if recty1!= None and recty2!= None and rectx1!= None and rectx2!= None:
        # Otherwise the AOI is all the image
        Dim_y = recty2-recty1 
        Dim_x = rectx2-rectx1


    h = int(Dim_y/(temp_dim+2*b+c)-1)# qui c'era -1
    w = int(Dim_x/(temp_dim+2*d+c)-1)#
    grid_points = H*V * w * h

    for l in range(levels):

            stop_index = start_index + image_time_sampling
            msg = ''
            msg = msg + '-----------'
            msg = msg +  "\nAnalysed images"
            msg = msg +  "\n1: "+ os.path.basename(str(img_names[start_index]))
            msg = msg +  "\n2: "+ os.path.basename(str(img_names[stop_index]))
            print (msg)

            img1 = cv2.imread(img_names[start_index], 0)
            crop_img1 = img1[recty1:recty2, rectx1:rectx2]
            img2 = cv2.imread(img_names[stop_index], 0)
            crop_img2 = img2[recty1:recty2, rectx1:rectx2]
                       
            # If statement to define width and height in ubuntu
            if windows == True:
                img1 = cv2.imread(img_names[start_index], 0)
                crop_img1 = img1[recty1:recty2, rectx1:rectx2]
                img2 = cv2.imread(img_names[stop_index], 0)
                crop_img2 = img2[recty1:recty2, rectx1:rectx2]
            else:
                cv2.flip(crop_img1, 0, crop_img1)
                cv2.flip(crop_img2, 0, crop_img2)

                crop_img1 = crop_img1.T.copy()
                crop_img2 = crop_img2.T.copy() 
            
            print ('-----------')
            print ("Camera resolution")
            print (Dim_x, "x", Dim_y, "px")
            print ('-----------')

            msg = msg + '\n-----------\n'
            msg = msg + "Camera resolution\n"
            msg = msg + str(Dim_x)+"x"+ str(Dim_y)
            msg = msg + '\n-----------\n'
            img1r = crop_img1.copy()
            img2r = crop_img2.copy()
            
            assert np.allclose(img1r.shape, img2r.shape)
            assert np.allclose(img1r.shape[0], Dim_y)
            assert np.allclose(img1r.shape[1], Dim_x)
            classe_matching = matching_class(img1r, img2r,temp_dim, b, c, d, Dim_y, Dim_x, h, w, H, V, step_y, step_x, num_cores, ml)
            classe_matching.run(num = range(grid_points))

            start_index = stop_index
            results_mm = classe_matching.results.copy()* dim_pixel

            dx = results_mm[:,2].copy()
            dy = results_mm[:,3].copy()
            results_mm[:,4] = np.sqrt((results_mm[:,2])**2 + (results_mm[:,3])**2) 
            mod = results_mm[:,4].copy()

            dx.shape = (h*V, H*w)
            dy.shape = (h*V, H*w)

            if os.path.exists(os.path.split(images_absolute_path[:-1])[0] + "\\Mask"):

                mask_read = cv2.imread(os.path.split(images_absolute_path[:-1])[0] + "\\Mask\\mask.png", 0)
                mask1 = np.array(mask_read, dtype=bool)
                mask = mask1[int(b+(temp_dim-1)/2):int(b+(temp_dim-1)/2+h*V), int(d+(temp_dim-1)/2):int(d+(temp_dim-1)/2+ H*w)]

                dx[mask] = float("NaN")
                dy[mask] = float("NaN")

            dispx_smoothed = gaussian_smooting_with_Nan(dx,9,21)
            dispy_smoothed = gaussian_smooting_with_Nan(dy,9,21)

            if defor == True:

                kernel_der = np.matrix([1.0/280.0,-4.0/105.0,1.0/5.0, -4.0/5.0,0.0,4.0/5.0,-1.0/5.0,4.0/105.0,-1.0/280.0])
                strain_xx = convolve2d(dispx_smoothed,-kernel_der/1.0,0.5,True)
                strain_yy = convolve2d(dispy_smoothed,-kernel_der.T/1.0,0.5,True)

                du_dx = strain_xx
                dv_dy = strain_yy
                du_dy = convolve2d(dispx_smoothed,-kernel_der.T/1.0,0.5,True)
                dv_dx =convolve2d(dispy_smoothed,-kernel_der/1.0,0.5,True)

                strain_xx2 = du_dx + .5*(np.power(du_dx,2) + np.power(dv_dx,2))
                strain_xy2 = .5*(du_dy + dv_dx + du_dx*du_dy + dv_dx*dv_dy)
                strain_yy2 = dv_dy + .5*(np.power(dv_dy,2) + np.power(du_dy,2))

                plt.figure('strain_yy2')
                plt.imshow(strain_yy2, cmap=cm.jet)#, cmap=cm.jet, norm=mcolors.Normalize(vmin=0.002, vmax=0.012), interpolation = 'None')
                plt.axis('off')
                cb=plt.colorbar()
                plt.gca().set_aspect('equal', adjustable='box')

                plt.figure('strain_xx2')
                plt.imshow(strain_xx2, cmap=cm.jet)#, cmap=cm.jet, norm=mcolors.Normalize(vmin=-0.015, vmax=0.00), interpolation = 'None')
                plt.axis('off')
                cb=plt.colorbar()
                plt.gca().set_aspect('equal', adjustable='box')

                plt.figure('strain_xy2')
                plt.imshow(strain_xy2, cmap=cm.jet)#, cmap=cm.jet, norm=mcolors.Normalize(vmin=-0.005, vmax=0.005), interpolation = 'None')
                plt.axis('off')
                cb=plt.colorbar()
                plt.gca().set_aspect('equal', adjustable='box')


            with open("OutputPlots/"+test_name+"_results_mm_dy_"+str(initial_start_index)+"_"+str(stop_index)+".txt","w") as file_stats:
                        for i in range(np.shape(dy)[0]):
                            for j in range(np.shape(dy)[1]):
                                file_stats.write(str(dy[i,j])+'\t')
                            file_stats.write('\n')

            with open("OutputPlots/"+test_name+"_results_mm_dx_"+str(initial_start_index)+"_"+str(stop_index)+".txt","w") as file_stats:
                        for i in range(np.shape(dx)[0]):
                            for j in range(np.shape(dx)[1]):
                               file_stats.write(str(dx[i,j])+'\t')
                            file_stats.write('\n')

            fig = plt.figure("Displacements smoothed dx between img " + str(initial_start_index)+" and img "+str(stop_index))
            plt.title("Smoothed horizontal displacements: u")
            plt.gca().set_aspect('equal', adjustable='box')
            plt.imshow(crop_img1[int((temp_dim-1)/2+b)::, int((temp_dim-1)/2+d)::], cmap=plt.cm.gray, origin='upper')
            plt.imshow(dispx_smoothed, cmap=cm.jet, alpha = 0.5)#, norm=mcolors.Normalize(vmin=soglia_inf_x, vmax=soglia_sup_x), interpolation = 'None')

            cb=plt.colorbar()
            if dim_pixel == 1: 
                cb.set_label('px')
            else:
                cb.set_label('mm')                                

            fig = plt.figure("Displacements smoothed dy between img " + str(initial_start_index)+" and img "+str(stop_index)) 
            plt.title("Smoothed vertical displacements: v")
            plt.gca().set_aspect('equal', adjustable='box')
            plt.imshow(crop_img1[int((temp_dim-1)/2+b)::, int((temp_dim-1)/2+d)::], cmap=plt.cm.gray, origin='upper')
            plt.imshow(dispy_smoothed, cmap=cm.jet, alpha = 0.5)#, norm=mcolors.Normalize(vmin=soglia_inf_y, vmax=soglia_sup_y), interpolation = 'None')
            cb=plt.colorbar() 
            if dim_pixel == 1: 
                cb.set_label('px')
            else:
                cb.set_label('mm')  

            ########### QUIVER PLOT  ###########  
            results2 =  results_mm.copy()
            X = results2[:,0]
            Y = results2[:,1]
            disp_U = results2[:,2]
            disp_V = results2[:,3]

            #http://stackoverflow.com/questions/11970186/matplotlib-quiver-and-imshow-superimposed-how-can-i-set-two-colorbars
            #http://stackoverflow.com/questions/23964856/plotting-streamlines-with-matplotlib-python
            nz = mcolors.Normalize()
            nz.autoscale(results2[:,4]) 
            fig = plt.figure("img " + str(initial_start_index)+" - img "+str(stop_index))
            ax = fig.add_subplot(111)
            plt.imshow(crop_img1, cmap=plt.cm.gray, origin='upper')
            plt.title("img " + str(initial_start_index)+" - img "+str(stop_index))
            plt.gca().set_aspect('equal', adjustable='box') # Same scale on the x and y axes

            plt.ylabel('pixels')
            plt.xlabel('pixels')

            # Plot quiver
            plt.quiver((X[::100]/dim_pixel), (Y[::100]/dim_pixel), disp_U[::100]/mod[::100], disp_V[::100]/mod[::100], angles='xy', scale = 30, color=cm.jet(nz(mod[::100])))

            cax,_ = mcolorbar.make_axes(plt.gca())

            # soglia_sup_prova_mm = np.nanmax(results2[:,4])
            # soglia_inf_prova_mm = np.nanmin(results2[:,4])

            # # vmin and vmax should be symmetric? ex: - 6 ,6
            #cb = mcolorbar.ColorbarBase(cax, cmap=cm.jet)#, norm=mcolors.Normalize(vmin= soglia_inf_prova_mm, vmax= soglia_sup_prova_mm))
            #cb.set_clim(soglia_inf_prova_mm, soglia_sup_prova_mm)

            cb = mcolorbar.ColorbarBase(cax, cmap=cm.jet, norm=nz)

            ########### GRAFICO SEZIONE ########### 
            plt.figure("Vertical central section")
            plt.plot(dispy_smoothed[:, int(np.shape(dispy_smoothed)[1]/2)])

            plt.xlabel('y position (px)')
            if dim_pixel == 1: 
                plt.ylabel('y displacements (px)')
            else:
                plt.ylabel('y displacements (mm)')

    plt.show()
    
    stop = timeit.default_timer()
    print('Time [min]', (stop - start)/60) 

    return msg






