# -*- coding: utf-8 -*-
'''
http://www.admin-magazine.com/HPC/Articles/Parallel-Python-with-Joblib
https://stackoverflow.com/questions/45155781/numba-jit-slower-that-pure-python
 py2DIC
 2D Digital Image Correlation software
 developed by Geodesy and Geomatics Division   
 University of Rome "La Sapienza"
 
 The information in this file is
 Copyright(c) 2017, 
 Andrea Nascetti    <andrea.nascetti@uniroma1.it>,  
 Martina Di Rita    <martina.dirita@uniroma1.it>, 
 Roberta Ravanelli  <roberta.ravanelli@uniroma1.it>
 Valeria Belloni    <valeria.belloni@uniroma1.it>
 
 and is subject to the terms and conditions of the
 GNU Lesser General Public License Version 2.1
 The license text is available from
 http://www.gnu.org/licenses/lgpl.html
 
 If you use this software, please cite the following scientific paper:

 R. Ravanelli, A. Nascetti, M. Di Rita, V. Belloni, D. Mattei, N. Nisticò, and M. Crespi. (2017),
 A new Digital Image Correlation software for displacements field measurement in structural applications.
 The International Archives of the Photogrammetry, Remote Sensing and Spatial Information Sciences, XLII-4/W2:139–145,
 doi: 10.5194/ isprsarchives-XLII-4-W2-139-2017
 
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
#from images2gif import writeGif
#import prova_mask1 as c
from PyQt5.QtCore import pyqtRemoveInputHook
import tempfile
import shutil
from joblib import Parallel, delayed
from joblib import load, dump
import timeit
import multiprocessing
import pdb
from scipy import signal as sg
from scipy.ndimage import gaussian_filter
from astropy.convolution import convolve
from astropy.convolution import Gaussian2DKernel
import matplotlib.ticker as ticker


def unwrap_self(arg, **kwarg):
    #pdb.set_trace()
    #print (arg, kwarg)
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

    #--------------Get mask for missings--------------
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

    #--------------------No missing--------------------
    if not has_missing:
        result=sciconvolve(slab2,kernel,mode='constant',cval=0.)
    else:
        H,W=slab.shape
        hh=int((kernel.shape[0]-1)/2)  # half height
        hw=int((kernel.shape[1]-1)/2)  # half width
        min_valid=(1-max_missing)*kernel.shape[0]*kernel.shape[1]

        # dont forget to flip the kernel
        kernel_flip=kernel[::-1,::-1]

        result=sciconvolve(slab2,kernel,mode='constant',cval=0.)
        slab2=np.where(slabmask==1,0,slab2)

        #------------------Get nan holes------------------
        miss_idx=zip(*np.where(slabmask==1))

        if missing_as=='mask':
            mask=np.zeros([H,W])

        for yii,xii in miss_idx:

            #-------Recompute at each new nan in result-------
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

                    #----Fill with missing if not enough valid data----
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

class matching_class:
    def __init__(self, A, B,temp_dim, b, c, d, Dim_y, Dim_x, h, w, H, V, step_y, step_x,number_of_processors):
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
        
    def run(self, num):
        #output = []#backend="threading"o multiprocess?
        # sommo sulf.results perche' potrei avrere i risultati del precedente livello temporale
        self.results[num] +=  Parallel(n_jobs= self.number_of_processors, backend="threading")\
            (delayed(unwrap_self)(i) for i in zip([self]*len(num), num))
        #print( self.results)
        
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
                    
                    #print 'controllo template'
                    shape_template = np.shape( self.A [  start_y_template_slice : stop_y_template_slice , start_x_template_slice : stop_x_template_slice ])
                   
                    #controlli sulle dimensioni del template
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
                    
                    #print 'controllo search'
                    shape_search = np.shape( self.B [  start_y_search_slice : stop_y_search_slice , start_x_search_slice : stop_x_search_slice ])

                    # controllo sulle dimensioni dell'area di ricerca
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
                        # questo e' il pezzo critico che puo' causare l'outofbounds
                        temp  = self.A[start_y_template_slice :   stop_y_template_slice, start_x_template_slice : stop_x_template_slice]
                        search_area = self.B[start_y_search_slice : stop_y_search_slice,   start_x_search_slice : stop_x_search_slice]

                        indx, indy, maxcc = template_match(temp.astype('uint8'), search_area.astype('uint8'), mlx = 10, mly =10, show = False)

                        TP_search_x = Delta_X + self.c + indx + n       # end point x [pixel]
                        TP_search_y = Delta_Y + self.c + indy + m       # end point y [pixel]

                        return np.array([TP_temp_x,TP_temp_y,TP_search_x-TP_temp_x,TP_search_y-TP_temp_y,0,0])


def DIC(images_absolute_path,format, vel, dim_pixel, frame_rate, start_index, levels, image_time_sampling , temp_dim, b, d, recty1, recty2, rectx1, rectx2, H, V, c =0):
    start = timeit.default_timer()
    num_cores = multiprocessing.cpu_count()
    print ('used CPU cores:', num_cores)
    plt.close("all")
    plt.switch_backend("Qt5Agg")

    # find the last folder, which we assume being the name of the test
    test_name = os.path.basename(os.path.normpath(images_absolute_path))

    initial_start_index = start_index
    print ('dentro DIC',recty1, recty2, rectx1, rectx2, images_absolute_path)

    print('Backend: {}'.format(plt.get_backend()))
    print ('Test:', test_name)
    plt.rc('text', usetex=False)

    msg = test_name + "\n"
    msg = msg + '-----------\n'
    msg = msg + "Imposed deformation velocity\n"
    msg = msg + str(vel)+ " mm/min\n"

    vel = vel / 60 # mm/s

    format = '.'+format

    img_names = glob.glob( images_absolute_path +"/*"+format)
    img_names.sort()

    # Read the first image to obtain information about the camera resolution 
    img1 = cv2.imread(img_names[0], 0)
    crop_img1 = img1[:,:]
    Dim_x = np.shape(crop_img1)[1] # width 
    Dim_y = np.shape(crop_img1)[0] # height 

    # If statement to define width and height in ubuntu
    windows = False
    if Dim_x <= Dim_y:
        print ("Windows or MAC")
        windows =True
    else:
        print ("UBUNTU")
        Dim_x = np.shape(crop_img1)[0] # width 
        Dim_y = np.shape(crop_img1)[1]

    print ("Dim x:", Dim_x)
    print ("Dim y:", Dim_y)

    ########### TEMPLATE PARAMETERS ################
    # matchZoneWidth should be greater than the maximum displacement expected
    H = 2*d + temp_dim # numero di griglie che si ripetono sulla singola riga della griglia grande
    V = 2*b + temp_dim # numero di griglie che si ripetono sulla singola colonna della griglia grande
    ltot = H*V
    step_x = 1
    step_y = 1
    
    if recty1!= None and recty2!= None and rectx1!= None and rectx2!= None:
        # otherwise the AOI is all the image
        Dim_y = recty2-recty1 #img1r.shape[0]# sto supponenedo che tuttw le immagini abbiano la stessa shape
        Dim_x = rectx2-rectx1#img1r.shape[1]

    # dimensioni singola griglia
    h = np.int(Dim_y/(temp_dim+2*b+c)-1)# qui c'era -1
    w = np.int(Dim_x/(temp_dim+2*d+c)-1)#
    grid_points = H*V * w * h
    
    # Array where to store the results
    # we cumulate the displacements computed for each level
    results_array = np.zeros((grid_points, 6))
    print (results_array.shape)
    for l in range(levels):

            stop_index = start_index + image_time_sampling
            print (stop_index, start_index, image_time_sampling, initial_start_index)
            msg = ''
            print ("step", l)
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
            print ("Camera frame rate")
            print (frame_rate, "fps (one shot every", 1/frame_rate, "seconds)")

            msg = msg + '\n-----------\n'
            msg = msg + "Camera resolution\n"
            msg = msg + str(Dim_x)+"x"+ str(Dim_y)
            msg = msg + '\n-----------\n'
            msg = msg + "Camera frame rate\n"
            msg = msg + str(frame_rate)+ " fps (one shot every "+str(1/frame_rate)+" seconds)"
            img1r = crop_img1.copy()
            img2r = crop_img2.copy()
            
            assert np.allclose(img1r.shape, img2r.shape)
            assert np.allclose(img1r.shape[0], Dim_y)
            assert np.allclose(img1r.shape[1], Dim_x)
            classe_matching = matching_class(img1r, img2r,temp_dim, b, c, d, Dim_y, Dim_x, h, w, H, V, step_y, step_x, num_cores)
            classe_matching.run(num = range(grid_points))

            start_index = stop_index
            results_mm = classe_matching.results.copy()* dim_pixel

            dx = results_mm[:,2].copy()
            dy = results_mm[:,3].copy()

            dx.shape = (h*V, H*w)
            dy.shape = (h*V, H*w)

            # A=dx.shape
            # radius = 55
            # mask = create_circular_mask(A[0], A[1], radius=radius)

            # dx[mask] = float("NaN")
            # dy[mask] = float("NaN")

            dispx_smoothed = gaussian_smooting_with_Nan(dx,9,21)
            dispy_smoothed = gaussian_smooting_with_Nan(dy,9,21)

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
            # strain_xx2[mask] = float("NaN")
            # strain_xy2[mask] = float("NaN")
            # strain_yy2[mask] = float("NaN")

            # ####  PLOTS
            # soglia_inf_y = -12.15
            # soglia_sup_mm_y = -4.8#np.inf#levels*spost_atteso_pixel*dim_pixel*1.1
            # soglia_inf_x = -1.3
            # soglia_sup_mm_x =  0

            # ########### DISPLACEMENTS PLOT ########### 
            # indici_inf = np.where(dy < soglia_inf_y) # Consider only positive dy
            # ix = indici_inf [1] # x index
            # iy = indici_inf [0] # y index
            # dy[ iy, ix ] = float('NaN')

            # indici_sup = np.where(dy > soglia_sup_mm_y) # Consider only positive dy
            # ix = indici_sup [1] # x index
            # iy = indici_sup [0] # y index
            # dy[ iy, ix ] = float('NaN')
            
            # indici_inf = np.where(dx < soglia_inf_x) # Consider only positive dy
            # ix = indici_inf [1] # x index
            # iy = indici_inf [0] # y index
            # dx[ iy, ix ] = float('NaN')

            # indici_sup = np.where(dx > soglia_sup_mm_x) # Consider only positive dy
            # ix = indici_sup [1] # x index
            # iy = indici_sup [0] # y index
            # dx[ iy, ix ] = float('NaN')


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

            fig = plt.figure("Displacementsdx between img " + str(initial_start_index)+" and img "+str(stop_index))
            plt.title("Horizontal displacements: u")
            plt.gca().set_aspect('equal', adjustable='box')
            plt.imshow(dx, cmap=cm.jet)#, norm=mcolors.Normalize(vmin=soglia_inf_x, vmax=soglia_sup_mm_x), interpolation = 'None')

            cb=plt.colorbar() 
            #cb.set_clim(soglia_inf_x, soglia_sup_mm_x)
            cb.set_label('mm')

            plt.savefig("OutputPlots/"+test_name+"_displacements_dx_img_" + str(initial_start_index)+"_img_"+str(stop_index)+".png")                    

            fig = plt.figure("Displacementsdy between img " + str(initial_start_index)+" and img "+str(stop_index)) 
            plt.title("Vertical displacements: v")
            plt.gca().set_aspect('equal', adjustable='box')
            plt.imshow(dy, cmap=cm.jet)#, norm=mcolors.Normalize(vmin=soglia_inf_y, vmax=soglia_sup_mm_y), interpolation = 'None')
            cb2=plt.colorbar() 
            #cb2.set_clim(soglia_inf_y, soglia_sup_mm_y)
            cb2.set_label('mm')
            #mng = plt.get_current_fig_manager()
            plt.savefig("OutputPlots/"+test_name+"_displacements_img_" + str(initial_start_index)+"_img_"+str(stop_index)+".png")
            plt.savefig("GIF/"+test_name+"_displacements_dy_img_" + str(initial_start_index)+"_img_"+str(stop_index)+".png")

            fig = plt.figure("Displacements smoothed dx between img " + str(initial_start_index)+" and img "+str(stop_index))
            plt.title("Smoothed horizontal displacements: u")
            plt.gca().set_aspect('equal', adjustable='box')
            plt.imshow(dispx_smoothed, cmap=cm.jet)#, norm=mcolors.Normalize(vmin=soglia_inf_x, vmax=soglia_sup_mm_x), interpolation = 'None')
            cb=plt.colorbar() 
            #cb.set_clim(soglia_inf_x, soglia_sup_mm_x)
            cb.set_label('mm')                    

            fig = plt.figure("Displacements smoothed dy between img " + str(initial_start_index)+" and img "+str(stop_index)) 
            plt.title("Smoothed vertical displacements: v")
            plt.gca().set_aspect('equal', adjustable='box')
            plt.imshow(dispy_smoothed, cmap=cm.jet)#, norm=mcolors.Normalize(vmin=soglia_inf_y, vmax=soglia_sup_mm_y), interpolation = 'None')
            cb2=plt.colorbar() 
            #cb2.set_clim(soglia_inf_y, soglia_sup_mm_y)
            cb2.set_label('mm')

            plt.figure('strain_yy2')
            plt.imshow(strain_yy2, cmap=cm.jet)
            plt.gca().invert_yaxis()
            plt.clim(0.002,0.012)
            plt.axis('off')
            v1 = np.linspace(0.002,0.012, 8, endpoint=True)
            cb3=plt.colorbar(ticks=v1, format=ticker.FuncFormatter(myfmt2))
            plt.gca().set_aspect('equal', adjustable='box')

            plt.figure('strain_xx2')
            plt.imshow(strain_xx2, cmap=cm.jet)
            plt.gca().invert_yaxis()
            plt.clim(-0.015,0.0)
            plt.axis('off')
            v1 = np.linspace(-0.015,0.00, 8, endpoint=True)
            cb4=plt.colorbar(ticks=v1, format=ticker.FuncFormatter(myfmt))
            plt.gca().set_aspect('equal', adjustable='box')

            plt.figure('strain_xy2')
            plt.imshow(strain_xy2, cmap=cm.jet)
            plt.gca().invert_yaxis()
            plt.clim(-0.005,0.005)
            plt.axis('off')
            v1 = np.linspace(-0.005,0.005, 10, endpoint=True)
            cb4=plt.colorbar(ticks=v1, format=ticker.FuncFormatter(myfmt2))
            plt.gca().set_aspect('equal', adjustable='box')

            ########### QUIVER PLOT  ###########  
            results2 =  results_mm.copy()

            # Remove strains whose absolute value is greater than soglia_sup and lower than soglia_inf
            # limite_inf = np.where(results_mm2[:,3]<soglia_inf)[0]
            # results_mm2 = np.delete(results_mm2, limite_inf, axis=0)
            # limite_sup = np.where(results_mm2[:,3]>soglia_sup_mm)[0]
            # results_mm2 = np.delete(results_mm2, limite_sup, axis=0)

            #http://stackoverflow.com/questions/11970186/matplotlib-quiver-and-imshow-superimposed-how-can-i-set-two-colorbars
            #http://stackoverflow.com/questions/23964856/plotting-streamlines-with-matplotlib-python
            # nz = mcolors.Normalize()
            # nz.autoscale(results2[:,4]) 
            # fig = plt.figure("img " + str(initial_start_index)+" - img "+str(stop_index))
            # ax = fig.add_subplot(111)
            # plt.imshow(crop_img1, cmap=plt.cm.gray, origin='upper')
            # plt.title("img " + str(initial_start_index)+" - img "+str(stop_index))
            # #ax.set_color_cycle(['red', 'black', 'yellow'])
            # # Same scale on the x and y axes
            # plt.gca().set_aspect('equal', adjustable='box')

            # plt.ylabel('pixels')
            # plt.xlabel('pixels')

            # # Plot quiver
            # plt.quiver(results2[:,0], results2[:,1], results2[:,2]/results2[:,4], results2[:,3]/results2[:,4],angles='xy', scale=30,color=cm.jet(nz(results2[:,4])))

            # # Colorbar
            # cax,_ = mcolorbar.make_axes(plt.gca())

            # soglia_sup_prova_mm = np.nanmax(results2[:,4])
            # soglia_inf_prova_mm = np.nanmin(results2[:,4])

            # # vmin and vmax should be symmetric? ex: - 6 ,6
            # cb = mcolorbar.ColorbarBase(cax, cmap=cm.jet, norm=mcolors.Normalize(vmin= soglia_inf_prova_mm, vmax= soglia_sup_prova_mm))
            # cb.set_clim(soglia_inf_prova_mm, soglia_sup_prova_mm)

            # #cb = mcolorbar.ColorbarBase(cax, cmap=cm.jet, norm=nz)
            # #cb = mcolorbar.ColorbarBase(cax, cmap=cm.jet)#, norm=mcolors.Normalize(vmin=soglia_inf, vmax=soglia_sup_mm))
            # #cb.set_clim(soglia_inf, soglia_sup_mm)# it doesn't work
            # cb.set_label('mm')
            # plt.savefig("GIFfrec/"+test_name+"_freccette_img_" + str(initial_start_index)+"_img_"+str(stop_index)+".png")
            # mng = plt.get_current_fig_manager()
            # #mng.window.showMaximized()


    ########### SECTION PLOT ########### 

    # plt.figure("Sample central section")
    # plt.plot(dy[:,np.shape(dy)[1]/2.0])
    # plt.ylabel('y displacements (mm)')
    # plt.xlabel('y axis on the grid nodes (mm)')
    # plt.savefig("OutputPlots/"+test_name+"_sezione_img_" + str(initial_start_index)+"_img_"+str(stop_index)+".png")
    plt.show()
    
    stop = timeit.default_timer()
    print ('Time', stop - start) 

    return msg




########### GIF FUNCTION ###########

# def gif(path1,filename):
#     #import imageio as io
#     imgs = []
#     #images = []
#     img_names = glob.glob(path1 + "/*.png")
#     img_names.sort()
#     for img_name in img_names:
#            img = cv2.imread(img_name) 
#            imgs.append(cv2.cvtColor(img, cv2.COLOR_BGR2RGB))
#            #images.append(io.imread(img_name))

#     #print writeGif.__doc__
#     #io.mimsave(path1+'surface1.gif', images, duration = 1)
#     # .png not .jpeg
#     writeGif(path1+filename, imgs, duration=1)



ciao =5






