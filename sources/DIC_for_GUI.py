# -*- coding: utf-8 -*-
'''

 py2DIC
 2D Digital Image Correlation software
 developed by Geodesy and Geomatics Division   
 University of Rome "La Sapienza"
 
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
import time
import pdb

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

    # Visualization of matching results
    if show == True:
        # Scale images for visualization
        img_master_scaled = cv2.convertScaleAbs(img_master, alpha=(255.0/500))
        img_slave_scaled = cv2.convertScaleAbs(img_slave, alpha=(255.0/500))
        cv2.rectangle(img_slave_scaled,top_left, bottom_right, 255, 2*mlx) 
        plt.figure(figsize=(20,10))
        plt.subplot(131),plt.imshow(res,cmap = 'gray')
        plt.title('Matching Result'), plt.xticks([]), plt.yticks([])
        plt.subplot(132),plt.imshow(img_master_scaled,cmap = 'gray')
        plt.title('Detected Point'), plt.xticks([]), plt.yticks([])
        plt.subplot(133),plt.imshow(img_slave_scaled, cmap = 'gray')
        plt.suptitle(method)
        plt.show()
    
    return px, py, max_val

def template_match_cv_OCL_transp_API(img_master, img_slave, method = 'cv2.TM_CCOEFF_NORMED', mlx = 1, mly = 1, show=True):    
    #https://www.learnopencv.com/opencv-transparent-api/

    # width and height of resampled images
    w = img_master.shape[1] * mlx
    h = img_master.shape[0] * mly

    # Convert the input images in UMat images (GPU) thanks to Open Computing Language (OpenCL)
    img_slave  =  cv2.UMat(img_slave)
    img_master = cv2.UMat(img_master)

    # Apply image oversampling 
    img_master = cv2.resize(img_master,None,fx=mlx, fy=mly, interpolation = cv2.INTER_CUBIC)
    img_slave  = cv2.resize(img_slave, None,fx=mlx, fy=mly, interpolation = cv2.INTER_CUBIC)

    res = cv2.matchTemplate(img_slave, img_master,eval(method))

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

    # Visualization of matching results
    if show == True:
        # Scale images for visualization
        img_master_scaled = cv2.convertScaleAbs(img_master, alpha=(255.0/500))
        img_slave_scaled = cv2.convertScaleAbs(img_slave, alpha=(255.0/500))
        cv2.rectangle(img_slave_scaled,top_left, bottom_right, 255, 2*mlx) 

        # convert back from Umat (GPU) to numpy array for matplotlib
        img_master_scaled = cv2.UMat.get(img_master_scaled)
        img_slave_scaled = cv2.UMat.get(img_slave_scaled)

        plt.figure(figsize=(20,10))
        plt.subplot(131),plt.imshow(res,cmap = 'gray')
        plt.title('Matching Result'), plt.xticks([]), plt.yticks([])
        plt.subplot(132),plt.imshow(img_master_scaled,cmap = 'gray')
        plt.title('Detected Point'), plt.xticks([]), plt.yticks([])
        plt.subplot(133),plt.imshow(img_slave_scaled, cmap = 'gray')
        plt.suptitle(method)
        plt.show()
    
    return px, py, max_val


def DIC(images_absolute_path,format, vel, dim_pixel, frame_rate, start_index, levels, image_time_sampling , temp_dim, b, d, recty1, recty2, rectx1, rectx2, c =0):

    start_time = time.time()
    ml_x = 10
    ml_y = 10

    # threshold values set for the Plate Hole DIC Challenge image collection
    dx_thresh_inf = -1.35
    dx_thresh_sup = 0
    dy_thresh_inf = -12.4
    dy_thresh_sup = -4.5

    #checking that the template dimension is odd (dispari)
    # if not, we make it odd
    if int(temp_dim) % 2 == 0:
        temp_dim = temp_dim + 1

    plt.close("all")
    plt.switch_backend("Qt5Agg")

    # find the last folder, which we assume being the name of the test
    test_name = os.path.basename(os.path.normpath(images_absolute_path))

    initial_start_index = start_index
    print 'dentro DIC',recty1, recty2, rectx1, rectx2, images_absolute_path

    print('Backend: {}'.format(plt.get_backend()))
    print 'Test:', test_name
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

    if recty1!= None and recty2!= None and rectx1!= None and rectx2!= None:
        # otherwise the AOI is all the image
        print 'p2', rectx2, recty2
        print 'p1', rectx1, recty1
        Dim_y = recty2-recty1 
        Dim_x = rectx2-rectx1 

    else:
        Dim_x = np.shape(img1)[1] # width 
        Dim_y = np.shape(img1)[0] # height 

    if np.shape(img1)[1] > np.shape(img1)[0]:
        # In this case the images have a horizontal layout
        # The x axis corresponds to the image columns
        # The y axis corresponds to the image rows
        boh = d
        d = b
        b = boh
        boh = ml_x
        ml_x = ml_y
        ml_y = boh
        
        dy_thresh_inf = -1.35
        dy_thresh_sup = 0        
        dx_thresh_inf    = 4.5#np.min(dx)#- 1.35
        dx_thresh_sup = 12.4#np.max(dx)#0

    print "Dim x:", Dim_x
    print "Dim y:", Dim_y

    # Dimensions of the matching grid
    h = np.int(Dim_y/(temp_dim+2*b+c))# rows
    w = np.int(Dim_x/(temp_dim+2*d+c))# columns
    print 'h grid rows', h
    print 'w grid cols', w

    # Array where to store the results
    # we cumulate the displacements computed for each level
    results = np.zeros((h*w,6))
    results_mm = np.zeros((h*w,6))
    for l in range(levels):

            stop_index = start_index + image_time_sampling
            print stop_index, start_index, image_time_sampling, initial_start_index

            msg = ''
            print "level", l+1
            msg = msg + '-----------'
            msg = msg +  "\nAnalysed images"
            msg = msg +  "\n1: "+ os.path.basename(str(img_names[start_index]))
            msg = msg +  "\n2: "+ os.path.basename(str(img_names[stop_index]))
            print msg

            delta_index = stop_index - start_index
            tempo_passato = delta_index /frame_rate # seconds
            spost_atteso_mm = vel * tempo_passato   # mm
            spost_atteso_pixel = round(spost_atteso_mm/dim_pixel, 0) # pixel

            print '-----------'
            print "Expected displacement"
            print spost_atteso_mm, 'mm ', spost_atteso_pixel, ' pixel'   

            msg = msg + '\n-----------\n'
            msg = msg + "Expected displacement\n"
            msg = msg + str(spost_atteso_mm) + ' mm '+ str(spost_atteso_pixel)+ ' pixel\n' 

            img1 = cv2.imread(img_names[start_index], 0)
            img2 = cv2.imread(img_names[ stop_index], 0)

            if recty1 == None:
                crop_img1 = img1
                crop_img2 = img2
            else:
                crop_img1 = img1[recty1:recty2, rectx1:rectx2]
                crop_img2 = img2[recty1:recty2, rectx1:rectx2]


            print '-----------'
            print "Camera resolution"
            print Dim_x, "x", Dim_y
            print '-----------'
            print "Camera frame rate"
            print frame_rate, "fps (one shot every", 1/frame_rate, "seconds)"

            msg = msg + '\n-----------\n'
            msg = msg + "Camera resolution\n"
            msg = msg + str(Dim_x)+"x"+ str(Dim_y)
            msg = msg + '\n-----------\n'
            msg = msg + "Camera frame rate\n"
            msg = msg + str(frame_rate)+ " fps (one shot every "+str(1/frame_rate)+" seconds)"

            ########### TEMPLATE PARAMETERS ################
            # matchZoneWidth should be greater than the maximum displacement expected
            k=0

            img1r = crop_img1.copy()
            img2r = crop_img2.copy()

            # Cycle before along x and then y to fill in the rows and then the columns 	
            for j in range(h):# loop on rows (Y)
                for i in range(w):# loop on columns (X)

                    Delta_X = i*(c+2*d+temp_dim)
                    Delta_Y = j*(c+2*b+temp_dim)

                    TP_temp_x = Delta_X+c+d+(temp_dim-1)/2.0 # x start
                    TP_temp_y = Delta_Y+c+b+(temp_dim-1)/2.0 # y start

                    start_x_template_slice = c+d+Delta_X
                    stop_x_template_slice  = c+d+Delta_X+temp_dim
                    start_y_template_slice = c+b+Delta_Y
                    stop_y_template_slice  = c+b+Delta_Y+temp_dim

                    shape_template = np.shape( img1r [  start_y_template_slice : stop_y_template_slice , start_x_template_slice : stop_x_template_slice ])

                    #checking the template dimensions
                    assert np.allclose(shape_template[0], temp_dim)
                    assert np.allclose(shape_template[1], temp_dim)
                    assert np.allclose(stop_y_template_slice - start_y_template_slice, temp_dim)
                    assert np.allclose(stop_x_template_slice - start_x_template_slice, temp_dim)
                    assert np.allclose(TP_temp_x, (start_x_template_slice+ stop_x_template_slice)//2.0)
                    assert np.allclose(TP_temp_y, (start_y_template_slice+ stop_y_template_slice)//2.0)

                    start_x_search_slice = c + Delta_X
                    stop_x_search_slice  = c + Delta_X + 2*d + temp_dim
                    start_y_search_slice = c + Delta_Y
                    stop_y_search_slice  = c + Delta_Y + 2*b + temp_dim

                    shape_search = np.shape( img2r [  start_y_search_slice : stop_y_search_slice , start_x_search_slice : stop_x_search_slice ])

                    # checking the search area dimensions
                    assert np.allclose((shape_search[0] - temp_dim) /2.0, b)
                    assert np.allclose((shape_search[1] - temp_dim) /2.0, d)
                    assert np.allclose(TP_temp_x, (start_x_search_slice+ stop_x_search_slice)//2.0)
                    assert np.allclose(TP_temp_y, (start_y_search_slice+ stop_y_search_slice)//2.0)

                    temp =        img1r [  start_y_template_slice : stop_y_template_slice , start_x_template_slice : stop_x_template_slice ]
                    search_area = img2r [  start_y_search_slice : stop_y_search_slice , start_x_search_slice : stop_x_search_slice ]

                    indx,indy, maxcc = template_match(temp.astype('uint8'), search_area.astype('uint8'), mlx = ml_x, mly = ml_y, show = False)
                    #indx,indy, maxcc = template_match_cv_OCL_transp_API(temp, search_area, mlx = ml_x, mly = ml_y, show = False)

                    TP_search_x = Delta_X+c+indx - 0.5   # end point x 
                    TP_search_y = Delta_Y+c+indy - 0.5   # end point y 

                    # Store the results in an array: evaluate if this part  can be removed!
                    '''results[k,0] = TP_temp_x        # start point x [pixel]
                    results[k,1] = TP_temp_y        # start point y [pixel]
                    results[k,2] = TP_search_x-TP_temp_x + results[k,2]  # dx [pixel]
                    results[k,3] = TP_search_y-TP_temp_y + results[k,3]  # dy [pixel]
                    results[k,4] = np.sqrt((results[k,3])**2 + (results[k,2] )**2)  # displ. modulo [pixel]
                    results[k,5] = maxcc+results[k,5]'''

                    # Convert the pixel results into millimetres 
                    results_mm[k,0] = TP_temp_x     # start point x [pixel]
                    results_mm[k,1] = TP_temp_y     # start point y [pixel]
                    results_mm[k,2] = TP_search_x*dim_pixel-TP_temp_x*dim_pixel + results_mm[k,2]  # dx [mm]
                    results_mm[k,3] = TP_search_y*dim_pixel-TP_temp_y*dim_pixel + results_mm[k,3]  # dy [mm]
                    results_mm[k,4] = np.sqrt((results_mm[k,3])**2 + (results_mm[k,2] )**2)        # displ. modulo [mm]
                    results_mm[k,5] = maxcc+results_mm[k,5]

                    k=k+1

            start_index = stop_index

            dx = results_mm[:,2].copy()
            dy = results_mm[:,3].copy()

            dx.shape = (h, w) # the computed displacements have the dimensions of the matching grid
            dy.shape = (h, w)

            # Derivative kernel
            operatore = np.matrix([-1.,0,1.])
            print '-----------'
            print "Used derivative kernel"
            print "dx ", operatore
            print "dy ", operatore.T
            print '-----------'

            msg = msg + '-----------\n'
            msg = msg + "Used derivative kernel\n"
            msg = msg + "dx "+ str(operatore)+ "\n"
            msg = msg + "dy "+ str(operatore.T)
            msg = msg + '\n-----------\n'

            # Convolution for derivative computation
            # Divide by the grid step
            # The grid (dx-dy) comes from a bigger grid 
            # http://docs.scipy.org/doc/scipy-0.16.1/reference/generated/scipy.signal.convolve.html
            Gx = sg.convolve(dx, -operatore/[(2*(2*d+temp_dim+c))*dim_pixel],  mode='same')
            Gy = sg.convolve(dy, -operatore.T/[(2*(2*b+temp_dim+c))*dim_pixel], mode='same')# Multiply by the pixel dimensions to have a millimeters value

            # Printing the files for the single level, before adding the NaNs for visualization purposes
            np.savetxt("OutputPlots/"+test_name+"_results_mm_dx_"+str(initial_start_index)+"_"+str(stop_index)+".txt", dx, fmt = '%.5f')
            np.savetxt("OutputPlots/"+test_name+"_results_mm_dy_"+str(initial_start_index)+"_"+str(stop_index)+".txt", dy, fmt = '%.5f')
            np.savetxt("OutputPlots/"+test_name+"_results_mm_Gy_"+str(initial_start_index)+"_"+str(stop_index)+".txt", Gy, fmt = '%.7f')
            np.savetxt("OutputPlots/"+test_name+"_results_mm_Gx_"+str(initial_start_index)+"_"+str(stop_index)+".txt", Gx, fmt = '%.7f')

            ####  PLOTS 
            # threshold values set for the Plate Hole DIC Challenge image collection
            soglia_inf_y    = dy_thresh_inf#np.min(dy)#-12.4
            soglia_sup_mm_y = dy_thresh_sup#np.max(dy)#- 4.5#np.inf#levels*spost_atteso_pixel*dim_pixel*1.1
            soglia_inf_x    = dx_thresh_inf#np.min(dx)#- 1.35
            soglia_sup_mm_x = dx_thresh_sup#np.max(dx)#0
            soglia_inf_Gx   = - 4/100.0
            soglia_sup_Gx   =   4.5/100.0
            soglia_inf_Gy   = - 4/100.0
            soglia_sup_Gy   =   4/100.0

            ########### DISPLACEMENTS PLOT ########### 
            indici_inf = np.where(dy < soglia_inf_y) # Consider only positive dy
            ix = indici_inf [1] # x index
            iy = indici_inf [0] # y index
            dy[ iy, ix ] = float('NaN')

            indici_sup = np.where(dy > soglia_sup_mm_y) # Consider only positive dy
            ix = indici_sup [1] # x index
            iy = indici_sup [0] # y index
            dy[ iy, ix ] = float('NaN')
            
            indici_inf = np.where(dx < soglia_inf_x) # Consider only positive dy
            ix = indici_inf [1] # x index
            iy = indici_inf [0] # y index
            dx[ iy, ix ] = float('NaN')

            indici_sup = np.where(dx > soglia_sup_mm_x) # Consider only positive dy
            ix = indici_sup [1] # x index
            iy = indici_sup [0] # y index
            dx[ iy, ix ] = float('NaN')

            fig = plt.figure("Displacements dx between img " + str(initial_start_index)+" and img "+str(stop_index))
            plt.title("$Horizontal\;displacements: u$")
            plt.gca().set_aspect('equal', adjustable='box')
            plt.imshow(dx, cmap=cm.jet, norm=mcolors.Normalize(vmin=soglia_inf_x, vmax=soglia_sup_mm_x), interpolation = 'None')
            cb=plt.colorbar() 
            cb.set_clim(soglia_inf_x, soglia_sup_mm_x)
            cb.set_label('$mm$')
            plt.savefig("OutputPlots/"+test_name+"_displacements_dx_img_" + str(initial_start_index)+"_img_"+str(stop_index)+".png")	                

            fig = plt.figure("Displacements dy between img " + str(initial_start_index)+" and img "+str(stop_index)) 
            plt.title("$Vertical\;displacements: v$")
            plt.gca().set_aspect('equal', adjustable='box')
            plt.imshow(dy, cmap=cm.jet, norm=mcolors.Normalize(vmin=soglia_inf_y, vmax=soglia_sup_mm_y), interpolation = 'None')
            cb2=plt.colorbar() 
            cb2.set_clim(soglia_inf_y, soglia_sup_mm_y)
            cb2.set_label('$mm$')
            mng = plt.get_current_fig_manager()
            plt.savefig("OutputPlots/"+test_name+"_displacements_img_" + str(initial_start_index)+"_img_"+str(stop_index)+".png")
            plt.savefig("GIF/"+test_name+"_displacements_dy_img_" + str(initial_start_index)+"_img_"+str(stop_index)+".png")

            ########### STRAIN PLOT ###########
            # Remove the absolute value of strains greater than 4% 
            indici_sup_def = np.where(Gy > soglia_sup_Gy) 
            ix = indici_sup_def [1] # x index
            iy = indici_sup_def [0] # y index
            Gy[ iy, ix ] = float('NaN')
            indici_inf_def = np.where(Gy < soglia_inf_Gy)
            ix = indici_inf_def [1] # x index
            iy = indici_inf_def [0] # y index
            Gy[ iy, ix ] =float('NaN')

            # Remove the absolute value of strains greater than 4% 
            indici_sup_def = np.where(Gx > soglia_sup_Gx) 
            ix = indici_sup_def [1] # x index
            iy = indici_sup_def [0] # y index
            Gx[ iy, ix ] = float('NaN')
            indici_inf_def = np.where(Gx < soglia_inf_Gx)
            ix = indici_inf_def [1] # x index
            iy = indici_inf_def [0] # y index
            Gx[ iy, ix ] =float('NaN')

            fig =  plt.figure("Gx between img " + str(initial_start_index)+" and img "+str(stop_index))
            plt.title("{\partial u \over \partial x}")
            plt.imshow(Gx, cmap=cm.jet, interpolation = 'None')
            plt.colorbar()

            fig = plt.figure("Gy between img " + str(initial_start_index)+" and img "+str(stop_index))
            plt.title("{\partial v \over \partial y}")
            plt.imshow(Gy, cmap=cm.jet, interpolation = 'None')
            plt.colorbar()
            plt.show()


            ########### QUIVER PLOT  ########### 
            print  np.max(results_mm[:,3]), np.min(results_mm[:,3])

            # Remove strains whose absolute value is greater than soglia_sup and lower than soglia_inf
            # limite_inf = np.where(results_mm2[:,3]<soglia_inf)[0]
            # results_mm2 = np.delete(results_mm2, limite_inf, axis=0)
            # limite_sup = np.where(results_mm2[:,3]>soglia_sup_mm)[0]
            # results_mm2 = np.delete(results_mm2, limite_sup, axis=0)

            print  np.nanmax(results_mm[:,3]), np.nanmin(results_mm[:,3]), np.nanmedian(results_mm[:,3]), np.nanmean(results_mm[:,3])

            #http://stackoverflow.com/questions/11970186/matplotlib-quiver-and-imshow-superimposed-how-can-i-set-two-colorbars
            #http://stackoverflow.com/questions/23964856/plotting-streamlines-with-matplotlib-python
            nz = mcolors.Normalize()
            nz.autoscale(results_mm[:,4]) 
            fig = plt.figure("img " + str(initial_start_index)+" - img "+str(stop_index))
            ax = fig.add_subplot(111)
            plt.imshow(crop_img1, cmap=plt.cm.gray, origin='upper')
            plt.title("img " + str(initial_start_index)+" - img "+str(stop_index))
            ax.set_color_cycle(['red', 'black', 'yellow'])
            # Same scale on the x and y axes
            plt.gca().set_aspect('equal', adjustable='box')

            plt.ylabel('pixels')
            plt.xlabel('pixels')

            # Plot quiver
            plt.quiver( results_mm[:,0], #start x  [pixel]
                        results_mm[:,1], #start y  [pixel]
                        results_mm[:,2]/results_mm[:,4], #dx / |displ| [mm/mm]
                        results_mm[:,3]/results_mm[:,4], #dy / |displ| [mm/mm]
                        angles='xy', 
                        scale=30,
                        color=cm.jet(nz(results_mm[:,4])),
                        edgecolor='k', # edge color of the quivers
                        linewidth=.2)

            # Colorbar
            cax,_ = mcolorbar.make_axes(plt.gca())

            soglia_sup_prova_mm = np.nanmax(results_mm[:,4])
            soglia_inf_prova_mm = np.nanmin(results_mm[:,4])

            # vmin and vmax should be symmetric? ex: - 6 ,6
            cb = mcolorbar.ColorbarBase(cax, cmap=cm.jet, norm=mcolors.Normalize(vmin= soglia_inf_prova_mm, vmax= soglia_sup_prova_mm))
            cb.set_clim(soglia_inf_prova_mm, soglia_sup_prova_mm)

            #cb = mcolorbar.ColorbarBase(cax, cmap=cm.jet, norm=nz)
            #cb = mcolorbar.ColorbarBase(cax, cmap=cm.jet)#, norm=mcolors.Normalize(vmin=soglia_inf, vmax=soglia_sup_mm))
            #cb.set_clim(soglia_inf, soglia_sup_mm)# it doesn't work
            cb.set_label('mm')
            plt.savefig("GIFfrec/"+test_name+"_freccette_img_" + str(initial_start_index)+"_img_"+str(stop_index)+".png")
            mng = plt.get_current_fig_manager()
            #mng.window.showMaximized()

    ########### SECTION PLOT ########### 
    plt.figure("Sample central section")
    plt.plot(dy[:,np.shape(dy)[1]/2])
    plt.ylabel('y displacements (mm)')
    plt.xlabel('y axis on the grid nodes (mm)')
    plt.savefig("OutputPlots/"+test_name+"_sezione_img_" + str(initial_start_index)+"_img_"+str(stop_index)+".png")
    print "py2DIC execution time: %s s" % (time.time() - start_time)
    plt.show()

    return msg


########### GIF FUNCTION ###########

def gif(path1,filename):
    from images2gif import writeGif
    #import imageio as io
    imgs = []
    #images = []
    img_names = glob.glob(path1 + "/*.png")
    img_names.sort()
    for img_name in img_names:
           img = cv2.imread(img_name) 
           imgs.append(cv2.cvtColor(img, cv2.COLOR_BGR2RGB))
           #images.append(io.imread(img_name))

    #print writeGif.__doc__
    #io.mimsave(path1+'surface1.gif', images, duration = 1)
    # .png not .jpeg
    writeGif(path1+filename, imgs, duration=1)



ciao =5






