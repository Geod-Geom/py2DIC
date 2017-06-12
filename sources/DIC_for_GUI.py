# -*- coding: utf-8 -*-
"""

 py2DIC
 2D Digital Image Correlation software
 developed by Geodesy and Geomatics Division   
 University of Rome "La Sapienza"
 
 The information in this file is
 Copyright(c) 2017, 
 Andrea Nascetti    <andrea.nascetti@uniroma1.it>,  
 Martina Di Rita    <martina.dirita@uniroma1.it>, 
 Roberta Ravanelli  <roberta.ravanelli@uniroma1.it>
 Valeria Belloni    <valeria.belloni01@gmail.com>
 and is subject to the terms and conditions of the
 GNU Lesser General Public License Version 2.1
 The license text is available from
 http://www.gnu.org/licenses/lgpl.html
 
 More information in the following scientific paper:

 R. Ravanelli, A. Nascetti, M. Di Rita, V. Belloni, D. Mattei, N. Nisticò, M. Crespi. (2017 - in press),
 A new Digital Image Correlation software for displacements field measurement in structural applications.
 ISPRS - International Archives of the Photogrammetry, Remote Sensing and Spatial Information Sciences.

"""

import cv2
import numpy as np
from matplotlib import pyplot as plt
import glob
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import matplotlib.colorbar as mcolorbar
from scipy import signal as sg
import os
from images2gif import writeGif
import prova_mask1 as c

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
    img_slave_scaled = cv2.convertScaleAbs(img_slave, alpha=(255.0/500))
    
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


def DIC(prova,format, vel, dim_pixel, frame_rate, start_index, levels, image_time_sampling , temp_dim, b, d, recty1, recty2, rectx1, rectx2, c =0):
	plt.close("all")
	plt.switch_backend("Qt4Agg")

	initial_start_index = start_index
	print 'dentro DIC',recty1, recty2, rectx1, rectx2
	
	print('Backend: {}'.format(plt.get_backend()))

	plt.rc('text', usetex=True)

	msg = prova + "\n"
	msg = msg + '-----------\n'
	msg = msg + "Imposed deformation velocity\n"
	msg = msg + str(vel)+ " mm/min\n"
	
	vel = vel / 60 # mm/s
	
	format = '.'+format

	img_names = glob.glob("InputImages/" + prova +"/*"+format)
	img_names.sort()
	
	# Read the first image to obtain information about the camera resolution 
	img1 = cv2.imread(img_names[0], 0)
	crop_img1 = img1[recty1:recty2, rectx1:rectx2]
	#img1r = img1.T.copy() # Return the transpose (otherwise the quiver image is horizontal but the images are vertical)
	Dim_x = np.shape(crop_img1)[1] # width 
	Dim_y = np.shape(crop_img1)[0] # height 
	
	# If statement to define width and height in ubuntu
	windows = False
	if Dim_x < Dim_y:
	    print "Windows or MAC"
	    #Dim_x = np.shape(img1)[0] # width 
	    #Dim_y = np.shape(img1)[1] 
	    windows =True
	else:
	     print "UBUNTU"
	     Dim_x = np.shape(crop_img1)[0] # width 
	     Dim_y = np.shape(crop_img1)[1]
	     
	print "Dim x:", Dim_x
	print "Dim y:", Dim_y
	
	# Array where to store the results
	results = np.zeros(((np.int(Dim_x/(temp_dim+d+c))-1)*(np.int(Dim_y/(temp_dim+b+c))-1),6))
        results_mm = np.zeros(((np.int(Dim_x/(temp_dim+d+c))-1)*(np.int(Dim_y/(temp_dim+b+c))-1),6))
	for l in range(levels):
		
			stop_index = start_index + image_time_sampling
			print stop_index, start_index, image_time_sampling, initial_start_index
		
			msg = ''
			print "step", l
			msg = msg + '-----------'
			msg = msg +  "\nAnalysed images"
			msg = msg +  "\n1: "+ img_names[start_index][len ("InputImages/" + prova +"/"):] 
			msg = msg +  "\n2: "+ img_names[stop_index] [len ("InputImages/" + prova +"/"):] 
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
			
			'''cv2.namedWindow("a", cv2.WINDOW_NORMAL)
			cv2.imshow("a", img1)
			cv2.waitKey()'''	 

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

			'''cv2.namedWindow("b", cv2.WINDOW_NORMAL)
			cv2.imshow("b", img1r)
			cv2.waitKey()'''

			# Cycle before along x and then y to fill in the rows and then the columns 	
			for j in range(np.int(Dim_y/(temp_dim+b+c))-1):
				for i in range(np.int(Dim_x/(temp_dim+d+c))-1):
				  
					Delta_X = i*(c+d+temp_dim)
					Delta_Y = j*(c+b+temp_dim)
					
					temp =        img1r[c+b+Delta_Y : c+b+Delta_Y+temp_dim  , c+d+Delta_X :  c+d+Delta_X+temp_dim  ]
					search_area = img2r[c+Delta_Y :   c+Delta_Y+2*b+temp_dim, c+Delta_X   :  c+2*d+Delta_X+temp_dim  ]
					
					indx,indy, maxcc = template_match(temp.astype('uint8'), search_area.astype('uint8'), mlx = 20, mly =2, show = False)
					
					TP_temp_x = Delta_X+c+d+1.0*temp_dim/2.0 # x start
					TP_temp_y = Delta_Y+c+b+1.0*temp_dim/2.0 # y start
					
					TP_search_x = Delta_X+c+indx    # end point x 
					TP_search_y = Delta_Y+c+indy    # end point y 
										
					# Store the results in an array 
					results[k,0] = TP_temp_x        # start point x
					results[k,1] = TP_temp_y        # start point y
					results[k,2] = TP_search_x-TP_temp_x + results[k,2]  # dx
					results[k,3] = TP_search_y-TP_temp_y + results[k,3] # dy
					results[k,4] = np.sqrt((results[k,3])**2 + (results[k,2] )**2)  # modulo
					results[k,5] = maxcc+results[k,5]
					
					# Convert the pixel results into millimetres 
					results_mm[k,0] = TP_temp_x     # start point x pixel
					results_mm[k,1] = TP_temp_y     # start point y pixel
					results_mm[k,2] = TP_search_x*dim_pixel-TP_temp_x*dim_pixel + results_mm[k,2]  # dx mm
					results_mm[k,3] = TP_search_y*dim_pixel-TP_temp_y*dim_pixel + results_mm[k,3] # dy mm
					results_mm[k,4] = np.sqrt((results_mm[k,3])**2 + (results_mm[k,2] )**2)          # absolute value mm
					results_mm[k,5] = maxcc+results_mm[k,5]
					
					
					k=k+1
					 
			start_index = stop_index
			
	                dx = results_mm[:,2].copy()	
	                dy = results_mm[:,3].copy()
	
	                dx.shape = (np.int(Dim_y/(temp_dim+b+c))-1, np.int(Dim_x/(temp_dim+d+c))-1) # the displacements shold have the dimensions of the research grid
	                dy.shape = (np.int(Dim_y/(temp_dim+b+c))-1, np.int(Dim_x/(temp_dim+d+c))-1) 	
	
	                ####  PLOTS
	                soglia_inf = 0
	                soglia_sup_mm = 1.26
                        #soglia_sup_mm = levels*spost_atteso_pixel*dim_pixel*1.1
                        
                        soglia_inf_def = -4/100.0
	                soglia_sup_def = 4/100.0

                        ########### DISPLACEMENTS PLOT ########### 
	                #boh = results_mm[:,4]
                        #boh.shape = (np.int(Dim_y/(temp_dim+b+c))-1, np.int(Dim_x/(temp_dim+b+c))-1)
	
	                indici_inf = np.where(dy < soglia_inf) # Consider only positive dy
	                ix = indici_inf [1] # x index
	                iy = indici_inf [0] # y index
	                dy[ iy, ix ] =float('NaN')
	
	                indici_sup = np.where(dy > soglia_sup_mm) # Consider only positive dy
	                ix = indici_sup [1] # x index
	                iy = indici_sup [0] # y index
	                dy[ iy, ix ] =float('NaN')

    
	                fig = plt.figure("Displacementsdx between img " + str(initial_start_index)+" and img "+str(stop_index))
	                
	                plt.title("$Horizontal\;displacements: u$")
	                plt.gca().set_aspect('equal', adjustable='box')
	                plt.imshow(dx, cmap=cm.jet )

	                cb=plt.colorbar() 
	                cb.set_label('$mm$')
	                
                        plt.savefig("OutputPlots/"+prova+"displacementsdx_img_" + str(initial_start_index)+"_img_"+str(stop_index)+".png")	                
	                
	                fig = plt.figure("Displacementsdy between img " + str(initial_start_index)+" and img "+str(stop_index)) 
	                
	                plt.title("$Vertical\;displacements: v$")
	                plt.gca().set_aspect('equal', adjustable='box')
	                plt.imshow(dy, cmap=cm.jet, norm=mcolors.Normalize(vmin=soglia_inf, vmax=soglia_sup_mm))

	                cb2=plt.colorbar() 
	                cb2.set_clim(soglia_inf, soglia_sup_mm)
	                cb2.set_label('$mm$')
	                mng = plt.get_current_fig_manager()
	                
	                plt.savefig("OutputPlots/"+prova+"displacements_img_" + str(initial_start_index)+"_img_"+str(stop_index)+".png")
	                plt.savefig("GIF/"+prova+"displacementsdy_img_" + str(initial_start_index)+"_img_"+str(stop_index)+".png")
	
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

	                # Convolution 
	                # Divide by the grid step
	                # The grid (dx-dy) comes from a bigger grid 
	                # http://docs.scipy.org/doc/scipy-0.16.1/reference/generated/scipy.signal.convolve.html
	                Gx = sg.convolve(dx, -operatore/[(2*(2*d+temp_dim+c))*dim_pixel],  mode='same')
	                Gy = sg.convolve(dy, -operatore.T/[(2*(2*b+temp_dim+c))*dim_pixel], mode='same')# Multiply by the pixel dimensions to have a millimeters value

	                
	                ########### STRAINS PLOT ###########
	                
	                # Remove the absolute value of strains greater than 4% 
	                indici_sup_def = np.where(Gy > soglia_sup_def) 
	                ix = indici_sup_def [1] # x index
	                iy = indici_sup_def [0] # y index
	                Gy[ iy, ix ] = float('NaN')
	                indici_inf_def = np.where(Gy < soglia_inf_def)
	                ix = indici_inf_def [1] # x index
	                iy = indici_inf_def [0] # y index
	                Gy[ iy, ix ] =float('NaN')
	
	                # Remove the absolute value of strains greater than 4% 
	                indici_sup_def = np.where(Gx > soglia_sup_def) 
	                ix = indici_sup_def [1] # x index
	                iy = indici_sup_def [0] # y index
	                Gx[ iy, ix ] = float('NaN')
	                indici_inf_def = np.where(Gx < soglia_inf_def)
	                ix = indici_inf_def [1] # x index
	                iy = indici_inf_def [0] # y index
	                Gx[ iy, ix ] =float('NaN')
	                
	                fig =  plt.figure("Gx between img " + str(initial_start_index)+" and img "+str(stop_index))
	                
	                plt.title("{\partial u \over \partial x}")
	                plt.imshow(Gx, cmap=cm.jet)
	                plt.colorbar()
	               
	                fig =  plt.figure("Gy between img " + str(initial_start_index)+" and img "+str(stop_index))
	                plt.title("{\partial v \over \partial y}")
	                plt.imshow(Gy, cmap=cm.jet)
	                plt.colorbar()
	                    
	                plt.show()

                        #####  PLOTS  
	                #soglia_inf = 0
	                #soglia_sup = levels*spost_atteso_pixel*1.1 
                        #soglia_sup_mm = levels*spost_atteso_pixel*dim_pixel*1.1
	                #soglia_sup_mm = 0.06
	                
	                ########### QUIVER PLOT  ########### 

                        # Copy the result values in a new array 
                        results2 =  results_mm.copy()
                        print  np.max(results2[:,3]), np.min(results2[:,3])

                        # Remove strains whose absolute value is greater than soglia_sup and lower than soglia_inf
	                # limite_inf = np.where(results_mm2[:,3]<soglia_inf)[0]
	                # results_mm2 = np.delete(results_mm2, limite_inf, axis=0)
	                # limite_sup = np.where(results_mm2[:,3]>soglia_sup_mm)[0]
	                # results_mm2 = np.delete(results_mm2, limite_sup, axis=0)
	  
	                print  np.nanmax(results2[:,3]), np.nanmin(results2[:,3]), np.nanmedian(results2[:,3]), np.nanmean(results2[:,3])
	  
	  
	                #http://stackoverflow.com/questions/11970186/matplotlib-quiver-and-imshow-superimposed-how-can-i-set-two-colorbars
	                #http://stackoverflow.com/questions/23964856/plotting-streamlines-with-matplotlib-python
	                nz = mcolors.Normalize()
	                nz.autoscale(results2[:,4]) 
	                fig = plt.figure("img " + str(initial_start_index)+" - img "+str(stop_index))
	                ax = fig.add_subplot(111)
                        plt.imshow(crop_img1, cmap=plt.cm.gray, origin='upper')
	                plt.title("img " + str(initial_start_index)+" - img "+str(stop_index))
	                ax.set_color_cycle(['red', 'black', 'yellow'])
	                # Same scale on the x and y axes
	                plt.gca().set_aspect('equal', adjustable='box')

	                plt.ylabel('mm')
	                plt.xlabel('mm')
	
	                # Plot quiver
	                plt.quiver(results2[:,0], results2[:,1], results2[:,2]/results2[:,4], results2[:,3]/results2[:,4],angles='xy', scale=30,color=cm.jet(nz(results2[:,4])))

	                # Colorbar
	                cax,_ = mcolorbar.make_axes(plt.gca())
	                	                
	                soglia_sup_prova_mm = np.nanmax(results2[:,4])
	                soglia_inf_prova_mm = np.nanmin(results2[:,4])
	                
	                # vmin and vmax should be symmetric? ex: - 6 ,6
	                cb = mcolorbar.ColorbarBase(cax, cmap=cm.jet, norm=mcolors.Normalize(vmin= soglia_inf_prova_mm, vmax= soglia_sup_prova_mm))
	                cb.set_clim(soglia_inf_prova_mm, soglia_sup_prova_mm)
	                
	                #cb = mcolorbar.ColorbarBase(cax, cmap=cm.jet, norm=nz)
	                #cb = mcolorbar.ColorbarBase(cax, cmap=cm.jet)#, norm=mcolors.Normalize(vmin=soglia_inf, vmax=soglia_sup_mm))
	                #cb.set_clim(soglia_inf, soglia_sup_mm)# it doesn't work
	                cb.set_label('mm')
	                plt.savefig("GIFfrec/"+prova+"_freccette_img_" + str(initial_start_index)+"_img_"+str(stop_index)+".png")
	                mng = plt.get_current_fig_manager()
	                #mng.window.showMaximized()
	                
	                #with open("OutputPlots/"+prova+"_results_mm_dy"+str(initial_start_index)+"_"+str(stop_index)+".txt","w") as file_stats:
		               # for i in xrange(np.shape(dy)[0]):
			              # for j in xrange(np.shape(dy)[1]):
				             # file_stats.write(str(dy[i,j])+'\t')
			              # file_stats.write('\n')
                        with open("OutputPlots/"+prova+"_results_mm_dy"+str(initial_start_index)+"_"+str(stop_index)+".txt","w") as file_stats:
		                for i in xrange(np.shape(dy)[0]):
			                for j in xrange(np.shape(dy)[1]):
				               file_stats.write(str(dy[i,j])+'\t')
			                file_stats.write('\n')
			                
			                
			with open("OutputPlots/"+prova+"_results_mm_Gy"+str(initial_start_index)+"_"+str(stop_index)+".txt","w") as file_stats:
		                for i in xrange(np.shape(Gy)[0]):
			                for j in xrange(np.shape(Gy)[1]):
				               file_stats.write(str(Gy[i,j])+'\t')
	                                file_stats.write('\n')
	                                
	               	with open("OutputPlots/"+prova+"_results_mm_Gx"+str(initial_start_index)+"_"+str(stop_index)+".txt","w") as file_stats:
		                for i in xrange(np.shape(Gx)[0]):
			                for j in xrange(np.shape(Gx)[1]):
				              file_stats.write(str(Gx[i,j])+'\t')
				        file_stats.write('\n') 
			
			with open("OutputPlots/"+prova+"_results_mm_dx"+str(initial_start_index)+"_"+str(stop_index)+".txt","w") as file_stats:
		                for i in xrange(np.shape(dx)[0]):
			                for j in xrange(np.shape(dx)[1]):
				               file_stats.write(str(dx[i,j])+'\t')
			                file_stats.write('\n')	        
				        
        ########## OUTPUT FILES ###########
	
 
	with open("OutputPlots/"+prova+"_results_mm_"+str(initial_start_index)+"_"+str(stop_index)+".txt","w") as file_stats:
		for k in xrange(len(results_mm)):
			file_stats.write(str(results_mm[k,0])+'\t'+str(results_mm[k,1])+'\t'+str(results_mm[k,2])+'\t'+str(results_mm[k,3])+'\t'+str(results_mm[k,4])+'\t'+str(results_mm[k,5])+'\n')
   
	with open("OutputPlots/"+prova+"_results_mm_Gx"+str(initial_start_index)+"_"+str(stop_index)+".txt","w") as file_stats:
		for i in xrange(np.shape(Gx)[0]):
			for j in xrange(np.shape(Gx)[1]):
				file_stats.write(str(Gx[i,j])+'\n')
		
	with open("OutputPlots/"+prova+"_results_mm_Gy"+str(initial_start_index)+"_"+str(stop_index)+".txt","w") as file_stats:
		for i in xrange(np.shape(Gy)[0]):
			for j in xrange(np.shape(Gy)[1]):
				file_stats.write(str(Gy[i,j])+'\t')
	                file_stats.write('\n')
	                
	with open("OutputPlots/"+prova+"_results_mm_dy"+str(initial_start_index)+"_"+str(stop_index)+".txt","w") as file_stats:
		for i in xrange(np.shape(dy)[0]):
			for j in xrange(np.shape(dy)[1]):
				file_stats.write(str(dy[i,j])+'\t')
			file_stats.write('\n')
			
	with open("OutputPlots/"+prova+"_results_mm_dx"+str(initial_start_index)+"_"+str(stop_index)+".txt","w") as file_stats:
		for i in xrange(np.shape(dx)[0]):
			for j in xrange(np.shape(dx)[1]):
				file_stats.write(str(dx[i,j])+'\t')
			file_stats.write('\n')		
			
 
        ########### SECTION PLOT ########### 
        
	plt.figure("Sample central section")
	plt.plot(dy[:,np.shape(dy)[1]/2])
	plt.ylabel('y displacements (mm)')
	plt.xlabel('y axis on the grid nodes (mm)')
	plt.savefig("OutputPlots/"+prova+"_sezione_img_" + str(initial_start_index)+"_img_"+str(stop_index)+".png")
	plt.show()
	return msg

########### GIF ###########

def gif(path1,filename):
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






