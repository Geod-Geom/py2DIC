# -*- coding: utf-8 -*-
'''

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
 
 More information in the following scientific paper:

 Ravanelli Roberta, Nascetti Andrea, Di Rita Martina, Belloni Valeria, Mattei Domitilla, Nisticò Nicola, and Crespi Mattia: A new Digital
 Image Correlation software for displacements field measurement in structural applications, The International
 Archives of the Photogrammetry, Remote Sensing and Spatial Information Sciences, XLII-4/W2, 139-145,
 https://doi.org/10.5194/isprs-archives-XLII-4-W2-139-2017, 2017

'''

import sys
import os
from PyQt4.QtCore import SIGNAL
from PyQt4.QtGui import *
import DIC_for_GUI as DIC_roby
from matplotlib import pyplot as plt
import subprocess
import prova_mask1 as c
import pdb

path = os.path.dirname(os.path.abspath(__file__))+'/'

# If not already present, create the directory to store the results
results_directory = 'OutputPlots'
if not os.path.exists(path+'\\'+results_directory):
    os.makedirs(path+'\\'+results_directory)

results_directory = 'GIF'
if not os.path.exists(path+'\\'+results_directory):
    os.makedirs(path+'\\'+results_directory)

results_directory = 'GIFfrec'
if not os.path.exists(path+'\\'+results_directory):
    os.makedirs(path+'\\'+results_directory)

class Form(QDialog):

    def __init__(self, parent=None):

        super(Form, self).__init__(parent)

        global path
        self.setWindowTitle("py2DIC by AGG")

        # Setting the scroll area
        self.scrollArea = QScrollArea(self)
        self.scrollArea.setWidgetResizable(True)
        self.scrollAreaWidgetContents = QWidget(self.scrollArea)
        self.scrollArea.setWidget(self.scrollAreaWidgetContents)

        # The second layout continer, which contains all the other widget
        self.verticalLayoutScroll = QVBoxLayout(self.scrollAreaWidgetContents)

        # The main layout container, which contains the scroll area
        verticalLayout = QVBoxLayout(self)
        verticalLayout.addWidget(self.scrollArea)

        # Setting all the widgets
        label_image = QLabel(self.scrollAreaWidgetContents)
        pixmap = QPixmap(path+'..\\logo_sapienza.jpg')

        if pixmap.width() == 0:
            path = path[:-1]+'\\'
            pixmap = QPixmap(path+'..\\logo_sapienza.jpg')

        pixmap = pixmap.scaled(445, 90)
        label_image.setPixmap(pixmap)

        self.l1 = QLabel(self.scrollAreaWidgetContents)
        self.l1.setText("Pixel dimension [mm]")
        self.le_pixel_dimension = QLineEdit(self.scrollAreaWidgetContents)
        self.le_pixel_dimension.setObjectName("DimPi")
        self.le_pixel_dimension.setText("0.023752851")

        self.l_frame_rate = QLabel(self.scrollAreaWidgetContents)
        self.l_frame_rate.setText("Camera acquisition time [s]")
        self.le_frame_rate= QLineEdit(self.scrollAreaWidgetContents)
        self.le_frame_rate.setObjectName("time")
        self.le_frame_rate.setText("5")

        self.labformat= QLabel(self.scrollAreaWidgetContents)
        self.labformat.setText("Image format (jpg, png...) - case sensitive -")
        self.lformat = QLineEdit(self.scrollAreaWidgetContents)
        self.lformat.setObjectName("Imformat")
        self.lformat.setText("JPG")

        self.l_prova= QLabel(self.scrollAreaWidgetContents)
        self.l_prova.setText("Test path")
        self.le_prova = QLineEdit(self.scrollAreaWidgetContents)
        self.le_prova.setObjectName("prova")
        self.le_prova.setText(c.absolute_path_of_images)

        self.l2 = QLabel(self.scrollAreaWidgetContents)
        self.l2.setText("Imposed deformation velocity [mm/m]")
        self.le2 = QLineEdit(self.scrollAreaWidgetContents)
        self.le2.setObjectName("vel def macchina")
        self.le2.setText("0.5")

        self.lsi = QLabel(self.scrollAreaWidgetContents)
        self.lsi.setText("Start index")
        self.lesi = QLineEdit(self.scrollAreaWidgetContents)
        self.lesi.setObjectName("start_index")
        self.lesi.setText("1")

        self.lstopi = QLabel(self.scrollAreaWidgetContents)
        self.lstopi.setText("Levels")
        self.lestopi = QLineEdit(self.scrollAreaWidgetContents)
        self.lestopi.setObjectName("levels")
        self.lestopi.setText("3")

        self.lsamp = QLabel(self.scrollAreaWidgetContents)
        self.lsamp.setText("Image time sampling")
        self.lesamp = QLineEdit(self.scrollAreaWidgetContents)
        self.lesamp.setObjectName("image_time_sampling")
        self.lesamp.setText("4")        

        self.pb = QPushButton(self.scrollAreaWidgetContents)
        self.pb.setObjectName("Run")
        self.pb.setText("Run") 

        self.l_tem_width = QLabel(self.scrollAreaWidgetContents)
        self.l_tem_width.setText("Template width [pixel]")
        self.le_tem_width = QLineEdit(self.scrollAreaWidgetContents)
        self.le_tem_width.setObjectName("templateWidth")
        self.le_tem_width.setText("65")

        self.l_b = QLabel(self.scrollAreaWidgetContents)
        self.l_b.setText("Edge y [pixel]")
        self.le_b = QLineEdit(self.scrollAreaWidgetContents)
        self.le_b.setObjectName("bordoy")
        self.le_b.setText("28")

        self.l_b1 = QLabel(self.scrollAreaWidgetContents)
        self.l_b1.setText("Edge x [pixel]")
        self.le_b1 = QLineEdit(self.scrollAreaWidgetContents)
        self.le_b1.setObjectName("bordox")
        self.le_b1.setText("5")

        self.display = QTextBrowser(self.scrollAreaWidgetContents)
        self.display.verticalScrollBar().setValue(0)
        self.display.verticalScrollBar().maximum()

        # Widget Radiobutton gif
        self.gif=QLabel(self.scrollAreaWidgetContents)
        self.gif.setText("GIF")
        self.rdbUno=QCheckBox("",self)
        self.rdbUno.setChecked(True)

        # Adding all the widgets to the scroll layout container
        self.verticalLayoutScroll.addWidget(label_image)
        self.verticalLayoutScroll.addWidget(self.l_prova)
        self.verticalLayoutScroll.addWidget(self.le_prova)
        self.verticalLayoutScroll.addWidget(self.labformat)
        self.verticalLayoutScroll.addWidget(self.lformat)
        self.verticalLayoutScroll.addWidget(self.l1)
        self.verticalLayoutScroll.addWidget(self.le_pixel_dimension)
        self.verticalLayoutScroll.addWidget(self.l2)
        self.verticalLayoutScroll.addWidget(self.le2)
        self.verticalLayoutScroll.addWidget(self.l_frame_rate)
        self.verticalLayoutScroll.addWidget(self.le_frame_rate)
        self.verticalLayoutScroll.addWidget(self.lsi)
        self.verticalLayoutScroll.addWidget(self.lesi)
        self.verticalLayoutScroll.addWidget(self.lstopi)
        self.verticalLayoutScroll.addWidget(self.lestopi)
        self.verticalLayoutScroll.addWidget(self.lsamp)
        self.verticalLayoutScroll.addWidget(self.lesamp)
        self.verticalLayoutScroll.addWidget(self.l_tem_width)
        self.verticalLayoutScroll.addWidget(self.le_tem_width)
        self.verticalLayoutScroll.addWidget(self.l_b)
        self.verticalLayoutScroll.addWidget(self.le_b)
        self.verticalLayoutScroll.addWidget(self.l_b1)
        self.verticalLayoutScroll.addWidget(self.le_b1)

        #Widget radiobutton gif
        self.verticalLayoutScroll.addWidget(self.gif)
        self.verticalLayoutScroll.addWidget(self.rdbUno)
        self.verticalLayoutScroll.addWidget(self.pb)
        self.verticalLayoutScroll.addWidget(self.display)

        # Button clicked event 
        self.connect(self.pb, SIGNAL("clicked()"),self.button_click)


    # Button clicked event handler
    def button_click(self):
        try:
            self.dim_pixel = float (self.le_pixel_dimension.text())
            img_format = str( self.lformat.text())
            vel = float (self.le2.text())
            start_index = int(self.lesi.text())
            levels =  int (self.lestopi.text())
            image_time_sampling =  int (self.lesamp.text())
            templateWidth = int(self.le_tem_width.text()) 
            b = int(self.le_b.text())
            d =int(self.le_b1.text())
            frame_rate = 1.0/ float(self.le_frame_rate.text()) 
            prova = str(self.le_prova.text())
            recty1 = c.rectangley1
            recty2 = c.rectangley2
            rectx1 = c.rectanglex1
            rectx2 = c.rectanglex2

            print "Selected parameters:"
            print self.dim_pixel, img_format, vel, start_index, levels, templateWidth, image_time_sampling, b, frame_rate
            print 'ciao',recty1, recty2, rectx1, rectx2
            self.display.setText(    self.display.toPlainText() +
                                    "SELECTED PARAMETERS:"+
                                    "\nTest name = "+prova+
                                    "\nImage Format = "+img_format+
                                    "\nPixel Dimension = "+str(self.dim_pixel)+ ' mm'
                                    '\nImposed deformation velocity = '+str(vel)+' mm/s'+
                                    '\nCamera acquisition time = '+str(1.0/frame_rate)+' s'
                                    '\nStart index = ' + str(start_index)+ 
                                    '\nLevels= '  + str(levels)+ 
                                    '\nImage time sampling=' + str(image_time_sampling)+
                                    '\nTemplate width = '+ str(templateWidth)+ ' pixel'+
                                    '\nSearch zone width = ' + str(b)+ ' pixel\n\n')


        except:
            print "Error: input must be numbers"
            self.display.setText(self.display.toPlainText()+"Errore: gli input devono essere un numero\n")
            return False
        try:  
                
            log_message = DIC_roby.DIC(prova,img_format, vel, self.dim_pixel, frame_rate, start_index, levels, image_time_sampling, templateWidth, b, d, recty1, recty2, rectx1, rectx2)

            if self.rdbUno.isChecked():
                DIC_roby.gif(path+"GIF\\", "gif.GIF")
            else:
                print 'no gif'

            self.display.setText(self.display.toPlainText()+ "RESULTS\n"+ log_message)
            plt.ioff()
            return True
        except:
            print "Error during cross correlation computation"
            self.display.setText(self.display.toPlainText()+"Error during cross correlation computation\n")

            print "Unexpected error:", sys.exc_info()[0]
            raise
            return False

app = QApplication(sys.argv)
form = Form()

form.show()
form.resize(500,900)

sys.exit(app.exec_())


