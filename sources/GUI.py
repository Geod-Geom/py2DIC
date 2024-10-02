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

import sys
import numpy as np
import os
from PyQt5.QtCore import Qt
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
import Main as DIC_main
from matplotlib import pyplot as plt

path = os.path.dirname(os.path.abspath(__file__))+'/'

# If not already present, create the directory to store the results
results_directory = 'OutputPlots'
if not os.path.exists(path+'/'+results_directory):
    os.makedirs(path+'/'+results_directory)

window_width = 500
window_height = 900

rectangley1 = None
rectangley2 = None 
rectanglex1 = None
rectanglex2 = None
absolute_path_of_images = None


class ImageDrawPanel(QGraphicsPixmapItem):
    
    def __init__(self, pixmap=None, parent=None, scene=None, width_scale_ratio=1, height_scale_ratio=1):
        super(ImageDrawPanel, self).__init__()
        self.cont = 0
        self.x, self.y = -1, -1
        self.radius = 10

        self.pen = QPen(Qt.SolidLine)
        self.pen.setColor(Qt.black)
        self.pen.setWidth(2)

        self.brush = QBrush(Qt.yellow)
        self.width_scale_ratio = width_scale_ratio
        self.height_scale_ratio = height_scale_ratio

        self.rect = np.zeros((2,2), dtype=int)

    def paint(self, painter, option, widget=None):
        global rectangley1
        global rectangley2
        global rectanglex1
        global rectanglex2
        rectangley1 = None
        rectangley2 = None 
        rectanglex1 = None
        rectanglex2 = None
        painter.drawPixmap(0, 0, self.pixmap())
        painter.setPen(self.pen)
        painter.setBrush(self.brush)

        if self.x >= 0 and self.y >= 0  and self.x < window_width and self.y < window_height:
            painter.drawEllipse(self.x-self.radius, self.y-self.radius, 2*self.radius, 2*self.radius)
            print(self.cont, self.x,  self.y)
            self.rect[self.cont, 0] = int(self.x)
            self.rect[self.cont, 1] = int(self.y)
           
            self.x, self.y = -1, -1
            self.cont = self.cont+1
        if self.cont ==2:
            print(self.rect)
            painter.drawRect(self.rect[0, 0], self.rect[0, 1], self.rect[1, 0]-self.rect[0, 0], self.rect[1, 1]-self.rect[0, 1])
            self.cont = 0 

            rectangley1 = int(self.rect[0, 1]/self.height_scale_ratio)
            rectangley2 = int(self.rect[1, 1]/self.height_scale_ratio)
            rectanglex1 = int(self.rect[0, 0]/ self.width_scale_ratio)
            rectanglex2 = int(self.rect[1, 0]/ self.width_scale_ratio)

    def mousePressEvent (self, event):
        print ('mouse pressed')
        self.x=int(event.pos().x())
        self.y=int(event.pos().y())
        self.update()

    def mouseMoveEvent (self, event):
        print ('mouse moving')
        self.x = int(event.pos().x())
        self.y = int(event.pos().y())
        self.update()

class Second(QMainWindow):
    def __init__(self, parent=None):
        super(Second, self).__init__(parent)
        
        self.scene = QGraphicsScene()
        self.scene.setSceneRect(0, 0, window_width, window_height)

        pixmap=self.openImage() 
        self.imagePanel = ImageDrawPanel(scene = self.scene)
        self.imagePanel.setPixmap(pixmap)
        self.scene.addItem(self.imagePanel)

        original_width = pixmap.width()
        original_height= pixmap.height()

        self.width_scale_ratio = float(window_width) / original_width
        self.height_scale_ratio= float(window_height) / original_height
        print ('**** py2DIC by Geodesy and Geomatics Division ****')
        # issues with jpeg format at least on Windows
        # solved adding the path to imageformats  app.addLibraryPath('/path/to/plugins/imageformats') 
        #print ('Image original width', original_width)#, '\nwidth ratio',self.width_scale_ratio)
        #print ('Image original height', original_height)#,'\nheight ratio',self.height_scale_ratio)

        pixmap = pixmap.scaled(window_width, window_height, aspectRatioMode=Qt.KeepAspectRatio)
        self.imagePanel = ImageDrawPanel(scene = self.scene, width_scale_ratio=self.width_scale_ratio, height_scale_ratio=self.height_scale_ratio)
        self.imagePanel.setPixmap(pixmap)
        self.scene.addItem(self.imagePanel)

        self.view = QGraphicsView(self.scene)

        layout = QHBoxLayout()
        layout.addWidget(self.view)

        self.widget = QWidget()
        self.widget.setLayout(layout)

        self.setCentralWidget(self.widget)
        self.setWindowTitle("Draw the Area of Interest (AOI)")

    def openImage(self):
        global absolute_path_of_images
        fname = QFileDialog.getOpenFileName(self, "Open image", ".", "Image Files (*.bmp *.JPG *.png *.xpm *.TIF)")
        absolute_path_of_images =  str(fname[0][:-1*len(os.path.basename(str(fname[0])))])

        if len(absolute_path_of_images)==0:
            return None
        return QPixmap(fname[0])


class First(QDialog):
    def __init__(self, parent=None):
        super(First, self).__init__(parent)
           
        global path
        #global absolute_path_of_images
        self.setWindowTitle("Py2DIC")
        self.setGeometry(600, 50, 620, 940)
        # Setting the scroll area
        self.scrollArea = QScrollArea(self)
        self.scrollArea.setWidgetResizable(True)
        self.scrollAreaWidgetContents = QWidget(self.scrollArea)
        self.scrollArea.setWidget(self.scrollAreaWidgetContents)

        #app_icon = QIcon()
        #app_icon.addFile(path+'../Kuaternion.ico', QSize(600, 600))
        #self.setWindowIcon(app_icon)

        # The second layout continer, which contains all the other widget
        self.verticalLayoutScroll = QVBoxLayout(self.scrollAreaWidgetContents)

        # The main layout container, which contains the scroll area
        verticalLayout = QVBoxLayout(self)
        verticalLayout.addWidget(self.scrollArea)

        # Setting all the widgets
        label_image = QLabel(self.scrollAreaWidgetContents)
        pixmap = QPixmap(path+'../logo_sapienza.jpg')

        if pixmap.width() == 0:
            path = path[:-1]+'/'
            pixmap = QPixmap(path+'../logo_sapienza.jpg')

        #pixmap = pixmap.scaled(450, 120)
        label_image.setPixmap(pixmap)
        self.verticalLayoutScroll.addWidget(label_image, alignment=Qt.AlignCenter)

        #self.AOIbutton = QPushButton("Select AOI")
        self.l1 = QLabel(self.scrollAreaWidgetContents)
        self.l1.setText("Pixel dimension [mm]")
        self.le_pixel_dimension = QLineEdit(self.scrollAreaWidgetContents)
        self.le_pixel_dimension.setObjectName("DimPi")
        self.le_pixel_dimension.setText("1")

        #self.l_frame_rate = QLabel(self.scrollAreaWidgetContents)
        #self.l_frame_rate.setText("Camera acquisition time [s]")
        #self.le_frame_rate= QLineEdit(self.scrollAreaWidgetContents)
        #self.le_frame_rate.setObjectName("time")
        #self.le_frame_rate.setText("5")

        #self.labformat= QLabel(self.scrollAreaWidgetContents)
        #self.labformat.setText("Image format (jpg, png...) - case sensitive -")
        #self.lformat = QLineEdit(self.scrollAreaWidgetContents)
        #self.lformat.setObjectName("Imformat")
        #self.lformat.setText("png")
        myFont=QFont()
        myFont.setPointSize(9)
        myFont.setBold(True)
        self.settings= QLabel(self.scrollAreaWidgetContents)
        self.settings.setText("General settings")
        self.settings.setFont(myFont)

        #self.l_prova= QLabel(self.scrollAreaWidgetContents)
        #self.l_prova.setText("Test path")
        self.path_define = QPushButton("Select reference image in folder")
        self.le_prova = QLineEdit(self.scrollAreaWidgetContents)
        self.le_prova.setObjectName("prova")
        self.le_prova.setText(absolute_path_of_images)

        #self.l2 = QLabel(self.scrollAreaWidgetContents)
        #self.l2.setText("Imposed deformation velocity [mm/m]")
        #self.le2 = QLineEdit(self.scrollAreaWidgetContents)
        #self.le2.setObjectName("vel def macchina")
        #self.le2.setText("0.5")

        self.lsi = QLabel(self.scrollAreaWidgetContents)
        self.lsi.setText("Start index")
        self.lesi = QLineEdit(self.scrollAreaWidgetContents)
        self.lesi.setObjectName("start_index")
        self.lesi.setText("0")

        self.lstopi = QLabel(self.scrollAreaWidgetContents)
        self.lstopi.setText("Number of images")
        self.lestopi = QLineEdit(self.scrollAreaWidgetContents)
        self.lestopi.setObjectName("levels")
        self.lestopi.setText("1")

        self.lsamp = QLabel(self.scrollAreaWidgetContents)
        self.lsamp.setText("Image time sampling")
        self.lesamp = QLineEdit(self.scrollAreaWidgetContents)
        self.lesamp.setObjectName("image_time_sampling")
        self.lesamp.setText("1")        

        self.pb = QPushButton(self.scrollAreaWidgetContents)
        self.pb.setObjectName("Run")
        self.pb.setText("Run") 

        self.settings2= QLabel(self.scrollAreaWidgetContents)
        self.settings2.setText("Displacement computation - template matching")
        self.settings2.setFont(myFont)       

        self.l_tem_width = QLabel(self.scrollAreaWidgetContents)
        self.l_tem_width.setText("Template width [pixel]")
        self.le_tem_width = QLineEdit(self.scrollAreaWidgetContents)
        self.le_tem_width.setObjectName("templateWidth")
        self.le_tem_width.setText("11")

        self.l_b = QLabel(self.scrollAreaWidgetContents)
        self.l_b.setText("Edge y [pixel]")
        self.le_b = QLineEdit(self.scrollAreaWidgetContents)
        self.le_b.setObjectName("bordoy")
        self.le_b.setText("15")

        self.l_b1 = QLabel(self.scrollAreaWidgetContents)
        self.l_b1.setText("Edge x [pixel]")
        self.le_b1 = QLineEdit(self.scrollAreaWidgetContents)
        self.le_b1.setObjectName("bordox")
        self.le_b1.setText("2")
        

        self.ml1 = QLabel(self.scrollAreaWidgetContents)
        self.ml1.setText("Multilook [pixel]")
        self.ml2 = QLineEdit(self.scrollAreaWidgetContents)
        self.ml2.setObjectName("Multilook y")
        self.ml2.setText("10")

        self.display = QTextBrowser(self.scrollAreaWidgetContents)
        self.display.verticalScrollBar().setValue(0)
        self.display.verticalScrollBar().maximum()

        self.settings3= QLabel(self.scrollAreaWidgetContents)
        self.settings3.setText("Strain computation")
        self.settings3.setFont(myFont) 

        # Widget Radiobutton gif
        self.defor=QLabel(self.scrollAreaWidgetContents)
        self.defor.setText("Compute strain field")
        self.rdbUno=QCheckBox("",self)
        self.rdbUno.setChecked(True)

        # Adding all the widgets to the scroll layout container
        self.verticalLayoutScroll.addWidget(label_image)

        self.verticalLayoutScroll.addWidget(self.settings)

        #self.verticalLayoutScroll.addWidget(self.l_prova)
        self.verticalLayoutScroll.addWidget(self.path_define)
        self.verticalLayoutScroll.addWidget(self.le_prova)
        #self.verticalLayoutScroll.addWidget(self.AOIbutton)
        #self.verticalLayoutScroll.addWidget(self.labformat)
        #self.verticalLayoutScroll.addWidget(self.lformat)
        self.verticalLayoutScroll.addWidget(self.l1)
        self.verticalLayoutScroll.addWidget(self.le_pixel_dimension)
        #self.verticalLayoutScroll.addWidget(self.l2)
        #self.verticalLayoutScroll.addWidget(self.le2)
        #self.verticalLayoutScroll.addWidget(self.l_frame_rate)
        #self.verticalLayoutScroll.addWidget(self.le_frame_rate)
        self.verticalLayoutScroll.addWidget(self.lsi)
        self.verticalLayoutScroll.addWidget(self.lesi)
        self.verticalLayoutScroll.addWidget(self.lstopi)
        self.verticalLayoutScroll.addWidget(self.lestopi)
        self.verticalLayoutScroll.addWidget(self.lsamp)
        self.verticalLayoutScroll.addWidget(self.lesamp)
        self.verticalLayoutScroll.addWidget(self.settings2)
        self.verticalLayoutScroll.addWidget(self.l_tem_width)
        self.verticalLayoutScroll.addWidget(self.le_tem_width)
        self.verticalLayoutScroll.addWidget(self.l_b)
        self.verticalLayoutScroll.addWidget(self.le_b)
        self.verticalLayoutScroll.addWidget(self.l_b1)
        self.verticalLayoutScroll.addWidget(self.le_b1)
        self.verticalLayoutScroll.addWidget(self.ml1)
        self.verticalLayoutScroll.addWidget(self.ml2)      
        #self.verticalLayoutScroll.addWidget(self.l_V)
        #self.verticalLayoutScroll.addWidget(self.le_V)
        #self.verticalLayoutScroll.addWidget(self.l_H)
        #self.verticalLayoutScroll.addWidget(self.le_H)

        #Widget radiobutton gif
        self.verticalLayoutScroll.addWidget(self.settings3)
        self.verticalLayoutScroll.addWidget(self.defor)
        self.verticalLayoutScroll.addWidget(self.rdbUno)
        self.verticalLayoutScroll.addWidget(self.pb)
        self.verticalLayoutScroll.addWidget(self.display)



        #self.AOIbutton.clicked.connect(self.on_pushButton_clicked)
        self.pb.clicked.connect(self.button_click)
        self.path_define.clicked.connect(self.open1)
        self.le_prova.setText(absolute_path_of_images)

    def open1(self):
        self.dialog = Second(self)
        self.le_prova.setText(absolute_path_of_images)
        return absolute_path_of_images

    def on_pushButton_clicked(self):
        self.dialog = Second(self)
        self.dialog.show()

    def appExit(self):
        # Replace the next line with something that calls the QApplication's
        #   exit() or quit() function.
        #sys.exit()
        app.quit()

    # Button clicked event handler
    def button_click(self):

        try:
            self.dim_pixel = float (self.le_pixel_dimension.text())
            #img_format = str( self.lformat.text())
            #vel = float (self.le2.text())
            start_index = int(self.lesi.text())
            levels =  int (self.lestopi.text())
            image_time_sampling =  int (self.lesamp.text())
            templateWidth = int(self.le_tem_width.text()) 
            b = int(self.le_b.text())
            d = int(self.le_b1.text())
            H = 2*d + templateWidth # numero di griglie che si ripetono sulla singola riga della griglia grande
            V = 2*b + templateWidth
            ml = int (self.ml2.text())
            #frame_rate = 1.0/ float(self.le_frame_rate.text()) 
            prova = str(self.le_prova.text())

            self.display.setText(    self.display.toPlainText() +
                                    "SELECTED PARAMETERS:"+
                                    "\nTest name = "+prova+
                                    #"\nImage Format = "+img_format+
                                    "\nPixel Dimension = "+str(self.dim_pixel)+ ' mm'
                                    #'\nImposed deformation velocity = '+str(vel)+' mm/s'+
                                    #'\nCamera acquisition time = '+str(1.0/frame_rate)+' s'
                                    '\nStart index = ' + str(start_index)+ 
                                    '\nLevels= '  + str(levels)+ 
                                    '\nImage time sampling=' + str(image_time_sampling)+
                                    '\nTemplate width = '+ str(templateWidth)+ ' pixel'+
                                    '\nSearch zone width = ' + str(b)+ ' pixel'+
                                    '\nNumber of grids along x = ' + str(H)+
                                    '\nNumber of grids along y = ' + str(V)+
                                    '\nMultilook = ' + str(ml)
                                    )

        except Exception as e:
            print ("Error:", e)
            self.display.setText(self.display.toPlainText()+"Errore: gli input devono essere un numero\n")
            return False
        try:  

            if self.rdbUno.isChecked():
                defor = True
            else:
                defor = False

            log_message = DIC_main.DIC(prova, self.dim_pixel, start_index, levels, image_time_sampling, templateWidth, b, d, rectangley1, rectangley2, rectanglex1, rectanglex2, H, V, defor, ml)

            self.display.setText(self.display.toPlainText()+ "\n"+ log_message)
            plt.ioff()
            return True
        except:
            print ("Error during cross correlation computation")
            self.display.setText(self.display.toPlainText()+"Error during cross correlation computation\n")

            print ("Unexpected error:", sys.exc_info()[0])
            raise
            return False


if __name__ == '__main__':
    app = QApplication(sys.argv)
    main = First()
    app.setActiveWindow(main)
    main.show()
    
    sys.exit(app.exec_())