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
 Valeria Belloni    <valeria.belloni01@gmail.com>
 and is subject to the terms and conditions of the
 GNU Lesser General Public License Version 2.1
 The license text is available from
 http://www.gnu.org/licenses/lgpl.html
 
 More information in the following scientific paper:

 R. Ravanelli, A. Nascetti, M. Di Rita, V. Belloni, D. Mattei, N. Nisticò, M. Crespi. (2017 - in press),
 A new Digital Image Correlation software for displacements field measurement in structural applications.
 ISPRS - International Archives of the Photogrammetry, Remote Sensing and Spatial Information Sciences.

'''

from PyQt4.QtCore import Qt
from PyQt4.QtGui import QBrush, QHBoxLayout,QFileDialog, QMainWindow, QGraphicsView, QGraphicsScene, QPen, QGraphicsPixmapItem, QWidget, QApplication, QPixmap
import sys
import numpy as np
import os
import pdb

window_width = 900
window_height = 500

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
        
        self.rect = np.zeros((2,2))

    def paint(self, painter, option, widget=None):               
        global rectangley1
        global rectangley2
        global rectanglex1
        global rectanglex2
        painter.drawPixmap(0, 0, self.pixmap())                
        painter.setPen(self.pen)
        painter.setBrush(self.brush)
        
        
        if self.x >= 0 and self.y >= 0  and self.x < window_width and self.y < window_height:
            painter.drawEllipse(self.x-self.radius, self.y-self.radius, 2*self.radius, 2*self.radius)
            print self.cont, self.x,  self.y
            self.rect[self.cont, 0] = self.x
            self.rect[self.cont, 1] = self.y
           
            self.x, self.y = -1, -1
            self.cont = self.cont+1
        if self.cont ==2:
			print self.rect
			painter.drawRect(self.rect[0, 0], self.rect[0, 1], self.rect[1, 0]-self.rect[0, 0], self.rect[1, 1]-self.rect[0, 1])
			self.cont = 0 
    
			rectangley1 = int(self.rect[0, 1]/self.height_scale_ratio)
			rectangley2 = int(self.rect[1, 1]/self.height_scale_ratio)
			rectanglex1 = int(self.rect[0, 0]/ self.width_scale_ratio)
			rectanglex2 = int(self.rect[1, 0]/ self.width_scale_ratio)
          
    def mousePressEvent (self, event):

        print 'mouse pressed'
        self.x=event.pos().x()
        self.y=event.pos().y()            
        self.update()

    def mouseMoveEvent (self, event):
        print 'mouse moving'
        self.x = event.pos().x()
        self.y = event.pos().y()            
        self.update()

class MainWindow(QMainWindow):
    def __init__(self):
        super(MainWindow, self).__init__()   

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
        print '**** py2DIC by Geodesy and Geomatics Division ****'
        # issues with jpeg format at least on Windows
		# solved adding the path to imageformats  app.addLibraryPath('/path/to/plugins/imageformats') 
        print 'Image original width', original_width, '\nwidth ratio',self.width_scale_ratio
        print 'Image original height', original_height,'\nheight ratio',self.height_scale_ratio
        
        pixmap = pixmap.scaled(window_width, window_height)
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
		fname = QFileDialog.getOpenFileName(self, "Open image", ".", "Image Files (*.bmp *.JPG *.png *.xpm)")
		absolute_path_of_images = str(fname[:-1*len(os.path.basename(str(fname)))])
		
		if fname.isEmpty(): 
			return None
		return QPixmap(fname)
          
	 
app = QApplication(sys.argv)
path = os.path.dirname(os.path.abspath(__file__))+'/'#
app.addLibraryPath(path)
mainWindow = MainWindow()
mainWindow.resize(window_width,window_height)
mainWindow.show()
app.exec_()
#print 'fine', rectangley1, absolute_path_of_images



