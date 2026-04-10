
import sys
import os
import numpy as np
from PyQt5.QtCore import Qt
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib import cm
import Main as DIC_main  # la tua funzione DIC
import cv2

path = os.path.dirname(os.path.abspath(__file__)) + '/'

results_directory = 'OutputPlots'
if not os.path.exists(path + '/' + results_directory):
    os.makedirs(path + '/' + results_directory)

rectangley1 = None
rectangley2 = None
rectanglex1 = None
rectanglex2 = None
absolute_path_of_images = None

class ImageDrawPanel(QGraphicsPixmapItem):
    def __init__(self, width_scale_ratio=1, height_scale_ratio=1):
        super().__init__()
        self.cont = 0
        self.x, self.y = -1, -1
        self.radius = 5
        self.pen = QPen(Qt.red)
        self.pen.setWidth(2)
        color = QColor(255, 255, 150, 80)  
        self.brush = QBrush(color)
        self.width_scale_ratio = width_scale_ratio
        self.height_scale_ratio = height_scale_ratio
        self.rect = np.zeros((2, 2), dtype=int)

    def paint(self, painter, option, widget=None):
        global rectangley1, rectangley2, rectanglex1, rectanglex2
        painter.drawPixmap(0, 0, self.pixmap())
        painter.setPen(self.pen)
        painter.setBrush(self.brush)

        if self.x >= 0 and self.y >= 0:
            painter.drawEllipse(self.x - self.radius, self.y - self.radius,
                                2 * self.radius, 2 * self.radius)
            self.rect[self.cont, 0] = int(self.x)
            self.rect[self.cont, 1] = int(self.y)
            self.x, self.y = -1, -1
            self.cont += 1

        if self.cont == 2:
            painter.drawRect(
                self.rect[0, 0],
                self.rect[0, 1],
                self.rect[1, 0] - self.rect[0, 0],
                self.rect[1, 1] - self.rect[0, 1]
            )
            rectangley1 = int(self.rect[0, 1] * self.height_scale_ratio)
            rectangley2 = int(self.rect[1, 1] * self.height_scale_ratio)
            rectanglex1 = int(self.rect[0, 0] * self.width_scale_ratio)
            rectanglex2 = int(self.rect[1, 0] * self.width_scale_ratio)
            print("AOI:", rectanglex1, rectangley1, rectanglex2, rectangley2)
            self.cont = 0

    def mousePressEvent(self, event):
        scene_pos = self.mapToScene(event.pos())
        item_pos = self.mapFromScene(scene_pos)

        self.x = int(item_pos.x())
        self.y = int(item_pos.y())
        self.update()

class First(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)

        self.setWindowTitle("Py2DIC")
        self.setStyleSheet("""
        QLabel {
            font-size: 15px;
        }

        QCheckBox {
            font-size: 15px;
        }

        QPushButton {
            font-size: 15px;
        }
        """)
        self.setWindowFlags(Qt.Window | Qt.WindowMinimizeButtonHint |
                            Qt.WindowMaximizeButtonHint | Qt.WindowCloseButtonHint)
        self.setGeometry(100, 50, 1600, 900)

        self.colorbar_limits = {}  

        # ===== LEFT SIDE - controlli =====
        self.scrollArea = QScrollArea()
        self.scrollArea.setWidgetResizable(True)
        self.scrollWidget = QWidget()
        self.scrollArea.setWidget(self.scrollWidget)
        self.leftLayout = QVBoxLayout(self.scrollWidget)

        leftTitle = QLabel("Input parameters")
        leftTitle.setAlignment(Qt.AlignCenter)
        leftTitle.setStyleSheet("font-weight: bold; font-size: 16px;")

        # ===== LOGO =====
        self.logo_label = QLabel()
        pixmap_logo = QPixmap(path+'../logo_sapienza.jpg')
        pixmap_logo = pixmap_logo.scaled(300, 200, Qt.KeepAspectRatio, Qt.SmoothTransformation)
        self.logo_label.setPixmap(pixmap_logo)
        self.logo_label.setAlignment(Qt.AlignCenter)

        self.path_define = QPushButton("Select reference image")
        self.le_prova = QLineEdit()
        self.l1 = QLabel("Pixel dimension [mm]")
        self.le_pixel_dimension = QLineEdit("1")
        self.lsi = QLabel("Start index")
        self.lesi = QLineEdit("0")
        self.lstopi = QLabel("Number of images")
        self.lestopi = QLineEdit("1")
        self.lsamp = QLabel("Image time sampling")
        self.lesamp = QLineEdit("1")
        self.l_tem_width = QLabel("Template width [pixel]")
        self.le_tem_width = QLineEdit("11")
        self.l_b = QLabel("Edge y [pixel]")
        self.le_b = QLineEdit("15")
        self.l_b1 = QLabel("Edge x [pixel]")
        self.le_b1 = QLineEdit("2")
        self.ml1 = QLabel("Multilook [pixel]")
        self.ml2 = QLineEdit("10")
        self.rdbUno = QCheckBox("Compute strain field")
        self.rdbUno.setChecked(True)
        self.pb = QPushButton("Run")
        self.pb.setStyleSheet("""
        QPushButton {
            background-color: #2f6fbb;
            color: white;
            font-weight: bold;
            border-radius: 6px;
            padding: 6px;
        }
        QPushButton:hover {
            background-color: #295fa0;
        }
        QPushButton:pressed {
            background-color: #1f4a7a;
        }
        """)
        self.display = QTextBrowser()

        for w in [self.path_define, self.le_prova,
                  self.l1, self.le_pixel_dimension,
                  self.lsi, self.lesi,
                  self.lstopi, self.lestopi,
                  self.lsamp, self.lesamp,
                  self.l_tem_width, self.le_tem_width,
                  self.l_b, self.le_b,
                  self.l_b1, self.le_b1,
                  self.ml1, self.ml2,
                  self.rdbUno, self.pb,
                  self.display]:
            self.leftLayout.addWidget(w)

        self.scrollArea.setSizePolicy(QSizePolicy.Preferred, QSizePolicy.Expanding)

        # ===== CENTER - IMAGE =====
        self.scene = QGraphicsScene()
        self.view = QGraphicsView(self.scene)
        self.view.setRenderHint(QPainter.Antialiasing)
        self.view.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)

        imageTitle = QLabel("Reference image")
        imageTitle.setAlignment(Qt.AlignCenter)
        imageTitle.setStyleSheet("font-weight: bold; font-size: 16px;")

        # ===== RIGHT - PLOT =====
        self.figure = Figure(figsize=(5, 5))
        self.canvas = FigureCanvas(self.figure)
        self.canvas.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)

        self.prev_button = QPushButton("<")
        self.next_button = QPushButton(">")

        for btn in [self.prev_button, self.next_button]:
            btn.setFixedHeight(40)
            btn.setFixedWidth(60)

        self.prev_button.clicked.connect(self.prev_plot)
        self.next_button.clicked.connect(self.next_plot)

        button_layout = QHBoxLayout()
        button_layout.addWidget(self.prev_button)
        button_layout.addWidget(self.next_button)
        self.button_widget = QWidget()

        self.l_cb_min = QLabel("Colorbar min")
        self.le_cb_min = QLineEdit("")
        self.l_cb_max = QLabel("Colorbar max")
        self.le_cb_max = QLineEdit("")

        right_layout = QVBoxLayout()
        right_layout.addWidget(self.canvas)
        right_layout.addLayout(button_layout)
        right_layout.addWidget(self.l_cb_min)
        right_layout.addWidget(self.le_cb_min)
        right_layout.addWidget(self.l_cb_max)
        right_layout.addWidget(self.le_cb_max)

        self.button_widget = QWidget()
        self.button_widget.setLayout(right_layout)

        plotTitle = QLabel("DIC results")
        plotTitle.setAlignment(Qt.AlignCenter)
        plotTitle.setStyleSheet("font-weight: bold; font-size: 16px;")

        # ===== MAIN LAYOUT =====
        mainLayout = QHBoxLayout(self)

        leftWidget = QWidget()
        leftColLayout = QVBoxLayout(leftWidget)
        leftColLayout.setContentsMargins(0, 0, 0, 0)
        leftColLayout.setSpacing(2)
        
        leftColLayout.addWidget(leftTitle, stretch=0)
        leftColLayout.addWidget(self.scrollArea, stretch=1)
        leftColLayout.addWidget(self.logo_label, stretch=0)
        leftWidget.setLayout(leftColLayout)
        leftWidget.setSizePolicy(QSizePolicy.Preferred, QSizePolicy.Expanding)

        mainLayout.addWidget(leftWidget)

        centerWidget = QWidget()
        centerLayout = QVBoxLayout(centerWidget)
        centerLayout.setContentsMargins(0, 0, 0, 0)
        centerLayout.setSpacing(2)
        centerLayout.addWidget(imageTitle, stretch=0)
        centerLayout.addWidget(self.view, stretch=1)
        centerWidget.setLayout(centerLayout)
        centerWidget.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        mainLayout.addWidget(centerWidget)

        rightWidget = QWidget()
        rightLayout = QVBoxLayout(rightWidget)
        rightLayout.setContentsMargins(0, 0, 0, 0)
        rightLayout.setSpacing(2)
        rightLayout.addWidget(plotTitle, stretch=0)
        rightLayout.addWidget(self.canvas, stretch=1)
        rightLayout.addWidget(self.button_widget, stretch=0)
        rightWidget.setLayout(rightLayout)
        rightWidget.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        mainLayout.addWidget(rightWidget)

        mainLayout.setStretch(0, 1)
        mainLayout.setStretch(1, 2)
        mainLayout.setStretch(2, 2)

        # ===== CONNECTIONS =====
        self.path_define.clicked.connect(self.openImage)
        self.pb.clicked.connect(self.button_click)
        self.le_cb_min.editingFinished.connect(self.apply_colorbar_limits)
        self.le_cb_max.editingFinished.connect(self.apply_colorbar_limits)

    def openImage(self):
        fname = QFileDialog.getOpenFileName(self, "Open image", ".", "Images (*.png *.jpg *.bmp *.tif *.tiff)")
        if not fname[0]:
            return
        global absolute_path_of_images
        absolute_path_of_images = os.path.dirname(fname[0]) + "/"
        self.le_prova.setText(absolute_path_of_images)
        pixmap = QPixmap(fname[0])
        self.scene.clear()
        view_w = self.view.width()
        view_h = self.view.height()
        scaled = pixmap.scaled(view_w, view_h, Qt.KeepAspectRatio)
        width_ratio = pixmap.width() / scaled.width()
        height_ratio = pixmap.height() / scaled.height()
        self.imagePanel = ImageDrawPanel(width_ratio, height_ratio)
        self.imagePanel.setPixmap(scaled)
        self.scene.addItem(self.imagePanel)

    def button_click(self):
        try:
            self.dim_pixel = float(self.le_pixel_dimension.text())
            start_index = int(self.lesi.text())
            levels = int(self.lestopi.text())
            image_time_sampling = int(self.lesamp.text())
            templateWidth = int(self.le_tem_width.text())
            self.templateWidth = templateWidth
            b = int(self.le_b.text())
            d = int(self.le_b1.text())
            self.b = b  
            self.d = d
            ml = int(self.ml2.text())
            H = 2 * d + templateWidth
            V = 2 * b + templateWidth
            prova = self.le_prova.text()
            defor = self.rdbUno.isChecked()

            log_message, results_for_plot = DIC_main.DIC(
                prova, self.dim_pixel, start_index, levels,
                image_time_sampling, templateWidth,
                b, d,
                rectangley1, rectangley2,
                rectanglex1, rectanglex2,
                H, V, defor, ml
            )
            self.display.append(log_message)

            self.plot_DIC_results_in_GUI(**results_for_plot, dim_pixel=self.dim_pixel)

        except Exception as e:
            self.display.append(f"Errore: {e}")

    def keyPressEvent(self, event):
        if not hasattr(self, "all_plots") or not self.all_plots:
            return
        if event.key() == Qt.Key_Right:
            self.next_plot()
        elif event.key() == Qt.Key_Left:
            self.prev_plot()

    # ===== PLOTTING =====
    def plot_DIC_results_in_GUI(self, crop_img1, dispx_smoothed, dispy_smoothed,
                                strain_xx2=None, strain_yy2=None, strain_xy2=None,
                                X=None, Y=None, disp_U=None, disp_V=None, mod=None, dim_pixel=1):
        self.dim_pixel = dim_pixel
        self.all_plots = []
        self.all_plots.append(('dx', crop_img1, dispx_smoothed))
        self.all_plots.append(('dy', crop_img1, dispy_smoothed))
        center_col = dispy_smoothed.shape[1] // 2
        section = dispy_smoothed[:, center_col]

        self.all_plots.append(('section_y', section, self.dim_pixel))

        if strain_xx2 is not None:
            self.all_plots.append(('strain_xx', strain_xx2))
            self.all_plots.append(('strain_yy', strain_yy2))
            self.all_plots.append(('strain_xy', strain_xy2))
        if X is not None:
            self.all_plots.append(('quiver', crop_img1, X, Y, disp_U, disp_V, mod, self.dim_pixel))

        self.current_plot_index = 0
        self.display_single_plot(0)

    def display_single_plot(self, idx):
        self.figure.clear()
        plot_data = self.all_plots[idx]
        kind = plot_data[0]

        self.figure.subplots_adjust(left=0.05, right=0.75, top=0.95, bottom=0.05)

        vmin, vmax = self.colorbar_limits.get(idx, (None, None))

        if kind in ('dx', 'dy'):
            img = plot_data[1]
            disp = plot_data[2]
            temp_dim = self.templateWidth
            b = self.b
            d = self.d
            ax = self.figure.add_subplot(111)
            ax.imshow(img[int((temp_dim-1)/2+b)::, int((temp_dim-1)/2+d)::], cmap='gray', origin='upper')
            im_disp = ax.imshow(disp, cmap=cm.jet, alpha=0.5, vmin=vmin, vmax=vmax)
            ax.set_title(f"{'Horizontal' if kind=='dx' else 'Vertical'} displacement")
            ax.axis('off')
            cbar = self.figure.colorbar(im_disp, ax=ax, orientation='vertical', shrink=1.0, pad=0.02)
            unit = "px" if self.dim_pixel == 1 else "mm"
            cbar.set_label(f"Displacement [{unit}]")

        elif kind.startswith('strain'):
            strain = plot_data[1]
            ax = self.figure.add_subplot(111)
            im = ax.imshow(strain, cmap=cm.jet, vmin=vmin, vmax=vmax)
            ax.set_title(f"Strain {kind.split('_')[1].upper()}")
            ax.axis('off')
            cbar = self.figure.colorbar(im, ax=ax, orientation='vertical', shrink=1.0, pad=0.02)
            cbar.set_label(f"Strain [-]")

        elif kind == 'section_y':
            self.figure.subplots_adjust(left=0.15, right=0.85, top=0.95, bottom=0.05)
            section = plot_data[1]
            self.dim_pixel = plot_data[2]
            ax = self.figure.add_subplot(111)
            ax.plot(section, color='blue', linewidth=2)
            ax.set_title("Vertical central section")
            ax.set_xlabel('y position (px)')
            unit = "px" if self.dim_pixel == 1 else "mm"
            ax.set_ylabel(f" y displacement [{unit}]")
            ax.grid(True)

        self.canvas.draw()

    def apply_colorbar_limits(self):
        try:
            vmin = float(self.le_cb_min.text()) if self.le_cb_min.text() else None
            vmax = float(self.le_cb_max.text()) if self.le_cb_max.text() else None
            self.colorbar_limits[self.current_plot_index] = (vmin, vmax)
            self.display_single_plot(self.current_plot_index)
        except ValueError:
            self.display.append("Colorbar limits must be numbers")

    def prev_plot(self):
        if not hasattr(self, "all_plots") or not self.all_plots:
            return
        self.current_plot_index = (self.current_plot_index - 1) % len(self.all_plots)
        self.display_single_plot(self.current_plot_index)

    def next_plot(self):
        if not hasattr(self, "all_plots") or not self.all_plots:
            return
        self.current_plot_index = (self.current_plot_index + 1) % len(self.all_plots)
        self.display_single_plot(self.current_plot_index)

if __name__ == '__main__':
    app = QApplication(sys.argv)
    main = First()
    main.showMaximized()
    sys.exit(app.exec_())