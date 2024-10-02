# Py2DIC 

<p align="justify"> 
<strong>Py2DIC</strong> is a free and open-source Python software for <strong>2D Digital Image Correlation (DIC)</strong> developed at the Geodesy and Geomatics Division of Sapienza University of Rome.
</p>

<p align="justify"> 
The software compares a series of images of a planar surface collected at different stages of deformation by tracking the pixel movement inside the Area of Interest (AOI) using matching algorithms (template matching method). 
At the end of the processing, py2DIC returns the <strong>displacement</strong> and <strong>strain fields</strong> (Green Lagrangian strains) inside the AOI. For strain computation the software applies smoothing techniques (Gaussian or Spline) to the displacement field to reduce the noise.
</p>

<p align="justify">   
The software allows users to set the main input parameters for displacement and strain computation such as template and search window dimensions (see Figure).
</p>

<p align="center">
  <img width="460" height="300" src="https://github.com/Geod-Geom/py2DIC/blob/master/template_matching2.png">
</p> 

## :wrench: Installation and usage 

To install and test the software, please follow the instructions in the manual. 

## :camera: Datasets

DIC datasets from the Society for Experimental Mechanics (SEM) can be found at https://idics.org/challenge/.

Our DIC datasets are:

- [Dataset 1](https://data.mendeley.com/datasets/dns97tfdjn/1)
- [Dataset 2](https://data.mendeley.com/datasets/z3yc9z84tk/2)

<p align="justify"> 
If you use the datasets in your research, please cite:
  
- *Sjölander, Andreas; Belloni, Valeria; Peterson, Viktor; Ledin, Jonatan* (2023). **Experimental dataset to assess the structural performance of cracked reinforced concrete using Digital Image Correlation techniques with fixed and moving cameras**. Data in Brief, Volume 51, https://doi.org/10.1016/j.dib.2023.109703
  
- *Sjölander, Andreas; Belloni, Valeria; Nascetti, Andrea* (2022), **Dataset to track concrete cracking using DIC with fixed and moving camera**, Mendeley Data, V1, doi: 10.17632/dns97tfdjn.1, https://data.mendeley.com/datasets/dns97tfdjn/1

- *Sjölander, Andreas; Belloni, Valeria; Peterson, Viktor; Ledin, Jonatan* (2023), **Dataset to assess the structural performance of cracked reinforced concrete using FEM, DIC and CMfM**, Mendeley Data, V2, doi: 10.17632/z3yc9z84tk.2, https://data.mendeley.com/datasets/z3yc9z84tk/3
<p>
  
## :pushpin: References

<p align="justify"> 
If you use Py2DIC in your research, please cite the following papers:

- *Belloni V., Ravanelli, R., Nascetti, A., Di Rita, M., Mattei, D., and Crespi, M.*: **py2dic: A new free and open source software for displacement and strain measurements in the field of experimental mechanics**. Sensors, 19(18):3832, https://doi.org/10.3390/s19183832, 2019

- *Belloni V., Ravanelli, R., Nascetti, A., Di Rita, M., Mattei, D., and Crespi, M.*: **Digital Image Correlation from commercial to FOS software: a mature technique for full-field displacement measurements**, The International Archives of the Photogrammetry, Remote Sensing and Spatial Information Sciences, XLII-2, 91-95, https://doi.org/10.5194/isprs-archives-XLII-2-91-2018, 2018

- *Ravanelli R., Nascetti A., Di Rita M., Belloni V., Mattei D., Nisticò N., and Crespi M.*: **A new Digital Image Correlation software for displacements field measurement in structural applications**, The International Archives of the Photogrammetry, Remote Sensing and Spatial Information Sciences, XLII-4/W2, 139-145, https://doi.org/10.5194/isprs-archives-XLII-4-W2-139-2017, 2017


If you want to process images collected with not-fixed cameras refer to our novel methodology Crack Monitoring from Motion (CMfM). Check out the following paper and Github repository:

- *Valeria Belloni and Andreas Sjölander and Roberta Ravanelli and Mattia Crespi and Andrea Nascetti* (2023). **Crack Monitoring from Motion (CMfM): Crack detection and measurement using cameras with non-fixed positions**. In: Automation in Construction, vol 156, https://doi.org/10.1016/j.autcon.2023.105072
  
- https://github.com/Geod-Geom/CMfM
</p>

## License

<p align="justify">
Code is released for non-commercial and research purposes only. For commercial purposes, please contact the authors (valeria.belloni@uniroma1.it).
</p>
