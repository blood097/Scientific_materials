# An Automated Instance Segmentation Pipeline for Assessing Nitrogen Fixation Zone Heterogeneity Using YOLO26 Model
<br>`выходные данные статьи`
<br><ссылка>
<p align="center">
  <img width="1200" height="900" src="https://github.com/blood097/Scientific_materials/blob/69295b004c4ef9165bd85867d81578cd587ad3fe/YOLO-learning/Fig-Abstract.png">
</p>

**Data**
<br>`Inf_cell_YOLOv26L-seg_vers2__learning_results.csv` - trained instance segmentation model performance
<br>`Symbiotic-nodule-occupancy.v1-nodule_micro_v1.yolo26.zip` - dataset and annotations in the YOLO26 format
<br>`data_2624.csv`, `data_3841.csv`, `data_AMPS22.csv` - raw data

**Analysis**
<br>`YOLO_nodule_infection.ipynb` - Jupyter Notebook script to detect infected cells on unannotated microphotographs
<br>`Cell-coverage comparison.ipynb` - Jupyter Notebook script to implement the pairwise two-sided permutation tests and Cohen's d, bootstrap confidence intervals for the mean difference
<br>`YOLOv26_cell.ipynb` - training environment on Google Colab

___
Created by Kirichek Evgenii - e.kirichek@arriam.ru

