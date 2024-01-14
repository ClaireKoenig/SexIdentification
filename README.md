# Sex identification application

This open-sourced pipeline enable acurate sex identification from dental enamel based on PRM mass spectrometric data. 

# Input
3 custom made reports from Skyline are requiered as an input.
1. Peptide intensities data (Peptide, Protein, Replicate, Peptide Peak Found Ratio, Normalised Area, Library Dot Product, Isotope Dot Product, Peptide Modified Sequence, File Name, Sample Type, Analyte Concentration and Concentration Multiplier)
2. Raw intensities data (Modified Sequence, Replicate, Precursor, Transition, Raw Times, Raw Intensities, Sample Type)
3. Experimental desing (File Name, Experiment) - not from Skyline

# Principle
The data is filtered based on the isotop dot product (≥ 0.9), the library dot product (≥ 0.9) and the peptide peak found ratio (≥ 0.9). The Standards are subseted from the samples and the limit of detection and quantification are calculated based on linear regression. The limit of detection per target is calculated as 3.3 x (Residual Standard Error / Slope) and the limit of quantification is calculated as 10 x (Residual Standard Error / Slope). The targets with measured intensities below their respective limits of quantification are filtered out. AMELX and AMELY intensities are calculated as the sum of the intensities of their targets. Samples measured with a non-zero AMELY intensity are identified as males. The AMELY/AMELX intensity relationship is modelled by (1) generating a linear regression based on the prior male identifications from the dataset studied or by (2)  using a predefined model based on the analysis of modern material generated in this study. For the first option, outliers can be controlled by allowing a maximum percentage of the data points to be removed. If a samples residual is > 2 * standard error of the residuals in the model, that outlier is removed and the model is recalculated. This process is performed as long as there are either no outliers left in the dataset or the maximum number of points allowed has been removed. The confidence interval of the model is calculated with a 95% confidence as ± 2 * standard error of the fitted values. For each sample presenting a null AMELY intensity, the predicted AMELY intensity is calculated based on the model with the corresponding confidence intervals. AMELX and AMELY LOQ and LOD are calculated as the highest LOQ or LOD of the targets corresponding to AMELX or AMELY respectively. A sample is identified as female if the predicted AMELY intensity is above AMELY LOQ and if the lower limit of the predicted AMELY is above 90% of AMELY LOQ and above the AMELY LOD. If the sample does not fit one of these conditions, it will be annotated as “Nonconclusive”.

# Analysis of the data using the Shiny interface

a.	Open the User_Interface_Sex_Identification.R file using R studio. 
b.	Click on “Run App”. 
c.	The Shiny app window opens. 
![image](https://github.com/ClaireKoenig/SexIdentification/assets/134442809/d6cbfad4-3637-4a59-9068-ed4fa0c84a00)

d.	Load the tables and generate the standard curves: 
  1.	Load the “PeptideIntensitiesData” report.
  2.	Load the “ExperimentalDesign” reports.
  3.	Click on “Plot STD curves”. 
  4.	The standard curves per target, per experiment are plotted and the LOD and LOQ are calculated.

<img width="471" alt="image" src="https://github.com/ClaireKoenig/SexIdentification/assets/134442809/8bbef7a1-2846-4e4c-a872-b6bd5c6f1b36">

e.	Get the sex identification: 
  1.	Click on the “Summary” tab.
  2.	Choose the model (Experimental or Pre-defined).
  3.	If the Experimental model is checked, the maximum percentage of data points can be adjusted to remove outliers. By default, the value is set at 0. 
  4.	Click on Plot model.
  5.	On the left, a summary plot with AMELY intensity = f(AMELX intensity) is displayed. On the right, the model is shown. 
<img width="446" alt="image" src="https://github.com/ClaireKoenig/SexIdentification/assets/134442809/04039796-936e-4eb6-a27b-9867245c8d4a">

6.	Click on Plot table. 
7.	The table summarizes the identification per sample.
8.	The table can be downloaded by clicking on “Download the data”. 
<img width="475" alt="image" src="https://github.com/ClaireKoenig/SexIdentification/assets/134442809/652bb3b5-9494-447c-a58f-03032a642f78">

f.	Go back to the MS signal on the Shiny app.
  1.	Click on the “Signal per sample” tab.
  2.	Load the “RawIntensitiesData” report.
  3.	Select one sample.
  4.	Click on “Plot XIC”.
  5.	The App is plotting the extracted ion current of the targets for the selected sample. The sex is displayed in the first panel. On the top row are the AMELX targets and on the bottom raw are the AMELY targets. 

<img width="476" alt="image" src="https://github.com/ClaireKoenig/SexIdentification/assets/134442809/5d9f44c7-bae1-4bdd-b0c6-081d6452eb53">








Automated data analysis for sex identification

I.	Extraction of the ion signal using Skyline:
a.	Download skyline from the website: https://skyline.ms/project/home/software/Skyline/begin.view
b.	Open the “SexIdentificationTemplate.sky” document in Skyline. Once opened, the document should look like that and contain 2 proteins, 5 peptides, 5 precursors and 57 transitions















c.	Check that the defined parameters are preserved. 
1.	Peptide settings (Settings -> Peptide settings)
i.	Digestion: Enzyme -> Unspecific (create it if needed by adding the peptide alphabet (ACDEFGHIKLMNPQRSTVWY) as being cleaved at the C-terminal). Max missed cleavages: 0. Background proteome: None.
ii.	Prediction: (keep default) Retention time predictor: None. 
iii.	Filter: (keep default) Min length: 6. Max length: 15. Exclude N-terminal AAs: 25.
iv.	Library: Check “Sex_ID_targets”. If the library is not present, generate the library using the following files: 20230807_EXPL2_Evo0_CLK_Enamel_Peptides_STDcurve_STD_A_1.raw file and Sex_ID_Targets.ssl. The files can be downloaded from Pride and should be saved in the same folder. (Click on “build”, name the library, select the wanted output path, make sure the selected data source is “Files” and click on “Next”. Click on “Add Files” and select the .ssl file. Click on “Finish”).
v.	Modifications: Make sure oxidation (M) is checked. Allow maximum one variable modification. 
vi.	Quantification: Regression fit: None. Normalization method: None. Regression weighting: None. MS level: 2.
2.	Transition settings (Settings -> Transition Settings):
i.	Prediction: keep default
ii.	Filter: Precursor charges: 2, 3. Ion charges: 1, 2. Ion types: y, b, p. Product ion selection: From “m/z > precursor” to “6 ions”. No special ions need to be checked. Make sure that “auto-select all matching transitions” is checked. 
iii.	Library: Ion march tolerance: 0.5 m/z. 
iv.	Instrument: keep default.
v.	Full-Scan: MS1 filtering: Isotop peaks included: “count”. Precursor mass analyzer “Centroided”. Peaks “3”. Mass Accuracy: “10 ppm”. Isotope labeling enrichment: “Default”. MS/MS filtering: Acquisition method: “PRM”. Product mass analyzer “Centroided”. Mass Accuracy: “10 ppm”. Retention time filtering: “Include all matching scans”.
vi.	Ion mobility: None.
d.	Build the custom-made reports:
1.	Click on File -> Export -> Report -> Edit List -> Add
2.	Name the report in the Report Name. Select the columns as indicated. The order of the columns has to be preserved. The columns can be found faster using CTRL+F.
i.	PeptideIntensitiesData: 

ii.	RawIntensitiesData:











Click on the “Filter” tab. Apply a filter on the Sample Type to remove the standards from the report.









e.	Load the raw files in Skyline: 
1.	File -> Import -> Results. 
2.	Keep the default parameters (add single-injection replicates in files with “none” for the optimizing.) Click “OK”.
3.	Select the raw files to load and click on “open”. Click “Ok” to remove the common pattern from the raw files.


4.	When this window appears, the files are loading. It will take a few minutes depending on the number of files loaded and the performance of the computer. 

f.	Manually check the peak integration. If a peak is not properly integrated, proceed to manual integration by moving the mouse under the x-axis and pressing on the click button. On the left, the chromatogram shows a peak before reintegration and on the right, the chromatogram shows the same peak after reintegration. 

g.	Add the information about the standards: 
1.	Click on View -> Document grid. 
2.	Pick the “PeptideIntensitiesData” in reports. 
3.	In the Sample Type column, pick the right sample type per file: Standards, Blanks or Unknown for a sample. 
4.	In the Analyte Concentration column, add the concentration of the standards as shown below: 
Standard	A	B	C	D	E	F
Analyte Concentration	10	7.5	5	2.5	1	0.1

5.	In the Concentration Multiplier column, add the multiplying factors as shown below: 
Peptide Modified Sequence	SMIRPPY	SM[+16]IRPPY	SM[+16]IRPPYS	SIRPPYPSY	SIRPPYPSYG
Concentration Multiplier	1	0.2	1	1	1

h.	Export the two custom made reports: 
1.	File -> Export -> Report. 
2.	Select PeptideIntensitiedData and click on “Export”. 
3.	Proceed the same way for “RawIntensitiesData”. 
II.	Generate experimental design table: 
a.	Create a table with 2 columns as in the example: 
1.	“File Name”: Copy paste the names of the raw files from the PeptideIntensitiesData report (STD and samples).
2.	“Experiment”: Name your experiment.
3.	Save the table as a txt file.

b.	If there is only one experiment, the experiment name should be the same for all the files. 
c.	If there is more than one experiment, specify different experiment names per experiment. Also provide to the script a unique PeptideIntensitiesData and RawIntensitiesData files. Different reports can be exported from Skyline separately and combined manually afterwards.

III.	Analysis of the data using a Shiny interface
a.	Open the User_Interface_Sex_Identification.R file using R studio. 
b.	Click on “Run App”. 













c.	The Shiny app window opens. 









d.	Load the tables and generate the standard curves: 
1.	Load the “PeptideIntensitiesData” report.
2.	Load the “ExperimentalDesign” reports.
3.	Click on “Plot STD curves”. 
4.	The standard curves per target, per experiment are plotted and the LOD and LOQ are calculated. 




















e.	Get the sex identification: 
1.	Click on the “Summary” tab.
2.	Choose the model (Experimental or Pre-defined).
3.	If the Experimental model is checked, the maximum percentage of data points can be adjusted to remove outliers. By default, the value is set at 0. 
4.	Click on Plot model.
5.	On the left, a summary plot with AMELY intensity = f(AMELX intensity) is displayed. On the right, the model is shown. 











6.	Click on Plot table. 
7.	The table summarizes the identification per sample.
8.	The table can be downloaded by clicking on “Download the data”. 












f.	Go back to the MS signal on the Shiny app.
1.	Click on the “Signal per sample” tab.
2.	Load the “RawIntensitiesData” report.
3.	Select one sample.
4.	Click on “Plot XIC”.
5.	The App is plotting the extracted ion current of the targets for the selected sample. The sex is displayed in the first panel. On the top row are the AMELX targets and on the bottom raw are the AMELY targets. 




