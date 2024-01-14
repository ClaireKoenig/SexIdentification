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




