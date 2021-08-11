# AutoQuake
This tool automatically extracts important data relevant to a specific earthquake.
The user provides the country, start date and end date of the earthquake, or only the USGS ID, and the tool will automatically extract the earthquake hazard shakemap from USGS, then overlap onto (SEDAC) population data and (Kummu) GDP data and generate spatial predictions of the displaced population using the IIDIPUS model. 
# Files
The most important file is 'Main.R', the only file the user really needs in order to generate a (modified) raster object of the earthquake. All that is required is that the 'input' variable be modified. A description of each file is given below:
- Functions.R: Basic functions that are required to run some of the scripts, for example, converting a date ('2021-03-15') into a year (2021).
- GetData.R: Real-time component of the tool to look out for recent earthquakes from the GDACS database
- GetDisaster.R: This script will become more useful when AutoQuake is extended to include a range of other natural hazards.
- GetDisplacements.R: Extract displacement estimate data from IDMC, including historical data (pre-2017) from the GIDD
- GetINFORM.R: Automatically extract information from the JRC-INFORM index
- GetInitialValues.R: Values required for the statistical model applied to predict the displaced population
- GetODDPackages.R: Load the required R modules (e.g. magrittr) and also source all the files that form AutoQuake
- GetOSM.R: Extract information (e.g. buildings) data from OpenStreetMaps
- GetPopDemo.R: Extract spatial population data (count & density), as well as age and gender demographic data from SEDAC. Note that the data has to already have been extracted from SEDAC, this is not automated through an API.
- GetSocioEconomic.R: Extract spatial GDP (Kummu) data and automatically access country-specific income distributions from the World Inequality Database
**- GetUSGS.R: Automated extraction of earthquake shakemap raster data from USGS
- HAZARDobj.R: class structure for a hazard object, used to store each USGS earthquake**
- INSTALL_INSTRUCTIONS.txt: All the current information required to be able to self-install AutoQuake onto a Linux server/computer
**- Main.R: Simple script that unifies all files to allow the user to quickly & simply extract and form a (modified) raster object of a specific earthquake event**
- Method.R: Purpose-built Markov Chain Monte Carlo (MCMC) algorithm used to optimise the IIDIPUS statistical model parameters.
- Model.R: Statistical model used to generate spatial predictions of the displaced population, trained on historical information.
**- ODDobj.R: class structure for an Oxford Disaster Displacement (ODD) object, used to store earthquake and demographic information about the affected area.**

