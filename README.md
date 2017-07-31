# cloud-regime-error-metric
This is a set of metrics based on RMS error which indicate the overall performance of a model to simulate the cloud radiative effects (CREs) of cloud regimes in annual mean climatology. The cloud regimes are defined as grid-boxes sharing similar daily mean cloud top pressures, optical depths and cloud covers from the ISCCP data. It is documented in Williams and Webb 2009 (for CFMIP1) and Tsushima et al., 2013 (for CMIP5).

References
----------
Williams KD, Webb MJ (2009) A quantitative performance assessment of cloud regimes in climate models. Clim Dyn 33(1):141–157

Tsushima Y, Ringer MA, Webb MJ, Williams KD (2012) Quantitative evaluation of the seasonal variations in climate model cloud regimes. Clim Dyn 41(9-10):2679-2696

Input
----------

| Frequency | Variable | CMOR labels | Unit | File Format |
|:----------|:-----------------------------|:-------------|:------|:------------|
| daily mean | ISCCP mean cloud top pressure | pctisccp     | Pa    | nc
|  | ISCCP mean cloud albedo | albisccp     |  1    | nc
|  | Total cloud cover simulating ISCCP | cltisccp     |  %    | nc
|  | Surface snow area fraction | snc     |  %    | nc
|  | Sea ice area fraction  | sic     |  1    | nc
|  | Outgoing shortwave flux at the top-of-the-atmosphere(TOA)  | rsut     |  W/m^-2    | nc
|  | TOA outgoing longwave flux  | rlut     |  W/m^-2    | nc
|  | TOA outgoing shortwave flux assuming clear-sky | rsutcs     |  W/m^-2    | nc
|  | TOA outgoing longwave flux assuming clear-sky | rsutcs     |  W/m^-2    | nc
|  | TOA Incident Shortwave Radiation*  | rsdt     |  W/m^-2    | nc

*rsdt is necessary only for CREM for the climatological annual variation in Tsushima et al.,(2013).

Sample input data are provided for CanAM4 during 19790101-19791231 in data dir

The reference properties of observational cloud regimes were estimated using the following data:

Daily ISCCP data for clouds and snow/ice : https://eosweb.larc.nasa.gov/project/isccp/isccp_d1_table

ISCCP-FD Daily radiative flux data: https://isccp.giss.nasa.gov/outgoing/FLUX/TOA/ 

Programs for Diagnostics Calculation and Outputs
----------

### A. Annual Mean Climatology

**Programs**: cloud_regime_error_metric_calc.py (python)

**Output**: Single value texts of Statistical mean of daily mean CREMpd [Wm-2]

Sample input model data files are provided under 'model_data' directory.

For the provided sample input data, the code should print the following output. 

CREMpd:  1.097

**Note**: Please set fixmdis=False, in case the program fails.

**A script to draw a figure in the paper is not included**

### B. Climatological Annual Variation

#### 1. Calculate cloud regime properties for climatological months: 
**Program**: A script 'run_amip' let 'script_for_2dproject_cmip5_mon.pro' run 'fuzzy_project_onto_clusters_2d_ncdf_mon.pro' (idl)

**Input data**: CMIP model data
    
**Auxiliary data**: 

    1.Observational cloud regime property data: global_2d_isccp_nclusters$number_$region.data (e.g. global_2d_isccp_nclusters6_icey.data) are provided in code/idl/auxiliary_data)
    
    2.Reference data with a grid size for model data to regrid to: As an example, 2.5degree grid ISCCP data in Obs4MIPs (cltisccp_obs4MIPs_ISCCP_L3_V1.0_200511-200511.nc) is in code/idl/auxiliary_data
    
**Output**: $region$_onto_isccp.$variable (e.g. 2020_onto_isccp.swcffrac)

**Note**: You need to edit data period in fuzzy_project_onto_clusters_2d_ncdf_mon.pro (e.g. 'from=[1979,1,1],to=[1979,12,31]')

For the provided sample input data, the sample output files are under code/idl)

#### 2. Calculate errors (break down into amplitude error and covariance error): 
**Program**: seasonal_rmserr_breakdown_sqrt.pro (idl)

**Input**: Model data: $region$_onto_isccp.$variable (Obtained from the Calculation 1)

**Auxiliary data**: $region$_onto_isccp.$variable from the observations (provided in code/idl/auxiliary_data)    

**Output**: rmserr_breakdown_sqrt_$variable_$region_onto_isccp

For the provided sample input data, sample output files are under code/idl. 

**A script to draw a figure in the paper is not included**



