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
|  | Outgoing shortwave flux at the top-of-the-atmosphere(TOA)  | rsut     |  Wm-2    | nc
|  | TOA outgoing longwave flux  | rlut     |  Wm-2    | nc
|  | TOA outgoing shortwave flux assuming clear-sky | rsutcs     |  Wm-2    | nc
|  | TOA outgoing longwave flux assuming clear-sky | rsutcs     |  Wm-2    | nc

Link to the observations (if they are expected in the code):

Monthly ISCCP data: http://climserv.ipsl.polytechnique.fr/cfmip-obs/

Output
----------
Single value texts of Statistical mean of daily mean CREMpd [Wm-2]

 
Is a script to draw a figure in the paper included ?: No
