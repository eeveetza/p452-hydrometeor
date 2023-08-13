# MATLAB/Octave Implementation of Recommendation ITU-R P.452 (Hydrometeor Scatterer)


This code repository contains a MATLAB/Octave software implementation of Recommendation [ITU-R P.452-18](https://www.itu.int/rec/R-REC-P.452/en) with a prediction procedure for the evaluation of transmission loss due to hydrometeor scatterer for frequencies above about 100 MHz.  

This is a development code. It contains two versions:
* `tl_p452_hydrometeor_vec.m` which is based on the equations in a vector form (from Recommendation [ITU-R P.452-18](https://www.itu.int/rec/R-REC-P.452/en)  and the corresponding [Fascicle ITU-R 3M/FAS/12](https://www.itu.int/oth/R0A04000092/en)) and 
* `tl_p452_hydrometeor.m` which is based on the equations in a scalar form from  Recommendation [ITU-R P.452-18](https://www.itu.int/rec/R-REC-P.452/en)


Comparison between the scalar form in `tl_p452_hydrometeor.m` and the vector form in  `tl_p452_hydrometeor_vec.m` is performed in `test_scalar_vs_vector_eqs.m` providing numerically identical values.



## Function Call

~~~
Lpq = tl_p452_hydrometeor(d, f, lat1, lon1, lat2, lon2, h1loc, h2loc, alpha1loc, alpha2loc, eps1loc, eps2loc, G1, G2, BW1, BW2, p, q, press, wat_vap_dens, Rm, T, M1, M2, M3)
~~~


## Required input arguments of function `tl_p452_hydrometeor`

| Variable          | Type   | Units | Limits       | Description  |
|-------------------|--------|-------|--------------|--------------|
| `d`               | scalar double |  km  |        |   Great-circle path distance        |
| `f`               | scalar double | GHz   | ~0.1 ≤ `f` ≤ ~50 | Frequency   |
| `lat1`, `lon1`               | scalar double | deg   |  | Latitude and Longitude for Station 1  |
| `lat2`, `lon2`               | scalar double | deg   |  | Latitude and Longitude for Station 2  |
| `h1loc`, `h2loc`               | scalar double | km  |  | Local heights above mean sea level of Stations 1 and 2 |
| `alpha1loc`, `alpha2loc`               | scalar double | rad |  | Local bearings (azimuth) of Station 1 from Station 2 and Station 2 from Station 1, in the clockwise sense |
| `eps1loc`, `eps2loc`               | scalar double | rad |  | Local horizon angles (elevation for Station 1 and 2) |
| `G1`, `G2`               | scalar double | - |  | Gain for each antenna as a function of both antenna boresight angle and antenna polarization (linear!) |
| `BW1`, `BW2`               | scalar double | rad |  | Antenna beam widths (either main beam or side lobes) depending on the required coupling  |
| `p`, `q`               | scalar int | - |  | Polarization of Station 1 (p) and Station 2 (q) antenna (1 = vertical, 2 = horizontal)  |
| `press`               | scalar double | hPa |  | Surface pressure (default 1013.25 hPa) |
| `wat_vap_dens`               | scalar double | g/m^3 |  | Surface water-vapour density (default 8g/m^3) |
| `Rm`               | scalar double | mm/h |  | Peak rain rate at the centre of the rain cell (> 0.4 mm/hr) |
| `T`               | scalar double | deg C |  | Surface temperature (default 15 deg C) |
| `M1`, `M2`, `M3`               | scalar int | -|  | Number of integration points in the numerical integration over rho, phi, and z, respectively |


## Outputs ##

| Variable   | Type   | Units | Description |
|------------|--------|-------|-------------|
| `Lpq`    | double | dB    | Transmission loss due to hydrometeor scatterer |


## Software Versions
The code was tested and runs on:
* MATLAB version 2022a 


## References

* [Recommendation ITU-R P.452-18](https://www.itu.int/rec/R-REC-P.452/en)
* [Fascicle ITU-R 3M/FAS/12](https://www.itu.int/oth/R0A04000092/en)
<!-- * [ITU-R SG 3 Software, Data, and Validation Web Page](https://www.itu.int/en/ITU-R/study-groups/rsg3/Pages/iono-tropo-spheric.aspx) -->
