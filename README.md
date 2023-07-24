# MATLAB/Octave Implementation of Recommendation ITU-R P.452 (Hydrometeor Scatterer)


This code repository contains a MATLAB/Octave software implementation of draft revised Recommendation [ITU-R P.452-17](https://www.itu.int/md/R19-SG03-C-0127/en) with a prediction procedure for the evaluation of transmission loss due to hydrometeor scatterer for frequencies above about 0.1 GHz.  

This is a development code. It contains two versions:
* `tl_p452_hydrometeor_vec.m` which is based on the equations in a vector form (from draft revised Recommendation [ITU-R P.452-17](https://www.itu.int/md/R19-SG03-C-0127/en)  and the corresponding [draft Fascicle](https://www.itu.int/md/R19-WP3M-230522-TD-0166/en)) and 
* `tl_p452_hydrometeor.m` which is based on the equations in a scalar form from draft revised Recommendation [ITU-R P.452-17](https://www.itu.int/md/R19-SG03-C-0127/en)

Note that several scalar equations in draft revised Recommendation [ITU-R P.452-17](https://www.itu.int/md/R19-SG03-C-0127/en) have typos or are incorrect. Their direct implementation results in either incorrect or non-physical values (such as complex valued transmission loss in dB). A non-exhaustive list of issues and  related corrections is given below (they are already implemented in `tl_p452_hydrometeor.m`):

- $D_{A2}$  in (97), $d_{B2}$ in (101) have a wrong sign in front of the factor $\cos(\alpha_2-\varphi)$
- $\varphi_s$ in (115) has a wrong sign in the argument of $\cos^{-1}(\cdot)$
- $\alpha_{vi}$ in (124) is incorrect and has several missing factors. Using a corrected equation based on equations (85) and (86) from the [draft Fascicle](https://www.itu.int/md/R19-WP3M-230522-TD-0166/en) 
- The scalar product in equation (126) has a multiplicative factor $d$ missing. The same scalar product in (127) is incorrect and not used. Instead the corrected equation (126) is implemented
- The scalar product in (129) is incorrect and not used. Instead (128) is implemented 

Comparison between the vector form in `tl_p452_hydrometeor_vec.m` and the corrected scalar form in  `tl_p452_hydrometeor.m` is performed in `test_scalar_vs_vector_eqs.m` and those two versions result in numerically identical values.

A version of the draft revised Recommendation ITU-R P.452-17 with the above corrections introduced in track changes (together with other improvements and comments) is available in [this document](https://github.com/eeveetza/p452-hydrometeor/matlab/C0127-RevIS.docx).



## Function Call

~~~
Lpq = tl_p452_hydrometeor(d, f, lat1, lon1, lat2, lon2, h1loc, h2loc, alpha1loc, alpha2loc, eps1loc, eps2loc, G1, G2, BW1, BW2, p, q, press, wat_vap_dens, Rm, T, M1, M2, M3)
~~~


## Required input arguments of function `tl_p452_hydrometeor`

| Variable          | Type   | Units | Limits       | Description  |
|-------------------|--------|-------|--------------|--------------|
| ``               | scalar double |    | |   TBD |




 
## Outputs ##

| Variable   | Type   | Units | Description |
|------------|--------|-------|-------------|
| `Lpq`    | double | dB    | Transmission loss due to hydrometeor scatterer |


## Software Versions
The code was tested and runs on:
* MATLAB version 2022a 


## References

* [Recommendation ITU-R P.452](https://www.itu.int/rec/R-REC-P.452/en)
* [Draft revised Recommendation ITU-R P.452-17](https://www.itu.int/md/R19-SG03-C-0127/en) 
* [Draft Fascicle on Hydrometeor Scatterer](https://www.itu.int/md/R19-WP3M-230522-TD-0166/en)
* [ITU-R SG 3 Software, Data, and Validation Web Page](https://www.itu.int/en/ITU-R/study-groups/rsg3/Pages/iono-tropo-spheric.aspx)
