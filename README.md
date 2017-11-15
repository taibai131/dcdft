# DCDFT
## Date-compensated Discrete Fourier Transfrom

The unevenness in the spacing between consecutive dates and the small size of the series makes use of the techniques normally used for uniformly spaced time series (e.g. annual crops, rainfall, etc.) very hazardous. Therefore a date-compensated technique must be adopted.

The technique called DateCompensated Discrete Fourier Transform (DCDFT), corresponds to a curve-fitting approach using a sinusoid-plus-constant model, and is summarized below. For each trial frequency w, one coefficient of spectral correlation S is obtained by the following formulae:

![equation](http://www.sciweavers.org/download/Tex2Img_1510748170.jpg)