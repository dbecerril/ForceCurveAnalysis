# ForceCurveAnalysis
Set of functions written in Julia to obtain  Young Modulus from Force curves obtained with an AFM. 
Filters and smoothing functions are also included. A typical approach curve measured with the AFM micrsocope is of the form
![ac_1_](https://github.com/dbecerril/ForceCurveAnalysis/assets/22774966/1071db9d-8c45-4782-9d23-7ec1f1e43430)

By obtaining an approach curve at each pixel and analyzing each approach curve it is possible to obtain the elasticity at each pixel such as 
![testOut2](https://github.com/dbecerril/ForceCurveAnalysis/assets/22774966/1d024d3a-c94c-4e03-a993-62c9d7d6ef04)

Where we can cleary differantate the hard substrate (glass) in blue and the softer nucleus in red.

# Acknowledgment
D. Becerril acknowledges financial support from the Mexico City Ministry of Education, Science, Technology and Innovation (2023).
