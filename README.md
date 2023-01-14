# How to make an affine coordinate shallow water model
mkdir build_affine
cd build_affine
cmake .. -DEquationOfState=shallow_water -DUseAffine=ON

# cartesian shallow water mode
mkdir build_cart
cd build_cart
cmake .. -DEquationOfState=shallow_water
