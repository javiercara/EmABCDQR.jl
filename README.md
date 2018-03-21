# emABCDQR

`EmABCDQR.jl` is a Julia package to estimate the state space model

x_{t+1} = Ax_{t} + Bu_{t} + w_{t}

y_{t} = Cx_{t} + Du_{t}  + v_{t}

where w_{t} -> N(0,Q), v_{t} -> N(0,R)

using the Expectation-Maximization algorithm.

## Installation

To install the package, from within Julia do

~~~
julia> Pkg.clone("git@github.com:javiercara/EmABCDQR.jl.git")
~~~

## Author

* **Javier Cara**, ETSI Industriales, Universidad Politecnica de Madrid (Spain)
