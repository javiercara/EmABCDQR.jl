# emABCDQR

`emABCDQR.jl` is a Julia package to estimate the state space model

$$x_{t+1} = Ax_t + Bu_t + w_t, \quad w_t \rightarrow N(0,Q)$$
$$y_t = Cx_t + Du_t + v_t, \quad v_t \rightarrow N(0,R)$$

using the Expectation-Maximization algorithm.

## Installation

To install the package, from within Julia do

~~~
julia> Pkg.add("git@github.com:javiercara/emABCDQR.jl.git")
~~~

## Author

* **Javier Cara**, ETSI Industriales, Universidad Politecnica de Madrid (Spain)

## Note

* Under development. To be used in a reasearch paper.