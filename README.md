# Computer-assisted proofs of localized patterns and branches of periodic solutions in the Swift-Hohenberg PDE.



Table of contents:


* [Introduction](#introduction)
* [The Swift Hohenberg equation](#the-swift-hohenberg-equation)
   * [Proof of a branch of periodic solutions limiting a localized pattern](#proof-of-a-branch-of-periodic-solutions-limiting-a-localized-pattern)
* [Utilisation and References](#utilisation-and-references)
* [License and Citation](#license-and-citation)
* [Contact](#contact)



# Introduction

This Julia code is a complement to the article 

#### [[1]](https://arxiv.org/abs/2302.12877) : "Stationary non-radial localized patterns in the planar Swift-Hohenberg PDE: constructive proofs of existence", M. Cadiot, J-P. Lessard and J-C. Nave, [ArXiv Link](https://arxiv.org/abs/2302.12877)

as it provides the necessary rigorous computations that are needed along the paper. The rigorous computations are performed using the package [IntervalArithmetic](https://github.com/JuliaIntervals/IntervalArithmetic.jl). The mathematical objects (spaces, sequences, operators,...) are built using the package [RadiiPolynomial](https://github.com/OlivierHnt/RadiiPolynomial.jl). 


# The Swift-Hohenberg equation

The Swift-Hohenberg equation
$$(I_d+\Delta)^2u +  \mu u + \nu_1 u^2 + \nu_2 u^3 =0$$
is known to have localized solutions on $\mathbb{R}^2$ that vanish at infinity. These solutions are called localized patterns (see [[1]](https://arxiv.org/abs/2302.12877) for an introduction to the subject). 

## Proof of a branch of periodic solutions limiting a localized pattern

In Section 3.6 in [[1]](https://arxiv.org/abs/2302.12877), we prove that, under the condition (67) in Theorem 3.7, localized patterns can be proven to be the limit of a branch of (spatially) periodic solutions as the period tends to infinity. In particular this condition involves the explicit computation of multiple bounds, which is achieved in the present code. Moreover, we verify that condition (67) is verified for 3 different localized patterns and we obtain a constructive proof of existence of a branch of periodic solutions limiting the localized pattern. These 3 patterns possess different visual symmetries : a square symmetry (associated to the group D4), an hexagonal symmetry (D6) and an octogonal symmetry (D8). The symmetries are not proven but simply indicative. In particular, each code proof_Dj_pattern.jl allows to demonstrate the branch of periodic solutions and the localized pattern associated. 

We provide as well candidate solutions for the proofs, which are given in the files .jld2. These correspond to the sequence U0 in Section 3.1 representing the approximate solution u0. In particular, $U_0$ has already been projected in the set of sequences representing trace zero functions (see Section 3.1). Consequently, the Fourier series associated to $U_0$ represents a smooth function on $\mathbb{R}^2$ with compact support on a square.

Given these approximate solution, each code proof_Dj_pattern.jl provides the explicit computation of the bounds in Theorem 3.7 and provides a value for $r$ is the proof is successful. In particular, the theorem states that there exists a smooth curve 
$$\{\tilde{u}(q) : q \in [d,\infty]\} \subset C^\infty(\mathbb{R}^2)$$
such that $\tilde{u}(q)$ is a periodic solution to the Swift-Hohenberg PDE with period $2q$ in both variables.  In particular, $\tilde{u}(\infty)$ is a localized pattern on $\mathbb{R}^2.$ Finally, the value of $r$ provides a uniform control on the branch of periodic solutions, making the proof constructive.

 
 # Utilisation and References

 The code in proof_Dj_pattern.jl can serve to prove other patterns than the one provided as illustration. The interested user can go directly to the line 440 of the code an modify the values of the parameters accordingly. In particular, the projection in the set of functions with null trace is commented by default but can be used if needed. 
 
 The code is build using the following packages :
 - [RadiiPolynomial](https://github.com/OlivierHnt/RadiiPolynomial.jl) 
 - [IntervalArithmetic](https://github.com/JuliaIntervals/IntervalArithmetic.jl)
 - [LinearAlgebra](https://docs.julialang.org/en/v1/stdlib/LinearAlgebra/)
 - [IntervalLinearAlgebra](https://github.com/JuliaIntervals/IntervalLinearAlgebra.jl)
 - [JLD2](https://github.com/JuliaIO/JLD2.jl)
 
 
 # License and Citation
 
This code is available as open source under the terms of the [MIT License](http://opensource.org/licenses/MIT).
  
If you wish to use this code in your publication, research, teaching, or other activities, please cite it using the following BibTeX template:

```
@software{ProofKawahara.jl,
  author = {Matthieu Cadiot},
  title  = {LocalizedPatternSH.jl},
  url    = {https://github.com/matthieucadiot/LocalizedPatternSH.jl},
  note = {\url{ https://github.com/matthieucadiot/LocalizedPatternSH.jl},
  year   = {2024},
  doi = {10.5281/zenodo.7656856}
}
```
DOI : [10.5281/zenodo.7656856](https://doi.org/10.5281/zenodo.10823084) 


# Contact

You can contact me at :

matthieu.cadiot@mail.mcgill.ca

