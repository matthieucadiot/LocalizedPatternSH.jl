# LocalizedPatternSH
Julia code for computer-assisted proofs of localized patterns in the Swift-Hohenberg PDE


# Computer-assisted proofs of solitons, eigencouples and orbital stability in the Kawahara equation.



Table of contents:


* [Introduction](#introduction)
* [The Kawahara equation](#the-kawahara-equation)
   * [Proof of solitons](#Proof-of-solitons)
   * [Proof of the first 3 eigencouples](#Proof-of-the-first-3-eigencouples)
   * [Proof of orbital stability](#Proof-of-orbital-stability)
* [Utilisation and References](#utilisation-and-references)
* [License and Citation](#license-and-citation)
* [Contact](#contact)



# Introduction

This Julia code is a complement to the article 

#### [[1]](https://arxiv.org/abs/2302.12877) : "Rigorous computation of solutions of semi-linear PDEs on unbounded domains via spectral methods", M. Cadiot, J-P. Lessard and J-C. Nave, [ArXiv Link](https://arxiv.org/abs/2302.12877)

as it provides the necessary rigorous computations that are needed in Section 6. The rigorous computations are performed using the package [IntervalArithmetic](https://github.com/JuliaIntervals/IntervalArithmetic.jl). The mathematical objects (spaces, sequences, operators,...) are built using the package [RadiiPolynomial](https://github.com/OlivierHnt/RadiiPolynomial.jl). 


# The Kawahara equation

The Kawahara equation
$$\lambda_2u'''' + \lambda_1u'' + u + \lambda_3u^2 = 0$$
is known to have solutions on $\mathbb{R}$ that decay to zero at infinity. These solutions are called solitary waves or soliton (see [[1]](https://arxiv.org/abs/2302.12877) for a complete description).

## Proof of solitons

The present code provides the rigorous numerics for the proof of solitons of the Kawahara equation using the analysis of [[1]](https://arxiv.org/abs/2302.12877) (specifically the Section 6). The user can choose, line 355 of the code, the values for N, d, T and c, that are described in [[1]](https://arxiv.org/abs/2302.12877). In particular, T and c need to be chosen such that
 - $0 \leq T < 0.397$ 
 - $c < 1- \frac{a(T)^2}{4b(T)}$    
where
- $a(T) = \frac{1-3T}{6}$
- $b(T) = \frac{19 - 30T - 45T^2}{360}$.   

The code will compute rigorously the needed bounds of Section 6 of [[1]](https://arxiv.org/abs/2302.12877) and validate of not the computer-assisted proof. If the computer-assisted proof succeeds, the radius for the smallest and biggest ball of contraction is displayed.

## Proof of the first 3 eigencouples

If the proof of the soliton is achieved, the code will then compute approximations for the first 3 eigencouples of the linearization around the proved soliton. Then, the needed bounds for the proof of eigencouples are computed, following the analysis of Section 6 of [[1]](https://arxiv.org/abs/2302.12877). In particular, the proof will be valid if T was chosen such that
 $\frac{1}{3} < T < 0.397$. 
  
   

 If the computer-assisted proof succeeds, the radius for the smallest and biggest ball of contraction is displayed for each eigencouple.
 
 
 ## Proof of orbital stability

If the proof of the first 3 eigencouples is achieved, the code will then try to prove Theorem 6.15 in [[1]](https://arxiv.org/abs/2302.12877). In particular, we want to prove that the 3 eigenvalues obtained beforehead, are actually the first 3 ones. The algorithm is explained in the proof of Theorem 6.15 and uses Proposition 6.14. 

 If the computer-assisted proof of Theorem 6.15 succeeds, the value for  <img src="https://latex.codecogs.com/gif.latex?\tau" /> (Proposition 6.14) is computed rigorously. In particular, we check that <img src="https://latex.codecogs.com/gif.latex?\tau" /> is striclty negative, and if that is the case, then we obtain that the proved soliton is orbitally stable.
 
 
 # Utilisation and References
 
 Go directly to line 355 of the code in order to enter the values for N, d, c and T. The values need to respect the requirements cited above. I would suggest to keep the by-default values of N and d. However, if the proof fails, you can try to increase the value of d and N. This can happen if you choose values of c and T that are close to singular values (bifurcation points) such as $c =1$ or $T = \frac{1}{3}$.
 
 The code is build using the following packages :
 - [RadiiPolynomial](https://github.com/OlivierHnt/RadiiPolynomial.jl) 
 - [IntervalArithmetic](https://github.com/JuliaIntervals/IntervalArithmetic.jl)
 - [LinearAlgebra](https://docs.julialang.org/en/v1/stdlib/LinearAlgebra/)
 - [FFTW](https://github.com/JuliaMath/FFTW.jl)
 - [PrettyTables](https://ronisbr.github.io/PrettyTables.jl/stable/).
 
 You will need to install these packages. After that, you should be able to run the code normally.
 
 # License and Citation
 
  This code is available as open source under the terms of the [MIT License](http://opensource.org/licenses/MIT).
  
If you wish to use this code in your publication, research, teaching, or other activities, please cite it using the following BibTeX template:

```
@software{ProofKawahara.jl,
  author = {Matthieu Cadiot},
  title  = {ProofKawahara.jl},
  url    = {https://github.com/matthieucadiot/ProofKawahara.jl},
  note = {\url{ https://github.com/matthieucadiot/ProofKawahara.jl},
  year   = {2023},
  doi = {10.5281/zenodo.7656856}
}
```
DOI : [10.5281/zenodo.7656856](https://zenodo.org/record/7656856#.Y_NqJ0NKhPZ) 


# Contact

You can contact me at :

matthieu.cadiot@mail.mcgill.ca

