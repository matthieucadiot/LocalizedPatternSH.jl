using MAT, JLD2, IntervalLinearAlgebra, FileIO, SparseArrays
using RadiiPolynomial, IntervalArithmetic, LinearAlgebra ; Interval(1)*Interval(1)

include("C:/Users/Matthieu/Desktop/code/julia/SwiftH/trace.jl")
include("C:/Users/Matthieu/Desktop/code/julia/SwiftH/isspd_interval.jl")
#
Base.abs2(x::Interval) = (y=abs(x);return y*y)
# # #
 N = 130 ;  d=interval(76.0) ; N0 = 90; 
 precision_param = 80 ;        setprecision(precision_param)
 μ = interval(0.28) ;  v1 = interval(-1.6) ; v2 = interval(1.0) 
μbig = interval(big(0.28)) ; v1big = interval(big(-1.6)) ; v2big =interval(big(1.0)) 
dbig = interval(big(76.0)) 

 fourier0 = CosFourier(N0, π/d)⊗CosFourier(N0, π/d)
 fourier = CosFourier(N, π/d)⊗CosFourier(N, π/d)
#
# W = load("Wo_028_70_130.jld2","W")
 
# prod_130=load("prod_130.jld2","prod_130") ; D_inv = load("D_inv.jld2","D_inv")

#  W = interval.(big.(coefficients(W)))
#  S = trace(N,80); C=S';
# U = (C.*D_inv)*solve(Matrix(prod_130),vec(S*W))
# display("U created")

# U = Sequence(fourier,vec(U))
# U2 = U*U
# jldsave("Wo_028_76_130_projected.jld2" ; U)
# jldsave("Wo_028_76_130_square.jld2" ; U2)

#
U = load("Wo_028_76_130_projected.jld2","U")
U2 = load("Wo_028_76_130_squared.jld2","U2")

# We build V0 in big float precision

V0_full = 2*v1big*U + 3*v2big*U2

# Conversion into usual Float64 precision
U = Interval.(Float64.(inf.(U),RoundDown),Float64.(sup.(U),RoundUp) )
V0_full = Interval.(Float64.(inf.(V0_full),RoundDown),Float64.(sup.(V0_full),RoundUp) )
U2 = Interval.(Float64.(inf.(U2),RoundDown),Float64.(sup.(U2),RoundUp) )
V0 = project(V0_full,CosFourier(2*N0, π/d)⊗CosFourier(2*N0, π/d))

display("computation of the square is done")

fourier0 = CosFourier(N0, π/d)⊗CosFourier(N0, π/d)
fourier = CosFourier(N, π/d)⊗CosFourier(N, π/d)


# Construction of the operator L and its inverse. Acutally, we only need to build their diagonal
∂1 = project(Derivative((2,0)), fourier0, fourier0,Interval{Float64})
∂2 = project(Derivative((0,2)), fourier0, fourier0,Interval{Float64})
Δ = copy(∂1) ;  radd!(Δ,∂2) ; 

∂1 = Nothing ; ∂2 = Nothing ; 

L = LinearOperator(fourier0,fourier0,Diagonal((diag(coefficients(Δ+I))).^2)) ;  Δ = Nothing
L = convert(Vector{Interval{Float64}},diag(coefficients(μ*I + L)))
Li = ones((N0+1)^2)./L

# Construction of the conversion operators from cosine to exponential and vice versa. In particular D_1 has terms (√α_n)_n on the diagonal and D2 has terms  (1/√α_n)_n on the diagonal. These operators allow to compute the norms defined on the spaces ℓᵖD₂.
D1 = convert(Vector{Interval{Float64}},exp2cos(N0))
D2 = ones((N0+1)^2)./D1

display("L, Li created")

DG = project(Multiplication(V0),fourier0, fourier0,Interval{Float64})
display("DG created")

# Construction of A
A = interval.(inv(I + mid.(DG).*mid.(Li)'))
display("A created")

 
# Computation of the norm of A
Aconj = LinearOperator(fourier0,fourier0,coefficients(A)')
norm_A = sqrt( opnorm(LinearOperator(coefficients(D1.*Aconj*A.*D2')),2) )
display(norm_A)

quant_fin = [0 norm_A]


############################### Z1 BOUND ##################################################################
# We present the computations of the different composents of Z1 

################### First part of Z1 ||I_d - A(I_d + DG(U_0)L^(-1))||_2 

lN = 1/(μ + (1-(π/d*(N0+1))^2)^2)
M1 = project(Multiplication(V0*V0),fourier0, fourier0, Interval{Float64})- DG*DG
n2 = opnorm(LinearOperator(coefficients(lN^2*D1.*A*M1*Aconj.*D2' )),2)
n1 = opnorm(LinearOperator(coefficients((I - D1.*A*(I +DG.*Li').*D2'))),2)

display("Ai created")

n_A = sqrt( n1^2 + n2 + lN^2*norm(V0,1)^2 + opnorm(LinearOperator(coefficients((D1.*Li).*M1.*(Li.*D2)')),2) )
display(n_A) ; M1 = Nothing; 

################## Second component of Z1 : computation of the Z1k's ###################

a = sqrt((-1+sqrt(1+μ)))/2
abig = sqrt((-1+sqrt(1+μbig)))/2
b = sqrt(interval(2))*a - im*sqrt(μ)/(2*sqrt(interval(2))*a)

############### Constants for μ = 0.28 ##########
#construction of the constant C0 (hat) described in Section 3.5.1
C0 = test_bessel_inequality(a,b,μ)

qq = [quant_fin C0]

fourierE = CosFourier(4*N0, π/d)⊗CosFourier(4*N0, π/d)

E1 = Sequence(fourierE ,interval.(big.(zeros((4*N0+1)^2))))
E2 = Sequence(fourierE ,interval.(big.(zeros((4*N0+1)^2))))
E12 = Sequence(fourierE ,interval.(big.(zeros((4*N0+1)^2))))
m1 =4;m2=4;m3=4;

for n1 = 1:m1*N0
    for n2 = 1:m1*N0
        if (n1 <= m3*N0)&&(n2 <= m3*N0)
            E12[(n1,n2)] = real((-1)^n1*( (1-exp(-4*abig*dbig))/(2*abig-im*n1*π/dbig) ))*real((-1)^n2*( (1-exp(-4*abig*dbig))/(2*abig-im*n2*π/dbig) ))  
        end 
    end
    E1[(n1,0)] = real((-1)^n1*( (1-exp(-4*abig*dbig))/(2*abig-im*n1*π/dbig) ))*(1-exp(-4*abig*dbig))
    
    if (n1 <= m3*N0)
    E12[(n1,0)] = real((-1)^n1*( (1-exp(-4*abig*dbig))/(2*abig-im*n1*π/dbig) ))*(1-exp(-4*abig*dbig))/(2*abig) 
    end
end

for n2 = 1:m1*N0
    E2[(0,n2)] = real((-1)^n2*( (1-exp(-4*abig*dbig))/(2*abig-im*n2*π/dbig) ))*(1-exp(-4*abig*dbig))
    
    if (n2 <= m3*N0)
    E12[(0,n2)] = real((-1)^n2*( (1-exp(-4*abig*dbig))/(2*abig-im*n2*π/dbig) ))*(1-exp(-4*abig*dbig))/(2*abig) 
    end
end

E1[(0,0)] = 1/(2*abig)*(1-exp(-4*abig*dbig)); E1 = E1/(2*dbig) ;
E2[(0,0)] = 1/(2*abig)*(1-exp(-4*abig*dbig)); E2 = E2/(2*dbig) ; 
E12[(0,0)] = 1/(4*abig^2)*(1-exp(-4*abig*dbig)) ; E12 = E12/(4*dbig^2)

E1 = Interval.(Float64.(inf.(E1),RoundDown),Float64.(sup.(E1),RoundUp) )
E2 = Interval.(Float64.(inf.(E2),RoundDown),Float64.(sup.(E2),RoundUp) )
E12 = Interval.(Float64.(inf.(E12),RoundDown),Float64.(sup.(E12),RoundUp) )
# ########################### Zu1 bound #####################################################

D12 = convert(Vector{Interval{Float64}},exp2cos(2*N0))
V = project((E1+E2)*V0,CosFourier(2*N0, π/d)⊗CosFourier(2*N0, π/d))
nV = coefficients(D12.*V0)'*coefficients(D12.*V) ; 

Zu1 = C0^2/a^2*4*d^2*nV
display(sqrt(Zu1))

qq = [qq Zu1]

# ########################### Zu2 bound #####################################################

b0, b1, b2 = computation_b(a,d)
b_normal = interval(sup(norm_A^2*4*d^2))
V = project((exp(-2*a*d)*b0*E1 + b1*E12 + exp(-2*a*d)*b2*E2)*(b_normal*V0),CosFourier(2*N0, π/d)⊗CosFourier(2*N0, π/d))
nV2 = maximum([0 coefficients(D12.*V0)'*coefficients(D12.*V)]) ; D12 = Nothing
 
Zu2 = Zu1+nV2*4*d^2/b_normal
display(sqrt(Zu2))

qq = [qq Zu2]
################################ Last component of Z1 : ||A(π_NV0)||_2

n_inf = norm(V0_full -V0,1)*norm_A ; V0_full = Nothing
display(n_inf)

################################### Computation of Z1

Z1 = n_A + 2*norm_A*sqrt(Zu1 + Zu2) + n_inf
display("Z1 has a value"); display(Z1)

qq = [qq n_A]
qq = [qq n_inf]
qq = [qq Z1]
################################### Z2 BOUND #####################################################

V = 2*(2*v1)*(6*v2*U) + 36*v2^2*U2
V[(0,0)] += 4*v1^2
GV02 = project(Multiplication(V),fourier0, fourier0, Interval{Float64}) 
GV02 = D1.*Aconj*GV02*A.*D2'
norm_VA = sqrt( opnorm(LinearOperator(coefficients(GV02)),2) ) ; GV02 = Nothing
r0 = interval(3e-4)

# β = interval(2.61)
# norm_l = sqrt( π*β/(β-1) + π^2*β/(4*sqrt(μ)) - π*β*atan(sqrt(μ)/(β-1))/(2sqrt(μ)) - π*β/(2*(μ*β^2/(β-1)^2 + β^2 )*(β-1) ) )

κ = sqrt( (2*sqrt(μ) + (1+ μ)*(2*π - 2*atan(sqrt(μ))))/(8*(1+μ)*μ^(3/2)) + 2*π^2/d*(interval(3)^(3/4)/(μ^(7/4)) + interval(3)/(μ^(5/2))) )

Z2 = [3*abs(v2)*norm_A*κ^2 ; κ*maximum([2*abs(v1) sqrt(norm_VA^2+norm(6*v2*U+2*v1,1)^2)])]

display("Z2 has a value"); display(Z2)
qq = [qq Z2[1] Z2[2]]

###################################### Y0 BOUND #######################################################

V = v1*U2 + v2*U2*U;

∂1 = project(Derivative((2,0)), fourier, fourier,Interval{Float64})
∂2 = project(Derivative((0,2)), fourier, fourier,Interval{Float64})
Δ = copy(∂1) ;  radd!(Δ,∂2) ; 

∂1 = Nothing ; ∂2 = Nothing

L = LinearOperator(fourier,fourier,Diagonal((diag(coefficients(Δ+I))).^2)) ;  Δ = Nothing
L = convert(Vector{Interval{Float64}},diag(coefficients(μ*I + L)))

 Y0 = 2*d*sup(sqrt( norm((A*(project(L.*U+V,fourier0))),2)^2 +  norm(L.*(U-project(U,fourier0)) + V - project(V,fourier0),2)^2 ) )
 display("Y0 has a value");display(Y0)

 qq = [qq Y0]


Z2test = sup(Z2[1])*1e-4 + sup(Z2[2])
 rmin=(1-sup(Z1) - sqrt((1-sup(Z1))^2-2*sup(Y0)*Z2test))/Z2test
 r0 = interval(rmin)
 Z2 = Z2[1]*r0 + Z2[2]

####################### Validation of the proof  #############################

if Z1 + Z2*r0 < 1
  if 1/2*Z2*r0^2 - (1-Z1) + Y0 < 0
    display("The proof was successful for r0 = ")
    display(sup(r0))  
  else
    display("failure: discriminant is negative")
  end
else
    display("failure: linear term is positive")
end



qq = [qq r0]

final_q_octo = qq
jldsave("final_q_octo.jld2" ; final_q_octo)






















# # b1 = 4*C0/a^2; b2 = 4*C1/a; b3 = 4*C2/(a*sqrt(a)) ; b4 = C0/a
# # E1 = project(E1*U0,CosFourier(3*N0, π/d)⊗CosFourier(3*N0, π/d))
# # E12 = project(E1,CosFourier(2*N0, π/d)⊗CosFourier(2*N0, π/d))
# # E1 = project(E1*U0,CosFourier(2*N0, π/d)⊗CosFourier(2*N0, π/d))

# # E1 = 2*v1*E12 + 3*v2*E1

# # E2 = project(E2*U0,CosFourier(3*N0, π/d)⊗CosFourier(3*N0, π/d))
# # E22 = project(E2,CosFourier(2*N0, π/d)⊗CosFourier(2*N0, π/d))
# # E2 = project(E2*U0,CosFourier(2*N0, π/d)⊗CosFourier(2*N0, π/d))

# # E2 = 2*v1*E22 + 3*v2*E2


# # D1 = convert(Vector{Interval{Float64}},diag(exp2cos(2*N0)))
# # D2 = ones((2*N0+1)^2)./D1
# # Cv0 = sqrt(coefficients((D1.*V0))'*coefficients(D1.*E1))+ sqrt(coefficients((D1.*V0))'*coefficients(D1.*E2))
# # D1 = convert(Vector{Interval{Float64}},diag(exp2cos(N0)))
# # D2 = ones((N0+1)^2)./D1


# # function g(x)
# #   return cosint(2*x)*sin(2*x) + (π/2 - sinint(2*x))*cos(2*x)
# # end

# # C_2_v0 = n_exp*sqrt(4*g(2*mid(a))*(g(4*mid(a))+π/2) + 8*g(mid(a))^2)
# # Cv0 = sqrt(Cv0^2 + C_2_v0^2)

# # c1 = π*b1^2/(4*d^2); c2 = d^2/(4*π^3)*(b2 + b3)^2
# # N_sum = floor( mid.((c2/c1)^(1/4)) )+1

# # S= interval(0.);
# # for n1 = -N_sum:N_sum
# #     for n2 = -N_sum:N_sum
# #         if n1^2+n2^2 <= N_sum^2
# #             S = S+1;
# #         end
# #     end
# # end

# # n_total = sqrt(b1^2/(4*d^2)*S + c2/interval(N_sum)^2 +  b4^2)*Cv0

# # display(n_total)














# #   W1 =mid.(U0)
# # # #
# # # # # # #
# #   a = sqrt((-1+sqrt(1+1/μ))/2)/sqrt(2)
# #  b = sqrt(2)*a - im/(2*sqrt(2)*a*μ)

# # # ############## μ = 0.32 ###############
# # # C0 =interval(0.703)
# # # C1 = interval(0.37)
# # # C2 = interval(0.375)

# # # ############## μ = 0.30 ###############
# # # # C0 =interval(0.68)
# # # # C1 = interval(0.335)
# # # # C2 = interval(0.38)

# # # # # ############# μ = 0.28 ################
# # C0 = interval(0.656)
# # C1 = interval(0.315)
# # C2 = interval(0.37)
# # # # #
# #    V = 2*v1*U0 + 3*v2*U0*U0;
# #    V1 = 2*v1*mid.(W) + 3*v2*mid.(W)*mid.(W)

# # # # #


# # fourierv = CosFourier(10, π/d)⊗CosFourier(10, π/d)

# # V = project(W,fourierv)

# #     @time C = compute_boundary_total_new(V,10,d,a,(μ),(C0),(C1),(C2))

# #      C = compute_boundary_total_new(mid.(big.(V1)),2*N0,big(d),mid(big(a)),mid(μ),mid(C0),mid(C1),mid(C2))
# #      C = compute_boundary_total_new(V0,2*N0,d,a,μ,C0,C1,C2)
# #     display((C))
# # # # # #
# # ∂1 = project(Derivative((2,0)), fourier0, fourier0,Interval{Float64})
# # ∂2 = project(Derivative((0,2)), fourier0, fourier0,Interval{Float64})
# # Δ = copy(∂1) ;  radd!(Δ,∂2) ; 

# # L = LinearOperator(fourier0,fourier0,Diagonal((diag(coefficients(Δ+I))).^2))
# # L = I + μ*L
# # Li = LinearOperator(fourier0,fourier0,Diagonal(ones((N0+1)^2)./diag(coefficients(L))))
# # # # #

# # G1 = project(Multiplication(U0),fourier0, fourier0,Interval{Float64})
# # G2 = project(Multiplication(U0*U0),fourier0, fourier0,Interval{Float64})
# # display("G1 G2 created")
# # # #
# #    D = LinearOperator(fourier0,fourier0,exp2cos(N0))
# #    Di = LinearOperator(fourier0,fourier0,Diagonal(ones((N0+1)^2)./diag(exp2cos(N0))))

# # # # # # # # #
# # # # # # # # #
# # # # # # # # #
# # # # # # # # # #
# # # # # # # # # # # # # ################ Z BOUND ######################################################
# # # # # # # # # #
# # # # #

# #    A = interval.(inv(I + (2*mid(v1)*mid.(G1) + 3*mid(v2)*mid.(G2))*mid.(Li))) ;
# # # #   display("A created")
# #  #   norm_A = maximum([1 opnorm(LinearOperator(coefficients(D*A*Di)),2)])
# #     B = D*A*Di; 
# #     B = coefficients(B)'*coefficients(B);
# #     eps_norm = 0.1
# #     c = interval.(opnorm(Hermitian(mid.(B)),2) + eps_norm)
    
# #     if isspd(Hermitian(c*I-B))
# #       norm_A = sqrt(c)
# #     else
# #       norm_A = false
# #     end 

# #      display(norm_A)
# # # # #   GC.gc()
# # # # # # # #
 
# #    dd = ones((N0+1)^2)./diag(coefficients(L));
# #    B = LinearOperator(fourier0,fourier0,Diagonal(dd))
# #  Df = I + (2*v1*G1+3*v2*G2)*B
# #  display("Df created")
# # # #   B= Nothing;  GC.gc();
# # # # # # #
# # # #   # 


# #  norm_I_A = sup( opnorm(D*(I  - A*Df)*Di ,2)) +  sup( norm(2*v1*U0+3*v2*U0*U0,1)*(1 + norm_A)/(1+μ*(1-((interval(pi)/d*interval(N+1)))^2)^2) )
# # display("value norm I A"); display(norm_I_A)
# # # #
# # Z = norm_A/(1-norm_I_A-norm_A*C)
# # display("value Z");display(Z)

# # # # # # # #
# # # # # # # # ################ Y BOUND ######################################################
# # # # # #
# # V = v1*W*W+v2*W*W*W;
# # ∂1 = project(Derivative((2,0)), fourier, fourier,Interval{Float64})
# # ∂2 = project(Derivative((0,2)), fourier, fourier,Interval{Float64})
# #   Δ = radd!(radd!(copy(∂1),I),∂2) ;

# #   d=76
# # V1 = 2*v1*W + 3*v2*W*W ; Q1 = mid.(project(Multiplication(V1),fourier2, fourier2,Interval{Float64}))
# # V2 = 2*v1*U2 + 3*v2*U2*U2 ; Q2 = mid.(project(Multiplication(V2),fourier2, fourier2,Interval{Float64}))
# # A1 = LinearOperator(fourier0,fourier0,convert(Matrix{Float64},Diagonal(ones((N0+1)^2))))
# # A1[(0:N2,0:N2),(0:N2,0:N2)] = coefficients(A)
# # C = norm_A*norm(V1-V2,1)/(1+μ*(1-((interval(pi)/d*interval(N2+1)))^2)^2) + opnorm(LinearOperator(coefficients(A*(Q1-Q2))),2)

# # L1 = ones((N+1)^2) + μ*diag(coefficients(Δ)).^2
# # Y = d*sup(norm((A*(L*U2+project(V,fourier2))),2) +  norm(L1.*W-L*U2 + V - project(V,fourier2),2) )
# # display("value Y");display(Y)
# # # # # # #
# # # # # # # # ################ Z2 BOUND ######################################################
# #  C0 =2.8
# # Z2 = 1/(1-norm_I_A-C)*sup(Z*C0*(abs(v1) +3*v2*norm(W,1) ))/2
# # # # # # #
# # # # # # # # # ############### Verification ###################################################
# # # # # # # # # r_star = 1.0;
# # # # # # # # # interval_proof = interval_of_existence(Y,Interval(0),Z2,r_star,C²Condition())
# # # # # # # # # display(interval_proof)
# # # # # # # # #
# # if 1- norm_I_A>0
# #   if 1-4*Y*Z2 > 0
# #     rmin=(1 - sqrt(1-4*Y*Z2))/(2*Z2);
# #     rmax=(1 + sqrt(1-4*Y*Z2))/(2*Z2);
# #     if rmin<rmax
# #       display(rmin)
# #       display("success")
# #     else
# #       display("rmin>=rmax")
# #     end
# #   else
# #     display("failure: discriminant is negative")
# #   end
# # else
# #     display("failure: linear term is positive")
# # end
