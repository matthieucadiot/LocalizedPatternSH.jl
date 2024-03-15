using JLD2, IntervalLinearAlgebra, RadiiPolynomial, IntervalArithmetic, LinearAlgebra


function trace(N,p)
    setprecision(p)
    M = (N+1)^2;
    D10 = interval.(big.(zeros(N+1,M))); D12 = interval.(big.(zeros(N+1,M))); D20 = interval.(big.(zeros(N+1,M)));
    D22 = interval.(big.(zeros(N+1,M)))
    w = interval.(big.(0:N));
    Y = w.^2;
    X = interval.(big.((-1).^(0:N)));
    k=2;

    for n2 in 0:N
        D10[n2+1,(n2*(N+1)+1):((n2+1)*(N+1))] = 2*X;
        D10[n2+1,n2*(N+1)+1] = interval.(big.(1));
        D12[n2+1,(n2*(N+1)+1):((n2+1)*(N+1))] = 2*Y.*X;
    end

    for n1 in 0:N
        D20[n1+1,(n1+1):(N+1):M] = 2*X;
        D20[n1+1,n1+1] = interval.(big.(1));
        D22[n1+1,(n1+1):(N+1):M] = 2*Y.*X;
    end
    # we remove some lignes to have a full rank operator

    f = ([D10[(k+1):(N+1),:];
                D12[(k+1):(N+1),:] ;
                D20;
                D22])
    return f
end

function compute_boundary_d1(V,N,l,a)
S = interval(.0);

@inbounds for n2 = 1:N
    @inbounds for n1 = 1:N
        @inbounds for m1 = 1:N
               b = 4*a*(1/((pi/l*(n1-m1))^2+4*a^2) + 1/((pi/l*(n1+m1))^2+4*a^2)) ;
               u= V[(n1,n2)];

               v= V[(m1,n2)];
               S = S+ 2*(-1)^(n1-m1)*u*v*b;
       end

        b =  4*a/((pi/l*(n1))^2+4*a^2);
        u= V[(n1,n2)];
        v= V[(0,n2)];
        S = S+ 2*(-1)^(n1)*u*v*b;
    end


    @inbounds for m1 = 1:N
               b = 4*a/((pi/l*(m1))^2+4*a^2) ;
               u= V[(0,n2)];
               v= V[(m1,n2)];
               S = S+ 2*(-1)^(m1)*u*v*b;
       end

        b =  1/(2*a) ;
        u= V[(0,n2)];
        v= V[(0,n2)];
        S = S+ 2*u*v*b;
end


@inbounds for n1 = 1:N
    @inbounds for m1 = 1:N
               b = 4*a*(1/((pi/l*(n1-m1))^2+4*a^2) + 1/((pi/l*(n1+m1))^2+4*a^2)) ;
               u= V[(n1,0)];
               v= V[(m1,0)];
               S = S+ (-1)^(n1-m1)*u*v*b;
       end

        b =  4*a/((pi/l*(n1))^2+4*a^2);
        u= V[(n1,0)];
        v= V[(0,0)];
        S = S+ (-1)^(n1)*u*v*b;
end


@inbounds for m1 = 1:N
               b = 4*a/((pi/l*(m1))^2+4*a^2) ;
               u= V[(0,0)];
               v= V[(m1,0)];
               S = S+ (-1)^(m1)*u*v*b;
   end

        b =  1/(2*a) ;
        u= V[(0,0)];
        v= V[(0,0)];
        S = S+ u*v*b;

f = (sqrt(abs(S)));
end


function compute_boundary_d2(V,N,l,a)
S = interval(.0);

@inbounds for n2 = 1:N
   @inbounds for n1 = 1:N
       @inbounds for m1 = 1:N
               b = 4*a*(1/((pi/l*(n1-m1))^2+4*a^2) + 1/((pi/l*(n1+m1))^2+4*a^2)) ;
               u= V[(n2,n1)];

               v= V[(n2,m1)];
               S = S+ 2*(-1)^(n1-m1)*u*v*b;
       end

        b =  4*a/((pi/l*(n1))^2+4*a^2);
        u= V[(n2,n1)];
        v= V[(n2,0)];
        S = S+ 2*(-1)^(n1)*u*v*b;
    end


      @inbounds for m1 = 1:N
               b = 4*a/((pi/l*(m1))^2+4*a^2) ;
               u= V[(n2,0)];
               v= V[(n2,m1)];
               S = S+ 2*(-1)^(m1)*u*v*b;
       end

        b =  1/(2*a) ;
        u= V[(n2,0)];
        v= V[(n2,0)];
        S = S+ 2*u*v*b;
end


@inbounds for n1 = 1:N
    @inbounds for m1 = 1:N
               b = 4*a*(1/((pi/l*(n1-m1))^2+4*a^2) + 1/((pi/l*(n1+m1))^2+4*a^2)) ;
               u= V[(0,n1)];
               v= V[(0,m1)];
               S = S+ (-1)^(n1-m1)*u*v*b;
       end

        b =  4*a/((pi/l*(n1))^2+4*a^2);
        u= V[(0,n1)];
        v= V[(0,0)];
        S = S+ (-1)^(n1)*u*v*b;
end


@inbounds for m1 = 1:N
               b = 4*a/((pi/l*(m1))^2+4*a^2) ;
               u= V[(0,0)];
               v= V[(0,m1)];
               S = S+ (-1)^(m1)*u*v*b;
   end

        b =  1/(2*a) ;
        u= V[(0,0)];
        v= V[(0,0)];
        S = S+ u*v*b;

f = (sqrt(abs(S)));
end



function compute_boundary_total_new(V,N,d,a,μ,C0,C1,C2)


######## computation of the constants A_0(w)
A0 =   sqrt(2*d)*(compute_boundary_d1(V,N,d,a) + compute_boundary_d2(V,N,d,a))          #A_0(v_0)
display(A0)

####### computation of C1 and C3
C_0 = 4*C0/a^2*A0
C_1 = 4*C1/a*A0
C_2 = 4*C2*sqrt(interval(π))/(a*sqrt(a))*A0
C_3 = C0/(a)*A0

######## computation of C

b = 4*d^2/(16*pi^3)*(C_1+C_2)^2
a = pi/(4*d^2)*C_0^2
N = floor((b/a)^(1/4))


I1 = sqrt( a*N^2 + b/N^2)
I2 = C_3

display(C_0)
display(C_1)
display(C_2)
display(C_3)
display(I1)
C = sqrt(I1^2 + I2^2)

    return C
end


function symmetry_SH(N,a,b)

f= zeros((N+1)^2,(N+1)^2);
for n1=0:N
    for n2=0:N

        for k1=1:N
            for k2=1:N
                u1 = sinc((a*k1+b*k2-n1))*sinc((a*k2-b*k1-n2));
                u2 = sinc((-a*k1+b*k2-n1))*sinc((a*k2+b*k1-n2));
                u3 = sinc((a*k1-b*k2-n1))*sinc((-a*k2-b*k1-n2));
                u4 = sinc((-a*k1-b*k2-n1))*sinc((-a*k2+b*k1-n2));
                f[n2*(N+1)+n1+1,k2*(N+1)+k1+1] = u1 + u2 + u3+ u4;
            end


                u1 = sinc((a*k1-n1))*sinc((-b*k1-n2));
                u2 = sinc((-a*k1-n1))*sinc((b*k1-n2));
                f[n2*(N+1)+n1+1,k1+1] = u1 + u2;


        end
        for k2=1:N
                u1 = sinc((b*k2-n1))*sinc((a*k2-n2));
                u3 = sinc((-b*k2-n1))*sinc((-a*k2-n2));
                f[n2*(N+1)+n1+1,k2*(N+1)+1] = u1 + u3;
        end

        u1 = sinc((-n1))*sinc((-n2));
        f[n2*(N+1)+n1+1,1] = u1;

    end
end
f= Matrix(I, (N+1)^2, (N+1)^2)-f;

return f
end



function exp2cos(N)

    d = 2*(interval.(ones((N+1)^2)))

    d[1] = interval(1);
    for n2=1:N
        d[n2+1] = interval(sqrt(2));
    end

    for n1 = 1:N
        d[n1*(N+1)+1] = interval(sqrt(2));
    end

    return d
end


function f_0(b,x)
    
    return ( b^2*(x^2+1) - 2*b*x +2 )/b^3
end

function f_1(b,x)
    
    return ( b^2*(x^2+1) - 2*b*x +2 )/b^3
end

function f_2(b,x)
   
    return ( b^2*(x^2+1) - 2*b*x +2 )/b^3
end

function f_0_int(a,n,d)
    n1 = n[1]
    n2 = n[2]
    if n2==0 
        return (-1)^(n1)*(2*d+2/3*d^3)*real(f_0(2*a+im*n1*π/d,d)-exp(-4*a*d)*f_0(2*a+im*n1*π/d,-d))
    else    
        return (-1)^(n2)*(-1)^(n1)*real(f_0(2*a+im*n1*π/d,d)-exp(-4*a*d)f_0(2*a+im*n1*π/d,-d))*real(f_0(im*n2*π/d,d)-exp(-4*a*d)*f_0(im*n2*π/d,-d)) 
    end
end

function f_1_int(a,n,d)
    n1 = n[1]
    n2 = n[2]
    if n1==0 
        return (-1)^(n2)*(2*d+2/3*d^3)*real(f_1(2*a+im*n2*π/d,d)-exp(-4*a*d)*f_1(2*a+im*n2*π/d,-d))
    else    
        return (-1)^(n2)*(-1)^(n1)*real(f_1(im*n1*π/d,d)-exp(-4*a*d)f_1(im*n1*π/d,-d))*real(f_1(2*a+im*n2*π/d,d)-exp(-4*a*d)*f_1(2*a+im*n2*π/d,-d)) 
    end
end

function f_2_int(a,n,d)
    n1 = n[1]
    n2 = n[2]
    return (-1)^(n2)*(-1)^(n1)*real(f_2(2*a+im*n1*π/d,d)-exp(-4*a*d)f_2(2*a+im*n1*π/d,-d))*real(f_2(2*a+im*n2*π/d,d)-exp(-4*a*d)*f_2(2*a+im*n2*π/d,-d)) 
end


function f_exp1(a,n1,d)
    return  (-1)^(n1)*real(f_2(2*a+im*n1*π/d,d)-exp(-4*a*d)*f_2(2*a+im*n1*π/d,-d))
end


function f_exp2(a,n2,d)
    return  (-1)^(n2)*real(f_2(2*a+im*n2*π/d,d)-exp(-4*a*d)*f_2(2*a+im*n2*π/d,-d))
end


function computation_b(a,d)
    x = 2*a*d
    c1 = (x+ 1 + exp(-x))/(a^2) + exp(-x)*(4*d + exp(-x)/a) + (2*exp(1)+1)/(a*(1-exp(x/2)))*(2*d+ (1+exp(-x))/a ) + (2*exp(1) +1)^2/(a^2*(1-exp(-x/2))^2)
    b1 = 2/a*( (1+exp(-x))/a + 2*d + exp(-x)*(4*d + exp(-x)/a) + (2*exp(1)+1)/(a*(1-exp(-x/2))) )

    C1 = 4*c1 + b1

    c2 = 2*d + (1+exp(-x))/(2*a) + (2*d + (1+exp(-x))/(2*a))*1/(1-exp(-x)) + (4*exp(1) + 1 + exp(-x))/(2*a*(1-exp(-x/2))^2)
    b2 = (2*d + 1/a)*8/(a*(1-exp(-x)))

    C12 = 8*c2*(2*d + 1/(2*a)) + b2

    C2 = 2/a*( (1+exp(-x))/a  + 2*d + exp(-x)*(4*d + exp(-x)/a) + (2*exp(1) + exp(-x))/(a*(1-exp(-x/2))))

    return C1, C12, C2
end


function bessel_conj(r,θ,μ)

    z = abs(r*b)

    γ = interval(0.57721566490153286060,0.57721566490153286061) #construction of the Euler constant
    Hn = interval(0.) #sum of the 1/n
    f = interval(0.) #initialization of f
    N = 10

    for m=0:N
        
        c = z^(2*m)/(interval(4.0)^m*interval(factorial(m))^2)*( (Hn - γ -log(z/2))*sin(2*m*θ) - θ*cos(2*m*θ) )
        f = f+c 
        Hn = Hn + 1/(interval(m)+1)
    end
    
     z = (4 + sup(abs(log(z/2))) + abs(θ))/(3*4^N*interval(factorial(N+1)))
     f = 1/sqrt(μ)*(f+z)
    
    return f
    
end
    

function test_bessel_inequality(a,b,μ)

    # Definition of the constants introduced in Section 3.5.1 and Proposition 3.1
    θ = -atan(sqrt(μ)/(4*a^2))
    β = π/(2*sqrt(2)*a*θ^2)
    δ = 0.0001
    ϵ = minimum([1/abs(b) sqrt(μ)*δ/((exp(interval(1/4))-1)*abs(b)*(4+exp(-interval(1))-θ))])

    #initialization of C0, in particular, by construction this constant deals with the case [0,ϵ] in the  verification of C0
    C0 = exp(sqrt(2)*a*ϵ)*(δ - θ/sqrt(μ)) 

    # Choice of the length of the intervals (In) used in the computer assisted verification (as described in Section 3.5.1)
    di = 0.001
    
    for i=inf(ϵ):di:sup(β)
        In = interval(i,i+di)
        f = bessel_conj(In,θ,μ) #construction of the rigorous enclosure of f_0(In)
        C0 = maximum([C0 sup(f*exp(sqrt(interval(2))*a*In))]) # we choose C0 as the max over the previous C0 and the supremum of g(In)
    end

    return C0

end



function change_domain(U,N1,N2,d1,d2)

    W = Sequence(CosFourier(N2, π/d2)⊗CosFourier(N2, π/d2),vec(zeros((N2+1)^2)))
    if d2<d1
        for m1 ∈ (0:N2)
            for m2 ∈ (0:N2)
                for n1 ∈ (-N1:N1)
                    for n2 ∈ (-N1:N1)
                        W[(m1,m2)] = W[(m1,m2)] + U[(abs(n1),abs(n2))]*sinc(d2*n1/d1 - m1)*sinc(d2*n2/d1-m2)
                    end
                end
            end
        end
    else
        for m1 ∈ (0:N2)
            for m2 ∈ (0:N2)
                for n1 ∈ (-N1:N1)
                    for n2 ∈ (-N1:N1)
                        W[(m1,m2)] = W[(m1,m2)] + d1^2/(d2^2)*U[(abs(n1),abs(n2))]*sinc(d2*n1/d1 - m1)*sinc(d2*n2/d1-m2)
                    end
                end
            end
        end

    end    
    return W
    
end



function PlotCoeffs2D(U0,a,b,c,d)
    #U0 is a sequence in 2D
    #a,b,c,d are the endpoints of the interval [a,b] × [c,d]
    y1 = a:0.1:b
    y2 = c:0.1:d
    m=length(y1)
    U = complex(zeros(m,m))
    for b₁ = 1:m
        for b₂ = 1:m
            U[b₁,b₂] = U0(y1[b₁],y2[b₂])
        end
    end
    U = real.(U)
    mat"
    h = surf($y1,$y2,$U)
    set(h,'LineStyle','none')"
end















###################################################################################################################################
                     MAIN CODE FOR THE PROOF
###################################################################################################################################



# Set up of the parameters
 N0 = 130 ;  d=interval(76.0) ; N = 90; 
 precision_param = 80 ;        setprecision(precision_param)
 μ = interval(0.28) ;  v1 = interval(-1.6) ; v2 = interval(1.0) 
μbig = interval(big(0.28)) ; v1big = interval(big(-1.6)) ; v2big =interval(big(1.0))  # using big float for computations in multi-precision
dbig = interval(big(76.0)) 

# Set up of the spaces for the sequences of Fourier coefficients
 fourier0 = CosFourier(N, π/dbig)⊗CosFourier(N, π/dbig)
 fourier = CosFourier(N0, π/dbig)⊗CosFourier(N0, π/dbig)

#############################################################
# The code below can be uncommented if one wants to compute the projection of an approximate solution W in the
# set of functions with zero trace (cf. Section 3.1). 


#   W = interval.(big.(coefficients(W)))
#   S = trace(N0,80); C=S';    # trace operator and its adjoint for a precision 80digits in big float

#   ∂1 = project(Derivative((2,0)), fourier, fourier,Interval{BigFloat})
#   ∂2 = project(Derivative((0,2)), fourier, fourier,Interval{BigFloat})
#   Δ = copy(∂1) ;  radd!(Δ,∂2) ; 
  
#   ∂1 = Nothing ; ∂2 = Nothing
  
#   L = LinearOperator(fourier,fourier,Diagonal((diag(coefficients(Δ+I))).^2)) ;  Δ = Nothing
#   L = convert(Vector{Interval{Float64}},diag(coefficients(I + μbig*L)))
#   Li = (ones((N+1)^2)./L).^2

#   prod_130 = S*(Li.*C)

#  U = W - (Li.*C)*solve(Matrix(prod_130),vec(S*W))   # construction of the projection using the fomula of Section 3.1
#  display("U created")

#  prod_130_2 = Nothing ; D_inv_2 = Nothing; W = Nothing; S= Nothing; C = Nothing;
 
#  U = Sequence(fourier,vec(U))


# Instead, we upload an approximate solution U which we already projected in order to have a null trace. U2 is the square of U, which we also 
# already computed in multi-precision for convenience

 U = load("Wo_028_76_130_projected.jld2","U")
 U2 = load("Wo_028_76_130_square.jld2","U2")

# We build V0 in big float precision
V0_full = 2*v1big*U + 3*v2big*U2

# Conversion into usual Float64 precision
U = Interval.(Float64.(inf.(U),RoundDown),Float64.(sup.(U),RoundUp) )
V0_full = Interval.(Float64.(inf.(V0_full),RoundDown),Float64.(sup.(V0_full),RoundUp) )
U2 = Interval.(Float64.(inf.(U2),RoundDown),Float64.(sup.(U2),RoundUp) )
V0 = project(V0_full,CosFourier(2*N, π/d)⊗CosFourier(2*N, π/d))

# The set of sequences are re-defined in standard precision now
fourier0 = CosFourier(N, π/d)⊗CosFourier(N, π/d)
fourier = CosFourier(N0, π/d)⊗CosFourier(N0, π/d)


# Construction of the operator L and its inverse. Acutally, we only need to build their diagonal
∂1 = project(Derivative((2,0)), fourier0, fourier0,Interval{Float64})
∂2 = project(Derivative((0,2)), fourier0, fourier0,Interval{Float64})
Δ = copy(∂1) ;  radd!(Δ,∂2) ; 

∂1 = Nothing ; ∂2 = Nothing ; 

L = LinearOperator(fourier0,fourier0,Diagonal((diag(coefficients(Δ+I))).^2)) ;  Δ = Nothing
L = convert(Vector{Interval{Float64}},diag(coefficients(μ*I + L)))
Li = ones((N+1)^2)./L

# Construction of the conversion operators from cosine to exponential and vice versa. In particular D_1 has terms (√α_n)_n on the diagonal
# and D2 has terms  (1/√α_n)_n on the diagonal. These operators allow to compute the norms defined on the spaces ℓᵖD₂.
D1 = convert(Vector{Interval{Float64}},exp2cos(N))
D2 = ones((N+1)^2)./D1

# Frechet derivative of G
DG = project(Multiplication(V0),fourier0, fourier0,Interval{Float64})


# Construction of B^N, the approximate inverse of (the finite projection) DF(U0)
B = interval.(inv(I + mid.(DG).*mid.(Li)'))

# Computation of the norm of B
Bconj = LinearOperator(fourier0,fourier0,coefficients(B)') # construction of the adjoint
norm_B = maximum( [1 sqrt( opnorm(LinearOperator(coefficients(D1.*Bconj*B.*D2')),2) )] )  # cf. equation (25)


###############################  𝒵1 BOUND ##################################################################
# We present the computations of the different composents of  𝒵1 

################### The bound Z1 (Lemma 3.5)

lN = 1/(μ + (1-(π/d*(N+1))^2)^2)
M1 = project(Multiplication(V0*V0),fourier0, fourier0, Interval{Float64})- DG*DG
n2 = opnorm(LinearOperator(coefficients(lN^2*D1.*B*M1*Bconj.*D2' )),2)
n1 = opnorm(LinearOperator(coefficients((I - D1.*B*(I +DG.*Li').*D2'))),2)

Z1 = sqrt( n1^2 + n2 + lN^2*norm(V0,1)^2 + opnorm(LinearOperator(coefficients((D1.*Li).*M1.*(Li.*D2)')),2) )
display(n_A) ; M1 = Nothing; 

################## Computation of 𝒵u (Section 3.5) ###################

# First define the constants of Section 3.5
a = sqrt((-1+sqrt(1+μ)))/2
abig = sqrt((-1+sqrt(1+μbig)))/2
b = sqrt(interval(2))*a - im*sqrt(μ)/(2*sqrt(interval(2))*a)

############### Constants for μ = 0.28 ##########
#construction of the constant C0 (hat) described in Section 3.5.1
C0 = test_bessel_inequality(a,b,μ)

fourierE = CosFourier(4*N, π/d)⊗CosFourier(4*N, π/d)

# Construction of the Fourier series of the exponential functions in Lemma 3.6
# Actually we redefine E1 as E1/(exp(2ad)), E2 as E2/(exp(2ad)) and E12 as E12/(exp(4ad)) to improve precision in the proof
E1 = Sequence(fourierE ,interval.(big.(zeros((4*N+1)^2))))
E2 = Sequence(fourierE ,interval.(big.(zeros((4*N+1)^2))))
E12 = Sequence(fourierE ,interval.(big.(zeros((4*N+1)^2))))
m1 =4;m2=4;m3=4;

for n1 = 1:m1*N
    for n2 = 1:m1*N
        if (n1 <= m3*N)&&(n2 <= m3*N)
            E12[(n1,n2)] = real((-1)^n1*( (1-exp(-4*abig*dbig))/(2*abig-im*n1*π/dbig) ))*real((-1)^n2*( (1-exp(-4*abig*dbig))/(2*abig-im*n2*π/dbig) ))  
        end 
    end
    E1[(n1,0)] = real((-1)^n1*( (1-exp(-4*abig*dbig))/(2*abig-im*n1*π/dbig) ))*(1-exp(-4*abig*dbig))
    
    if (n1 <= m3*N)
    E12[(n1,0)] = real((-1)^n1*( (1-exp(-4*abig*dbig))/(2*abig-im*n1*π/dbig) ))*(1-exp(-4*abig*dbig))/(2*abig) 
    end
end

for n2 = 1:m1*N
    E2[(0,n2)] = real((-1)^n2*( (1-exp(-4*abig*dbig))/(2*abig-im*n2*π/dbig) ))*(1-exp(-4*abig*dbig))
    
    if (n2 <= m3*N)
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

D12 = convert(Vector{Interval{Float64}},exp2cos(2*N))  # diagonal matrix with the square root of the coefficients alpha_n
V = project((E1+E2)*V0,CosFourier(2*N, π/d)⊗CosFourier(2*N, π/d))
nV = coefficients(D12.*V0)'*coefficients(D12.*V) ; 

Zu1 = sqrt(C0^2/a^2*4*d^2*nV)

# ########################### Zu2 bound #####################################################

C1, C12, C2 = computation_b(a,d)  #computation of the constants C1, C12, C2 in Lemma 3.6
b_normal = interval(sup(norm_B^2*4*d^2))  #we introduce a normalization factor (which is a numerical trick to gain precision)
V = project((exp(-2*abig*dbig)*C1*E1 + C12*E12 + exp(-2*abig*dbig)C2*E2)*(b_normal*V0),CosFourier(2*N, π/d)⊗CosFourier(2*N, π/d))
nV2 = maximum([0 coefficients(D12.*V0)'*coefficients(D12.*V)]) ; D12 = Nothing
 
Zu2 = sqrt(Zu1^2+nV2*4*d^2/b_normal)

################################ Last component of  𝒵1 : ||B(π_NV0)||_2

n_inf = norm(V0_full -V0,1)*norm_B ; V0_full = Nothing

################################### Computation of  𝒵1 (Lemma 3.4)

𝒵1 = n_A + 2*norm_B*sqrt(Zu1^2 + Zu2^2) + n_inf
# note the factor 2 in front of norm_B. This comes from Theorem 3.7

################################### Z2 BOUND (Lemma 3.3) #####################################################

# We compute V2 = V*V where V is defined in Lemma 3.3
V2 = 2*(2*v1)*(6*v2*U) + 36*v2^2*U2
V2[(0,0)] += 4*v1^2  

GV02 = project(Multiplication(V),fourier0, fourier0, Interval{Float64}) 
GV02 = D1.*Bconj*GV02*B.*D2'
 # computation of the term \|V(B^N)^*\| using that 
# \|V(B^N)^*\|^2 = \|B^N(V*V)(B^N)^*\| = \|B^NV2(B^N)^*\|, using the adjoint
norm_VB = sqrt( opnorm(LinearOperator(coefficients(GV02)),2) ) ;
GV02 = Nothing
r0 = interval(3e-4)

# Using kappa from Lemma 3.8
κ = sqrt( (2*sqrt(μ) + (1+ μ)*(2*π - 2*atan(sqrt(μ))))/(8*(1+μ)*μ^(3/2)) + 2*π^2/d*(interval(3)^(3/4)/(μ^(7/4)) + interval(3)/(μ^(5/2))) )

# This provides the two components a0 and a1 of Z2(r) = a1r + a0 in Lemma 3.3
Z2 = [3*abs(v2)*norm_B*κ^2 ; κ*maximum([2*abs(v1) sqrt(norm_VB^2+norm(6*v2*U+2*v1,1)^2)])]

###################################### Y0 BOUND (Lemma 3.2)#######################################################

V = v1*U2 + v2*U2*U;

∂1 = project(Derivative((2,0)), fourier, fourier,Interval{Float64})
∂2 = project(Derivative((0,2)), fourier, fourier,Interval{Float64})
Δ = copy(∂1) ;  radd!(Δ,∂2) ; 

∂1 = Nothing ; ∂2 = Nothing

L = LinearOperator(fourier,fourier,Diagonal((diag(coefficients(Δ+I))).^2)) ;  Δ = Nothing
L = convert(Vector{Interval{Float64}},diag(coefficients(μ*I + L)))

 Y0 = 2*d*sup(sqrt( norm((B*(project(L.*U+V,fourier0))),2)^2 +  norm(L.*(U-project(U,fourier0)) + V - project(V,fourier0),2)^2 ) )
 display("Y0 has a value");display(Y0)


######### Choice for r0 and definition of Z2
Z2test = sup(Z2[1])*1e-4 + sup(Z2[2])
 rmin=(1-sup(𝒵1) - sqrt((1-sup(𝒵1))^2-2*sup(Y0)*Z2test))/Z2test
 r0 = interval(rmin) # obtain a candidate value for r0
 Z2 = Z2[1]*r0 + Z2[2]  # we can compute Z2(r0) once r0 is fixed

####################### Validation of the proof  #############################
# We verify the condition of Theorem 3.7
if 𝒵1 + Z2*r0 < 1
  if 1/2*Z2*r0^2 - (1-𝒵1) + Y0 < 0
    display("The proof was successful for r0 = ")
    display(sup(r0))  
    display("There exists a branch of periodic solutions converging to a localized pattern as the period tends to infinity.")
  else
    display("failure: discriminant is negative")
  end
else
    display("failure: linear term is positive")
end
