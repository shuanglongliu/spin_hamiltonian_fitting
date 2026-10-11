using DelimitedFiles
using LinearAlgebra
using Rotations
using YAML
using Plots
using LaTeXStrings
using Statistics
using Printf

# Task: Minimize |A x - y|^2 subject to C x = z. 

# ys are changes in the dipole
# xs are coupling matrix elements

# \vec{p} = sum_{ij} S_\alpha \vec{P}^{ij}_{\alpha \beta} S_\beta or
# p_\gamma = sum_{ij \alpha \beta} S_{i,\alpha} P^{ij}_{\gamma \alpha \beta} S_{j,\beta}
# p_\gamma' = sum_{ij \gamma \alpha \beta} S_{i,\alpha} S_{j,\beta} \delta_{\gamma, \gamma'} P^{ij}_{\gamma \alpha \beta}

# There are 3*(n_conf-1) ys and 3*9*n_pair xs
# The factor 3 is the number of components of a vector, i.e. x, y, and z.

# P^{ij}_{\gamma \alpha \beta} are arranged as follows
# \beta (spin direction) is the fastest index
# \alpha (spin direction) is the second fastest index
# \gamma (dipole direction) is the third fastest (composite) index
# ij (pair) is slowest index

function cls_solve_kkt(A,b,C,d)
    m, n = size(A)
    p, n = size(C)
    G = A'*A  # Gram matrix
    KKT = [2*G C'; C zeros(p,p)]  # KKT matrix
    xzhat = KKT \ [2*A'*b; d]
    return xzhat[1:n,:]
end

function cls_solve(A,b,C,d)
    m, n = size(A)
    p, n = size(C)
    Q, R = qr([A; C])
    Q = Matrix(Q)
    Q1 = Q[1:m,:]
    Q2 = Q[m+1:m+p,:]
    Qtil, Rtil = qr(Q2')
    Qtil = Matrix(Qtil)
    w = Rtil \ (2*Qtil'*Q1'*b - 2*(Rtil'\d))
    return xhat = R \ (Q1'*b - Q2'*w/2)
end

function kronecker_delta(i::Integer, j::Integer)
    k = 0
    if i == j
        k = 1
    end
    return k
end

function sph2cart(r::Float64, theta::Float64, phi::Float64, deg::Bool)
    if deg 
        theta = deg2rad(theta)
        phi = deg2rad(phi)
    end 
    x = r * sin(theta) * cos(phi)
    y = r * sin(theta) * sin(phi)
    z = r * cos(theta)
    return [x, y, z]
end

function get_eij(positions::Vector{Vector{Float64}}, pairs::Vector{Vector{Int64}}, n_pair::Int64)
    eij = [Vector{Float64}(undef, 3) for _ in 1:n_pair]
    for i_pair=1:n_pair
        i, j = pairs[i_pair]
        v = positions[j]-positions[i]
        eij[i_pair] = v / norm(v)
    end
    return eij
end

function read_yaml()
    
    # Read input.yaml
    # Return spins, n_site, pairs, n_pair, eij, constraints, n_c_matrix, n_c

    data = YAML.load_file("input.yaml")
    
    spins = data["spins"]
    n_site = length(data["spins"])
    
    pairs = data["pairs"]
    n_pair = length(data["pairs"])
    
    positions = data["positions"]
    eij = get_eij(positions, pairs, n_pair)
    
    constraints = data["constraints"]
    if constraints == nothing
        n_c_matrix = 0
    else
        n_c_matrix = length(constraints)
    end
    n_c = 3 * 9 * n_c_matrix

    return spins, n_site, pairs, n_pair, eij, constraints, n_c_matrix, n_c
end

function read_dipoles(n_site::Int64,spins::Vector{Float64})
    # Return n_conf, S (spin configurations), and P (dipoles)

    data = readdlm("dipoles.dat")
    n_conf, _ = size(data)

    S = [[Vector{Float64}(undef,3) for _ in 1:n_site] for _ in 1:n_conf]
    P = [Vector{Float64}(undef,3) for _ in 1:n_conf]
    
    # Get spin configurations and dipoles
    for i = 1:n_conf
        for j = 1:n_site
            theta = data[i, j*2 - 1]
            phi = data[i, j*2]
            S[i][j] = sph2cart(spins[j], theta, phi, true)
        end
        P[i] = data[i, 2*n_site+1:2*n_site+3]
    end
    
    return n_conf, S, P
end

function get_y(P::Vector{Vector{Float64}}, n_y::Integer, n_conf::Integer, i_conf_ref::Integer)
    # Get changes in dipole, i.e. the vector y.
    # n_y = 3*(n_conf - 1): number of change of electric dipole 

    # y = y(i_conf, gamma)
    # (i_conf, gamma) is a composite index
    # \gamma is faster than i_conf
    
    y = Vector{Float64}(undef, n_y); y .= 0.0
    i = 1
    for i_conf = 1:n_conf
        if i_conf != i_conf_ref
            y[3*i-2: 3*i] = P[i_conf] - P[i_conf_ref]
            i += 1
        end
    end
    
    return y
end

function get_A(S::Vector{Vector{Vector{Float64}}}, n_p::Int64, n_x::Int64, n_y::Int64, n_conf::Int64, i_conf_ref::Int64, n_pair::Int64, pairs::Vector{Vector{Int64}})

    # S: spins
    # n_p = 3*n_conf: number of components of polarization for all spin configurations
    # n_x = 3*9*n_pair: number of coupling matrix elements
    # n_y = 3*(n_conf - 1): number of change of electric dipole 
    # n_conf: number of spin configurations
    # i_conf_ref: index of the reference spin configuration
    # n_pair: number of spin pairs
    # pairs: spin pairs

    # First, get the B matrix which contains the coefficients before subtraction
    # between any configuration and the reference configuration. 

    # p_\gamma = sum_{i_pair \alpha' \beta'} S_{i,\alpha'} S_{j,\beta'} P^{i_pair}_{\gamma \alpha' \beta'}

    B = Array{Float64}(undef, n_p, n_x); B .= 0.0
    for i_conf=1:n_conf
        for gamma = 1:3
            irow = 3*(i_conf-1)+gamma
            for i_pair=1:n_pair
                i_site, j_site = pairs[i_pair]
                for alpha_ = 1:3
                    for beta_ = 1:3
                        icol = 27 * (i_pair-1) + 9*(gamma-1) + 3*(alpha_-1) + beta_
                        B[irow, icol] = S[i_conf][i_site][alpha_] * S[i_conf][j_site][beta_]
                    end 
                end
            end
        end
    end

    # Then, get the A matrix 
    # A[row] = B[row] - B[ref row]
    
    A = Array{Float64}(undef, n_y, n_x); A .= 0.0
    i = 1
    for i_conf=1:n_conf
        if i_conf != i_conf_ref
            for gamma = 1:3
                A[3*(i-1)+gamma, :] = B[3*(i_conf-1)+gamma, :] - B[3*(i_conf_ref-1)+gamma, :]
            end
            i += 1
        end
    end

    return A
end

function get_constraints(constraints::Vector{Vector{Real}}, n_c_matrix::Int64, n_c::Int64, n_x::Int64)

    # P^{i_pair}: Si, Sj -> p
    # P^{j_pair}: R(Si) R(Sj) -> R(p)
    # ( R Si ) P^{j_pair} ( R Sj ) = R ( Si P^{i_pair} Sj )
    # R_{alpha alpha'} S_{i, alpha'} P^{j_pair}_{gamma alpha beta} R_{beta beta'} S_{j beta'} = R_{gamma gamma'} S_{i, alpha'} P^{i_pair}_{gamma' alpha' beta'} S_{j beta'}
    # For any gamma, sum_{alpha' beta'} ( sum_{alpha beta} { R_{alpha alpha'} P^{j_pair}_{gamma alpha beta} R_{beta beta'} } - sum_{gamma'} { R_{gamma gamma'} P^{i_pair}_{gamma' alpha' beta'}} ) S_{i, alpha'} S_{j beta'} = 0
    # The above equation is valid for any S_{i, alpha'} and S_{j beta'}, so
    # R_{alpha alpha'} P^{j_pair}_{gamma alpha beta} R_{beta beta'} - R_{gamma gamma'} P^{i_pair}_{gamma' alpha' beta'} = 0
    # The sum is over alpha, beta and gamma'. The indices gamma, alpha', and beta' are fixed.
    # Rewrite the above equation as
    # R_{alpha alpha'} R_{beta beta'} P^{j_pair}_{gamma alpha beta} - R_{gamma gamma'} P^{i_pair}_{gamma' alpha' beta'} = 0
    # This is equivalent to the following
    # P^{j_pair}_{gamma alpha beta} = R_{gamma gamma'} R_{alpha alpha'} R_{beta beta'} P^{i_pair}_{gamma' alpha' beta'} or
    # R_{gamma gamma'} R_{alpha alpha'} R_{beta beta'} P^{i_pair}_{gamma' alpha' beta'} - P^{j_pair}_{gamma alpha beta} = 0

    # Dimensions of the constraint matrix C: n_c by n_x
    C = Array{Float64}(undef, n_c, n_x); C .= 0.0
    
    for i_c_matrix=1:n_c_matrix
        # alpha, beta, gamma are Euler angles below
        i_pair, j_pair, alpha, beta, gamma = constraints[i_c_matrix]
        rot = RotZ(deg2rad(gamma)) * RotY(deg2rad(beta)) * RotZ(deg2rad(alpha))

        # alpha, beta, gamma are directions x/y/z below
        for gamma in 1:3
            for alpha in 1:3
                for beta in 1:3
                    irow = 27*(i_c_matrix-1) + 9*(gamma-1) + 3*(alpha-1) + beta

                    for gamma_ in 1:3
                        for alpha_ in 1:3
                            for beta_ in 1:3
                                icol = 27 * (i_pair-1) + 9*(gamma_-1) + 3*(alpha_-1) + beta_
                                C[irow, icol] = rot[gamma, gamma_] * rot[alpha, alpha_] * rot[beta, beta_]
                            end
                        end
                    end

                    icol = 27 * (j_pair-1) + 9*(gamma-1) + 3*(alpha-1) + beta
                    C[irow, icol] = -1.0
                end
            end
        end
    end
    
    z = Vector{Float64}(undef, n_c); z .= 0.0
    
    return C, z
end

function print_results(x::Vector{Float64}, n_pair::Int64, pairs::Vector{Vector{Int64}})
    f = open("P_total.log", "w")
    directions = ["x", "y", "z"]
    for i_pair = 1:n_pair
        @printf f "Pair %3d %3d\n\n" pairs[i_pair][1] pairs[i_pair][2]
        for gamma=1:3
            @printf f "The %s component of coupling matrix\n\n" directions[gamma]
            for alpha=1:3
                index = 27*(i_pair-1) + 9*(gamma-1) + 3*(alpha-1)
                @printf f "%12.6f %12.6f %12.6f\n" x[index+1] x[index+2] x[index+3]
            end
            @printf f "\n"
        end
        @printf f "\n"
    end
    close(f)
end

function print_results(x::Matrix{Float64}, n_pair::Int64, pairs::Vector{Vector{Int64}})
    f = open("P_total.log", "w")
    directions = ["x", "y", "z"]
    for i_pair = 1:n_pair
        @printf f "Pair %3d %3d\n\n" pairs[i_pair][1] pairs[i_pair][2]
        for gamma=1:3
            @printf f "The %s component of coupling matrix\n\n" directions[gamma]
            for alpha=1:3
                index = 27*(i_pair-1) + 9*(gamma-1) + 3*(alpha-1)
                @printf f "%12.6f %12.6f %12.6f\n" x[index+1] x[index+2] x[index+3]
            end
            @printf f "\n"
        end
        @printf f "\n"
    end
    close(f)
end

function analysis(x::Vector{Float64}, y::Vector{Float64}, A::Array{Float64}, n_y::Int64)
    y_fit = A * x
    dy = y_fit - y
    dy = abs.(dy)
    dy_max = maximum(dy)
    dy_avg = Statistics.mean(dy)
    
    f = open("fit.log", "w")
    
    @printf f "Maximal deviation: %12.6f e Ang.\n" dy_max
    @printf f "Average deviation: %12.6f e Ang.\n" dy_avg
    
    close(f)
    
    Plots.gr()
    
    x_ = range(1, n_y)
    p = plot(
        x_, 
        [y, y_fit], 
        seriestype=:scatter,
        label=["DFT" "Fit"], 
        ylabel = L"\delta P_\alpha\quad (e\;\AA)",
        xlabel = L"\mathrm{Index\; of\;} \delta P_\alpha\; \mathrm{for\; all\; spin\; configurations}", 
        foreground_color_legend = nothing, 
        legend=:bottomright, 
    )
    
    savefig(p, "fitting_quality.pdf")
end

function analysis(x::Matrix{Float64}, y::Vector{Float64}, A::Array{Float64}, n_y::Int64)
    y_fit = A * x
    dy = y_fit - y
    dy = abs.(dy)
    dy_max = maximum(dy)
    dy_avg = Statistics.mean(dy)
    
    f = open("fit.log", "w")
    
    @printf f "Maximal deviation: %12.6f e Ang.\n" dy_max
    @printf f "Average deviation: %12.6f e Ang.\n" dy_avg
    
    close(f)
    
    Plots.gr()
    
    x_ = range(1, n_y)
    p = plot(
        x_, 
        [y, y_fit], 
        seriestype=:scatter,
        label=["DFT" "Fit"], 
        ylabel = L"\delta P_\alpha\quad (e\;\AA)",
        xlabel = L"\mathrm{Index\; of\;} \delta P_\alpha\; \mathrm{for\; all\; spin\; configurations}", 
        foreground_color_legend = nothing, 
        legend=:bottomright, 
    )
    
    savefig(p, "fitting_quality.pdf")
end

function main()
    spins, n_site, pairs, n_pair, eij, constraints, n_c_matrix, n_c = read_yaml()

    n_conf, S, P = read_dipoles(n_site,spins)

    # Number of components of polarization for all spin configurations
    n_p = 3*n_conf
    
    # Dimensions of the A matrix: n_y by n_x
    n_y = 3*(n_conf-1)
    n_x = 3*9*n_pair # number of coupling matrix elements
    
    # Index of the reference configuration
    i_conf_ref = 1
        
    y = get_y(P, n_y, n_conf, i_conf_ref)

    A = get_A(S, n_p, n_x, n_y, n_conf, i_conf_ref, n_pair, pairs)

    #x = A \ y

    if n_c == 0
        x = A \ y
    else
        C, z = get_constraints(constraints, n_c_matrix, n_c, n_x)

        x = cls_solve_kkt(A, y, C, z) # x is a matrix
        #x = cls_solve(A, y, C, z) # x is a vector
    end
    
    print_results(x, n_pair, pairs)
    
    analysis(x, y, A, n_y)
end

main()


