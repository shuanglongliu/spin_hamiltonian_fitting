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
# xs are the coefficients for the dot and cross products of spins

# There are n_conf-1 ys and 4*n_pair xs

# The first n_pair xs are
# alpha_pair1, alpha_pair2, ..., alpha_pairn

# The rest 3*n_pair xs are
# pi_pair1_x, pi_pair1_y, pi_pair1_z
# pi_pair2_x, pi_pair2_y, pi_pair2_z
# ...
# pi_pairn_x, pi_pairn_y, pi_pairn_z


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
    # Return spins, n_site, pairs, n_pair, eij, c_alpha, n_c_alpha, c_pi, n_c_pi


    data = YAML.load_file("input.yaml")
    
    spins = data["spins"]
    n_site = length(data["spins"])
    
    pairs = data["pairs"]
    n_pair = length(data["pairs"])
    
    positions = data["positions"]
    eij = get_eij(positions, pairs, n_pair)
    
    c_alpha = data["constraints_alpha"]
    n_c_alpha = length(c_alpha)
    c_pi = data["constraints_pi"]
    n_c_pi = length(c_pi)
    n_c = n_c_alpha + 3*n_c_pi
    
    return spins, n_site, pairs, n_pair, eij, c_alpha, n_c_alpha, c_pi, n_c_pi, n_c
end

function read_dipoles(n_site::Int64,spins::Vector{Float64})
    # Read input.yaml
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

function get_A(S::Vector{Vector{Vector{Float64}}}, eij::Vector{Vector{Float64}}, n_p::Int64, n_x::Int64, n_y::Int64, n_conf::Int64, i_conf_ref::Int64, n_pair::Int64, pairs::Vector{Vector{Int64}})

    # First, get the B matrix which constains the coefficients before subtraction
    # between any configuration and the reference configuration. 

    B = Array{Float64}(undef, n_p, n_x); B .= 0.0
    for i_conf=1:n_conf
        for i_cart = 1:3
            for i_pair=1:n_pair
                i, j = pairs[i_pair]
                B[3*(i_conf-1)+i_cart, i_pair] = S[i_conf][i][i_cart] * dot(S[i_conf][j], eij[i_pair]) - dot(S[i_conf][i], eij[i_pair]) * S[i_conf][j][i_cart]
                B[3*(i_conf-1)+i_cart, n_pair + 3*(i_pair-1) + i_cart] = dot(S[i_conf][i], S[i_conf][j])
            end
        end
    end
    
    # Then, get the A matrix 
    
    A = Array{Float64}(undef, n_y, n_x); A .= 0.0
    i = 1
    for i_conf=1:n_conf
        if i_conf != i_conf_ref
            for i_cart = 1:3
                A[3*(i-1)+i_cart, :] = B[3*(i_conf-1)+i_cart, :] - B[3*(i_conf_ref-1)+i_cart, :]
            end
            i += 1
        end
    end
    
    return A
end

function get_constraints(c_alpha::Vector{Vector{Int64}}, c_pi::Vector{Vector{Real}}, n_c_alpha::Int64, n_c_pi::Int64, n_c::Int64, n_x::Int64, n_pair::Int64)

    # Dimensions of the constraint matrix C: n_c by n_x
    C = Array{Float64}(undef, n_c, n_x); C .= 0.0
    
    for i_c_alpha=1:n_c_alpha
        # @show c_alpha[i_c_alpha]
        i_pair, j_pair = c_alpha[i_c_alpha]
        C[i_c_alpha, i_pair] = 1.0; C[i_c_alpha, j_pair] = -1.0
    end
    
    i_c = n_c_alpha + 1
    for i_c_pi=1:n_c_pi
        # @show c_pi[i_c_pi]
        i_pair, j_pair, alpha, beta, gamma = c_pi[i_c_pi]
        rot = RotZ(deg2rad(gamma)) * RotY(deg2rad(beta)) * RotZ(deg2rad(alpha))
        # @show rot
    
        C[i_c, n_pair + 3*(i_pair) - 2: n_pair + 3*(i_pair)] = rot[1, 1:3]
        C[i_c, n_pair + 3*(j_pair) - 2] = -1.0
        i_c += 1
        
        C[i_c, n_pair + 3*(i_pair) - 2: n_pair + 3*(i_pair)] = rot[2, 1:3]
        C[i_c, n_pair + 3*(j_pair) - 1] = -1.0
        i_c += 1
        
        C[i_c, n_pair + 3*(i_pair) - 2: n_pair + 3*(i_pair)] = rot[3, 1:3]
        C[i_c, n_pair + 3*(j_pair)] = -1.0
        i_c += 1
    end
    
    z = Vector{Float64}(undef, n_c); z .= 0.0
    
    return C, z
end

function print_results(x::Vector{Float64}, n_pair::Int64, pairs::Vector{Vector{Int64}})
    f = open("log", "w")
    for i_pair = 1:n_pair
        @printf f "Pair %3d %3d\n\n" pairs[i_pair][1] pairs[i_pair][2]
        @printf f "alpha = %12.6f e Ang\n" x[i_pair]
        @printf f "pi = (%12.6f %12.6f %12.6f) e Ang\n" x[n_pair + 3*i_pair-2] x[n_pair + 3*i_pair-1] x[n_pair + 3*i_pair]
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
    
    f = open("log", "a")
    
    @printf f "\n\nMaximal deviation: %12.6f e Ang.\n" dy_max
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
    # cd("/Users/shlufl/Documents/Researches_iMac/Projects/SLL-Mn4Na/enotes/polarization")
    
    spins, n_site, pairs, n_pair, eij, c_alpha, n_c_alpha, c_pi, n_c_pi, n_c = read_yaml()
    
    n_conf, S, P = read_dipoles(n_site,spins)
    
    # Number of components of polarization for all spin configurations
    n_p = 3*n_conf
    
    # Dimensions of the A matrix: n_y by n_x
    n_y = 3*(n_conf-1)
    n_x = 4*n_pair
    
    # Index of the reference configuration
    i_conf_ref = 1
        
    y = get_y(P, n_y, n_conf, i_conf_ref)
    
    A = get_A(S, eij, n_p, n_x, n_y, n_conf, i_conf_ref, n_pair, pairs)
    
    C, z = get_constraints(c_alpha, c_pi, n_c_alpha, n_c_pi, n_c, n_x, n_pair)
        
    kkt_sol = [2*A'*A   C'; C  zeros(n_c,n_c)] \ [2*A'*y; z]
    x = kkt_sol[1:n_x]
    
    print_results(x, n_pair, pairs)
    
    analysis(x::Vector{Float64}, y::Vector{Float64}, A::Array{Float64}, n_y::Int64)
end

main()

