using DelimitedFiles
using LinearAlgebra
using Rotations
using YAML
using Plots
using LaTeXStrings
using Statistics
using Printf

# Task: Minimize |A x - y|^2

# ys are changes in the dipole
# xs are irr coupling matrix elements

# \vec{p} = sum_{ij} S_\alpha \vec{P}^{ij}_{\alpha \beta} S_\beta

# For a irr pair ij
# p_\gamma = sum_{\alpha \beta} S_{i,\alpha} P^{ij}_{\gamma \alpha \beta} S_{j,\beta}

# For a non-irr pair ab <- ij, R transforms the local reference frome for the pair ij to the local reference frame for the pair ab
# p_\gamma = sum_{\alpha \beta gamma' alpha' beta'} R_{alpha alpha'} R_{beta beta'} R_{gamma gamma'} S_{a,\alpha} S_{b,\beta} P^{ij}_{\gamma' \alpha' \beta'}

# There are 3*(n_conf-1) ys and 27*n_irrp xs

# P^{ij}_{\gamma \alpha \beta} are arranged as follows
# \beta (spin direction) is the fastest index
# \alpha (spin direction) is the second fastest index
# \gamma (dipole direction) is the third fastest (composite) index
# ij (pair) is slowest index

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
    # Return spins, n_site, pairs, n_pair, eij, transformations, n_tranfo, n_c

    data = YAML.load_file("input.yaml")
    
    spins = data["spins"]
    n_site = length(data["spins"])
    
    pairs = data["pairs"]
    n_pair = length(data["pairs"])

    positions = data["positions"]
    eij = get_eij(positions, pairs, n_pair)

    irr_pairs = data["irr_pairs"]
    n_irrp = length(data["irr_pairs"])

    dict_irrp = Dict()
    for i in 1:n_irrp
        dict_irrp[irr_pairs[i]] = i
    end

    # rmats[i] = emats[i]'
    emats = Vector{Matrix{Float64}}(undef, 2)
    rmats = Vector{Matrix{Float64}}(undef, 2)

    emats_vec = data["emats"]
    for i=1:n_irrp
        emat_vec = emats_vec[i]
        rmat = reshape( emat_vec, 3, 3 )
        rmats[i] = rmat
        emat = permutedims(rmat, (2, 1))
        emats[i] = emat
    end

    transformations = data["transformations"]
    if transformations == nothing
        n_tranfo = 0
    else
        n_tranfo = length(transformations)
    end

    dict_tranfo = Dict()
    for i in 1:n_tranfo
        dict_tranfo[transformations[i][2]] = i
    end

    return spins, n_site, pairs, n_pair, eij, irr_pairs, n_irrp, dict_irrp, emats, rmats, transformations, n_tranfo, dict_tranfo
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

function get_B(S::Vector{Vector{Vector{Float64}}}, n_p::Int64, n_x::Int64, n_conf::Int64, n_pair::Int64, pairs::Vector{Vector{Int64}}, irr_pairs::Vector{Int64}, dict_irrp::Dict{Any, Any}, transformations::Vector{Vector{Real}}, dict_tranfo::Dict{Any, Any})
    B = Array{Float64}(undef, n_p, n_x); B .= 0.0
    for i_conf in 1:n_conf
        for gamma in 1:3
            irow = 3*(i_conf-1)+gamma
            for k_pair=1:n_pair
                i_site, j_site = pairs[k_pair]
                if k_pair in irr_pairs
                    # p_\gamma = sum_{\alpha \beta} S_{i,\alpha} S_{j,\beta} P^{ij}_{\gamma \alpha \beta}
                    i_irrp = dict_irrp[k_pair]
                    for alpha = 1:3
                        for beta = 1:3
                            icol = 27 * (i_irrp-1) + 9*(gamma-1) + 3*(alpha-1) + beta
                            B[irow, icol] = B[irow, icol] + S[i_conf][i_site][alpha] * S[i_conf][j_site][beta]
                        end
                    end
                else
                    # p_\gamma = sum_{\alpha \beta gamma' alpha' beta'} R_{alpha alpha'} R_{beta beta'} R_{gamma gamma'} S_{a,\alpha} S_{b,\beta} P^{ij}_{\gamma' \alpha' \beta'}
                    i_tranfo = dict_tranfo[k_pair]
                    i_pair, j_pair, alpha_euler, beta_euler, gamma_euler = transformations[i_tranfo] # j_pair = k_pair
                    rot = RotZ(deg2rad(gamma_euler)) * RotY(deg2rad(beta_euler)) * RotZ(deg2rad(alpha_euler))
                    i_irrp = dict_irrp[i_pair]
                    for gamma_ = 1:3
                        for alpha_ = 1:3
                            for beta_ = 1:3
                                icol = 27 * (i_irrp-1) + 9*(gamma_-1) + 3*(alpha_-1) + beta_
                                for alpha = 1:3
                                    for beta = 1:3
                                        B[irow, icol] = B[irow, icol] + rot[alpha, alpha_] * rot[beta, beta_] * rot[gamma, gamma_] * S[i_conf][i_site][alpha] * S[i_conf][j_site][beta]
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    return B
end

function get_A(S::Vector{Vector{Vector{Float64}}}, n_p::Int64, n_x::Int64, n_y::Int64, n_conf::Int64, i_conf_ref::Int64, n_pair::Int64, pairs::Vector{Vector{Int64}}, irr_pairs::Vector{Int64}, dict_irrp::Dict{Any, Any}, transformations::Vector{Vector{Real}}, dict_tranfo::Dict{Any, Any})

    # S: spins
    # n_p = 3*n_conf: number of components of polarization for all spin configurations
    # n_x = 3*9*n_pair: number of coupling matrix elements
    # n_y = 3*(n_conf - 1): number of change of electric dipole 
    # n_conf: number of spin configurations
    # i_conf_ref: index of the reference spin configuration
    # n_pair: number of spin pairs
    # pairs: spin pairs
    # irr_pairs: indices of irreducible/unique pairs
    # dict_irrp: a dictionary of indices of irreducible/unique pairs
    # transformations: a list/vector of transformations, see input.yaml for an explanation.
    # dict_tranfo: a dictionary of indices of transformations, see read_yaml for its definition.

    # First, get the B matrix which contains the coefficients before subtraction
    # between any configuration and the reference configuration. 

    B = get_B(S, n_p, n_x, n_conf, n_pair, pairs, irr_pairs, dict_irrp, transformations, dict_tranfo)

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

function check_fitting_quality(x::Vector{Float64}, y::Vector{Float64}, A::Array{Float64}, n_y::Int64)
    y_fit = A * x
    dy = y_fit - y
    dy = abs.(dy)
    dy_max = maximum(dy)
    dy_avg = Statistics.mean(dy)
    dy_std = Statistics.std(dy)

    dy_per = Vector{Float64}(undef, 1)
    dy_per[1] = 100*dy[1]/y[1]

    for i=2:n_y
        if y[i] >= 0.001
            push!( dy_per, 100*dy[i]/y[i] )
        end
    end

    dy_max_per = maximum(dy_per)
    dy_avg_per = Statistics.mean(dy_per)
    dy_std_per = Statistics.std(dy_per)
    
    f = open("fit.log", "w")
    
    @printf f "Maximal deviation: %12.6f e Ang.\n" dy_max
    @printf f "Average deviation: %12.6f e Ang.\n" dy_avg
    @printf f " std of deviation: %12.6f e Ang.\n" dy_std

    @printf f "\n"

    @printf f "Maximal deviation in percentage: %12.6f \n" dy_max_per
    @printf f "Average deviation in percentage: %12.6f \n" dy_avg_per
    @printf f " std of deviation in percentage: %12.6f \n" dy_std_per

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

function save_P(x::Vector{Float64}, pairs::Vector{Vector{Int64}}, irr_pairs::Vector{Int64}, n_irrp::Int64, fname::String)
    f = open(fname, "w")
    directions = ["x", "y", "z"]
    for i_irrp = 1:n_irrp
        i_site, j_site = pairs[irr_pairs[i_irrp]]
        @printf f "Pair %3d %3d\n\n" i_site j_site
        for gamma=1:3
            @printf f "The %s component of coupling matrix\n\n" directions[gamma]
            for alpha=1:3
                index = 27*(i_irrp-1) + 9*(gamma-1) + 3*(alpha-1)
                @printf f "%12.6f %12.6f %12.6f\n" x[index+1] x[index+2] x[index+3]
            end
            @printf f "\n"
        end
        @printf f "\n"
    end
    close(f)
end

function decompose_P(x::Vector{Float64}, n_x::Int64, n_irrp::Int64)
    # P^J = ave( diag( (P + P^T) / 2 ) )
    #     P^J_{gamma x x} = P^J_{gamma y y} = P^J_{gamma z z} = (P_{gamma x x} + P_{gamma y y} + P_{gamma z z})/3
    # P^D = (P - P^T) / 2
    #     P^D_{gamma alpha beta} = (P_{gamma alpha beta} - P_{gamma beta alpha}) / 2
    # P^Gamma = P - P^J - P^D

    x_J     = Array{Float64}(undef, n_x); x_J     .= 0.0
    x_D     = Array{Float64}(undef, n_x); x_D     .= 0.0
    x_Gamma = Array{Float64}(undef, n_x); x_Gamma .= 0.0

    for i_irrp = 1:n_irrp
        offset = 27*(i_irrp-1)

        # Obtain P^J
        for gamma = 1:3
            i = 9*(gamma-1) + 3*(1 - 1) + 1 + offset
            j = 9*(gamma-1) + 3*(2 - 1) + 2 + offset
            k = 9*(gamma-1) + 3*(3 - 1) + 3 + offset
            x_J[i] = (x[i] + x[j] + x[k])/3
            x_J[j] = (x[i] + x[j] + x[k])/3
            x_J[k] = (x[i] + x[j] + x[k])/3
        end

        # Obtain P^D
        for gamma = 1:3
            for alpha = 1:3
                for beta = 1:3
                    i = 9*(gamma-1) + 3*(alpha - 1) + beta + offset
                    j = 9*(gamma-1) + 3*(beta - 1) + alpha + offset
                    x_D[i] = (x[i] - x[j])/2
                end
            end
        end
    end

    x_Gamma = x - x_J - x_D

    return x_J, x_D, x_Gamma
end

function get_pi_and_d(x_J::Vector{Float64}, x_D::Vector{Float64}, pairs::Vector{Vector{Int64}}, irr_pairs::Vector{Int64}, n_irrp::Int64, fname::String)
    # S_{i alpha} P^J_{gamma alpha beta} S_{j beta} = pi_{gamma} S_i \cdot S_j
    #     pi_{gamma} = P^J_{gamma x x} = P^J_{gamma y y} = P^J_{gamma z z}

    # S_{i alpha} P^D_{gamma alpha beta} S_{j beta} = d_{gamma} \cdot (S_i \cross S_j)
    #     d_{gamma x} = P^D_{gamma y z}
    #     d_{gamma y} = P^D_{gamma z x}
    #     d_{gamma z} = P^D_{gamma x y}

    # S_{i alpha} P^Gamma_{gamma alpha beta} S_{j beta} :: No simplified expression

    f = open(fname, "w")

    pi_vec = Array{Float64}(undef, 3)
    d_vecs = Array{Float64}(undef, 3, 3)

    for i_irrp = 1:n_irrp
        i_site, j_site = pairs[irr_pairs[i_irrp]]
        @printf f "Pair %3d %3d\n\n" i_site j_site

        offset = 27*(i_irrp-1)

        for gamma = 1:3
            i = 9*(gamma-1) + 3*(1 - 1) + 1 + offset
            pi_vec[gamma] = x_J[i]
        end
        @printf f "pi = %12.6f %12.6f %12.6f\n\n" pi_vec[1] pi_vec[2] pi_vec[3]

        for gamma = 1:3
            i = 9*(gamma-1) + 3*(2 - 1) + 3 + offset
            j = 9*(gamma-1) + 3*(3 - 1) + 1 + offset
            k = 9*(gamma-1) + 3*(1 - 1) + 2 + offset
            d_vecs[gamma, 1] = x_D[i]
            d_vecs[gamma, 2] = x_D[j]
            d_vecs[gamma, 3] = x_D[k]
        end
        @printf f "dx = %12.6f %12.6f %12.6f\n" d_vecs[1, 1] d_vecs[1, 2] d_vecs[1, 3]
        @printf f "dy = %12.6f %12.6f %12.6f\n" d_vecs[2, 1] d_vecs[2, 2] d_vecs[2, 3]
        @printf f "dz = %12.6f %12.6f %12.6f\n" d_vecs[3, 1] d_vecs[3, 2] d_vecs[3, 3]
        @printf f "\n\n"
    end

    close(f)
end

function get_P_in_local_reference_frames(x::Vector{Float64}, n_x::Int64, n_irrp::Int64, rmats::Vector{Matrix{Float64}})
    x_loc = Vector{Float64}(undef, n_x); x_loc .= 0.0
    for i_irrp = 1:n_irrp
        for gamma = 1:3
            for alpha = 1:3
                for beta = 1:3
                    i_x_loc = 27*(i_irrp-1) + 9*(gamma-1) + 3*(alpha-1) + beta
                    for gamma_ = 1:3
                        for alpha_ = 1:3
                            for beta_ = 1:3
                                i_x = 27*(i_irrp-1) + 9*(gamma_-1) + 3*(alpha_-1) + beta_
                                x_loc[i_x_loc] = x_loc[i_x_loc] + rmats[i_irrp][gamma, gamma_] * rmats[i_irrp][alpha, alpha_] * rmats[i_irrp][beta, beta_] * x[i_x]
                            end
                        end
                    end
                end
            end
        end
    end
    return x_loc
end

function save_results(x::Vector{Float64}, n_x::Int64, pairs::Vector{Vector{Int64}}, irr_pairs::Vector{Int64}, n_irrp::Int64, modifier::String)
    save_P(x, pairs, irr_pairs, n_irrp, "P_total" * modifier * ".log")

    x_J, x_D, x_Gamma = decompose_P(x, n_x, n_irrp)

    save_P(x_J, pairs, irr_pairs, n_irrp, "P_J" * modifier * ".log")
    save_P(x_D, pairs, irr_pairs, n_irrp, "P_D" * modifier * ".log")
    save_P(x_Gamma, pairs, irr_pairs, n_irrp, "P_Gamma" * modifier * ".log")

    get_pi_and_d(x_J, x_D, pairs, irr_pairs, n_irrp, "pi_and_d" * modifier * ".log")
end

function print_relative_contribution(y_per::Vector{Float64}, f::IOStream, tag::String)
    y_max_per = maximum(filter(!isnan,y_per))
    y_avg_per = NaNStatistics.nanmean(y_per)
    y_std_per = NaNStatistics.nanstd(y_per)

    @printf f "Maximal contribution from %8s: %6.1f percent.\n" tag  y_max_per
    @printf f "Average contribution from %8s: %6.1f percent.\n" tag  y_avg_per
    @printf f " std of contribution from %8s: %6.1f percent.\n" tag  y_std_per

    @printf f "\n"
end

function get_relative_contribution(x::Vector{Float64}, n_x::Int64, n_y::Int64, A::Array{Float64}, pairs::Vector{Vector{Int64}}, n_irrp::Int64, irr_pairs::Vector{Int64})

    # Get relative contributions from irreducible pairs
    y_irr = Vector{Vector{Float64}}(undef, n_irrp)
    y_irr_per = Vector{Float64}(undef, n_irrp)
    y_sum_irr = Vector{Float64}(undef, n_y); y_sum_irr .= 0.0
    for i_irrp = 1:n_irrp
        i = 27 * (i_irrp - 1) + 1
        j = 27 * i_irrp

        y_irr[i_irrp] = abs.( A[:, i:j] * x[i:j] )
        y_sum_irr = y_sum_irr .+ y_irr[i_irrp]
    end

    x_J, x_D, x_Gamma = decompose_P(x, n_x, n_irrp)

    f = open("relative_contribution.log", "w")
    
    for i_irrp = 1:n_irrp
        i_site, j_site = pairs[irr_pairs[i_irrp]]
        y_irr_per[i_irrp] = NaNStatistics.nanmean( 100 .* y_irr[i_irrp] ./ y_sum_irr )

        @printf f "Pair %3d %3d: relative contribution = %6.1f percent.\n\n" i_site j_site y_irr_per[i_irrp]

        i = 27 * (i_irrp - 1) + 1
        j = 27 * i_irrp

        y_J     = abs.( A[:, i:j] * x_J[i:j] )
        y_D     = abs.( A[:, i:j] * x_D[i:j] )
        y_Gamma = abs.( A[:, i:j] * x_Gamma[i:j] )
    
        y_sum   = y_J .+ y_D .+ y_Gamma

        y_J_per     = 100 .* ( y_J     ./ y_sum )
        y_D_per     = 100 .* ( y_D     ./ y_sum )
        y_Gamma_per = 100 .* ( y_Gamma ./ y_sum )
    
        print_relative_contribution(y_J_per, f, "P_J")
        print_relative_contribution(y_D_per, f, "P_D")
        print_relative_contribution(y_Gamma_per, f, "P_Gamma")
    end

    close(f)
end

function save_P_for_one_pair_for_latex(f::IOStream, x::Vector{Float64}, i_irrp::Int64)
    x = 10^5 .* x
    directions = ["x", "y", "z"]
    for gamma=1:3
        @printf f "The %s component of coupling matrix\n\n" directions[gamma]
        for alpha=1:3
            index = 27*(i_irrp-1) + 9*(gamma-1) + 3*(alpha-1)
            if alpha != 3
                @printf f "%6.1f & %6.1f &  %6.1f \\\\\n" x[index+1] x[index+2] x[index+3]
            else
                @printf f "%6.1f & %6.1f &  %6.1f\n" x[index+1] x[index+2] x[index+3]
            end
        end
        @printf f "\n"
    end
end

function save_P_for_latex(x::Vector{Float64}, n_x::Int64, pairs::Vector{Vector{Int64}}, irr_pairs::Vector{Int64}, n_irrp::Int64)
    x_J, x_D, x_Gamma = decompose_P(x, n_x, n_irrp)

    f = open("P_for_latex.log", "w")
    for i_irrp = 1:n_irrp
        i_site, j_site = pairs[irr_pairs[i_irrp]]
        @printf f "Pair %3d %3d\n\n" i_site j_site
        @printf f "Total \n\n"; save_P_for_one_pair_for_latex(f, x, i_irrp); @printf f "\n"
        @printf f "The isotropic term\n\n"; save_P_for_one_pair_for_latex(f, x_J, i_irrp); @printf f "\n"
        @printf f "The antisymmetric term\n\n"; save_P_for_one_pair_for_latex(f, x_D, i_irrp); @printf f "\n"
        @printf f "The anisotropic symmetric term\n\n"; save_P_for_one_pair_for_latex(f, x_Gamma, i_irrp); @printf f "\n"
    end
    close(f)
end

function main()
    spins, n_site, pairs, n_pair, eij, irr_pairs, n_irrp, dict_irrp, emats, rmats, transformations, n_tranfo, dict_tranfo = read_yaml()

    n_conf, S, P = read_dipoles(n_site,spins)

    # Number of components of polarization for all spin configurations
    n_p = 3*n_conf
    
    # Dimensions of the A matrix: n_y by n_x
    n_y = 3*(n_conf-1)
    n_x = 27*n_irrp
    
    # Index of the reference configuration
    i_conf_ref = 1

    y = get_y(P, n_y, n_conf, i_conf_ref)

    A = get_A(S, n_p, n_x, n_y, n_conf, i_conf_ref, n_pair, pairs, irr_pairs, dict_irrp, transformations, dict_tranfo)

    x = A \ y

    #check_fitting_quality(x, y, A, n_y)

    #save_results(x, n_x, pairs, irr_pairs, n_irrp, "")

    x_loc = get_P_in_local_reference_frames(x, n_x, n_irrp, rmats)

    #save_results(x_loc, n_x, pairs, irr_pairs, n_irrp, "_loc")

    #get_relative_contribution(x, n_x, n_y, A, pairs, n_irrp, irr_pairs)

    save_P_for_latex(x_loc, n_x, pairs, irr_pairs, n_irrp)
end

main()


