
using YAML
using LinearAlgebra
using Printf

# Base.show(io::IO, f::Float64) = @printf(io, "%8.3f", f)

function get_index_of_reference_frame(spin_pair::Vector, spin_pairs, n_spin_pair)
    i_reference_frame = 0
    for i=1:n_spin_pair
        if spin_pair[1] == spin_pairs[i][1] && spin_pair[2] == spin_pairs[i][2]
            i_reference_frame = i
            break
        end
    end
    return i_reference_frame
end

function transform_reference_frame(spin_vectors::Array, n_site::Int64, emat_T_from::Array{Float64, 2}, emat_T_to::Array{Float64, 2})
    spin_vectors_to = [Vector{Float64}(undef,3) for _ in 1:n_site]
    for i=1:n_site
        spin_vectors_to[i] = transpose(emat_T_to) * emat_T_from * spin_vectors[i]
    end
    return spin_vectors_to
end

function fatten_one_energy_pair(energy_pair::Dict, reference_frames::Array, n_site::Int64, spins::Array, spin_pairs::Array, n_spin_pair::Int64)
    energies = Dict()
    energies["conf1"] = energy_pair["conf1"][2*n_site + 1]
    energies["conf2"] = energy_pair["conf2"][2*n_site + 1]
    
    spin_vectors_conf1 = [Vector{Float64}(undef,3) for _ in 1:n_site]
    spin_vectors_conf2 = [Vector{Float64}(undef,3) for _ in 1:n_site]

    for i=1:n_site
        theta = energy_pair["conf1"][i*2 - 1]
        phi = energy_pair["conf1"][i*2]
        spin_vector = sph2cart(spins[i], theta, phi, true)
        spin_vectors_conf1[i] = spin_vector

        theta = energy_pair["conf2"][i*2 - 1]
        phi = energy_pair["conf2"][i*2]
        spin_vector = sph2cart(spins[i], theta, phi, true)
        spin_vectors_conf2[i] = spin_vector
    end
        
    configurations = Dict()
    configurations["conf1"] = [[Vector{Float64}(undef,3) for _ in 1:n_site] for _ in 1:n_spin_pair]
    configurations["conf2"] = [[Vector{Float64}(undef,3) for _ in 1:n_site] for _ in 1:n_spin_pair]
    
    this_frame = get_index_of_reference_frame(energy_pair["pair"], spin_pairs, n_spin_pair)
    for i_reference_frame=1:n_spin_pair
        spin_vectors = transform_reference_frame(spin_vectors_conf1, n_site, reference_frames[this_frame], reference_frames[i_reference_frame])
        for i=1:n_site
            configurations["conf1"][i_reference_frame][i] = spin_vectors[i]
        end
        spin_vectors = transform_reference_frame(spin_vectors_conf2, n_site, reference_frames[this_frame], reference_frames[i_reference_frame])
        for i=1:n_site
            configurations["conf2"][i_reference_frame][i] = spin_vectors[i]
        end
    end
    
    return (configurations, energies)
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

function get_one_row_of_coefficients(data::Dict, i::Int64, reference_frames::Array, n_site::Int64, spins::Array, spin_pairs::Array, n_spin_pair::Int64, factor::Float64)
    configurations, energies = fatten_one_energy_pair(data["energies"][i], reference_frames, n_site, spins, spin_pairs, n_spin_pair)
    n_J = 9*n_spin_pair
    Ci = Array{Float64, 1}(undef, n_J)
    ind = 1
    for i_spin_pair=1:n_spin_pair
        i = spin_pairs[i_spin_pair][1]
        j = spin_pairs[i_spin_pair][2]
        i_reference_frame = get_index_of_reference_frame([i, j], spin_pairs, n_spin_pair)
        for alpha=1:3
            for beta=1:3
                # Order of Js: pair-1_xx, pair-1_xy, ..., pair-n_zz
                y = configurations["conf2"][i_reference_frame][i][alpha] * configurations["conf2"][i_reference_frame][j][beta]
                x = configurations["conf1"][i_reference_frame][i][alpha] * configurations["conf1"][i_reference_frame][j][beta]
                Ci[ind] = factor*(y - x)
                ind += 1
            end
        end
    end
    return (Ci, energies)
end

function get_coefficient_matrix(data::Dict, reference_frames::Array, n_site::Int64, spins::Array, spin_pairs::Array, n_spin_pair::Int64, n_energy_pair::Int64, factor::Float64)
    n_J = 9*n_spin_pair
    Cs = Array{Float64, 2}(undef, n_energy_pair, n_J)
    Es = Array{Float64, 1}(undef, n_energy_pair)
    for i=1:n_energy_pair
        Ci, energies = get_one_row_of_coefficients(data, i, reference_frames, n_site, spins, spin_pairs, n_spin_pair, factor)
        for j=1:n_J
            Cs[i, j] = Ci[j]
        end
        Es[i] = 1000*(energies["conf2"] - energies["conf1"])
    end
    return (Cs, Es)
end


function get_Js(data::Dict, reference_frames::Array, n_site::Int64, spins::Array, spin_pairs::Array, n_spin_pair::Int64, n_energy_pair::Int64, factor::Float64, quantum_correction::Bool)
    Cs, Es = get_coefficient_matrix(data, reference_frames, n_site, spins, spin_pairs, n_spin_pair, n_energy_pair, factor)
    n_J = 9*n_spin_pair
    if n_energy_pair < n_J
        error(@sprintf("Underdetermined system: %d energy pair(s) for %d unknown(s).", n_energy_pair, n_J))
    end

    # Exact solve for a square system, least-squares fit for an overdetermined one
    Js = Cs \ Es
    if n_energy_pair > n_J
        @printf("Least-squares fit of %d energy pairs for %d unknowns: residual = %.3e meV\n", n_energy_pair, n_J, norm(Cs*Js - Es))
    end

    if quantum_correction
        # See Ruiz 2003
        # This is purely empirical, which doesn't make sense logically. 
        ind = 1
        for i_spin_pair=1:n_spin_pair
            i = spin_pairs[i_spin_pair][1]
            j = spin_pairs[i_spin_pair][2]
            for alpha=1:3
                for beta=1:3
                    Js[ind] *= 2*spins[i]*spins[j] / (2*spins[i]*spins[j] + min(spins[i], spins[j]))
                    ind += 1
                end
            end
        end
    end
    
    return Js
end


function print_Js(n_spin_pair::Int64, Js::Vector{Float64}, iformat)
    if iformat == 1
        for i=1:n_spin_pair
            @printf("Pair %3d\n", i)
            @printf("J_xx = %9.3f\n", Js[9*(i-1)+1])
            @printf("J_xy = %9.3f\n", Js[9*(i-1)+2])
            @printf("J_xz = %9.3f\n", Js[9*(i-1)+3])
            @printf("J_yx = %9.3f\n", Js[9*(i-1)+4])
            @printf("J_yy = %9.3f\n", Js[9*(i-1)+5])
            @printf("J_yz = %9.3f\n", Js[9*(i-1)+6])
            @printf("J_zx = %9.3f\n", Js[9*(i-1)+7])
            @printf("J_zy = %9.3f\n", Js[9*(i-1)+8])
            @printf("J_zz = %9.3f\n", Js[9*(i-1)+9])
        end
    elseif iformat == 2
        for i=1:n_spin_pair
            # @printf("Pair %3d\n", i)
            @printf(" %9.3f", Js[9*(i-1)+1])
            @printf(" %9.3f", Js[9*(i-1)+2])
            @printf(" %9.3f", Js[9*(i-1)+3])
            @printf(" %9.3f", Js[9*(i-1)+4])
            @printf(" %9.3f", Js[9*(i-1)+5])
            @printf(" %9.3f", Js[9*(i-1)+6])
            @printf(" %9.3f", Js[9*(i-1)+7])
            @printf(" %9.3f", Js[9*(i-1)+8])
            @printf(" %9.3f", Js[9*(i-1)+9])
        end
        println()
    end
end

function main(dirname::String)
    # Go to the working directory
    cd(dirname)

    # Comments on input file
    # The reference frame for each pair of spins must be specified.

    # Read input file
    data = YAML.load_file("input.yaml")
    n_site = length(data["spins"])
    n_spin_pair = length(data["exchange"])
    n_energy_pair = length(data["energies"])
    spins = data["spins"]

    # Gather spin pairs and reference frames
    spin_pairs = Any[]
    reference_frames = Any[]
    for i=1:n_spin_pair
        spin_pair = data["exchange"][i]["pair"]
        reference_frame = reshape(data["exchange"][i]["reference_frame"], 3, 3) # ! Transpose of emat
        push!(spin_pairs, spin_pair)
        push!(reference_frames, reference_frame)
    end

    # Solve exchange coupling constants
    factor = -2.0
    quantum_correction = false
    Js = get_Js(data, reference_frames, n_site, spins, spin_pairs, n_spin_pair, n_energy_pair, factor, quantum_correction)

    # Print out exchange coupling constants
    print_Js(n_spin_pair, Js, 2)
end

# Usage: julia main.jl [directory containing input.yaml]
if abspath(PROGRAM_FILE) == @__FILE__
    main(isempty(ARGS) ? pwd() : ARGS[1])
end

