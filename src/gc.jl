# GC calculator
# =====
#
# Calculate the fraction of nucleotide which are G or C
#
# By Benjamin Lozano

using BioSequences

function gc(seq::LongNuc; ambiguous::String="remove")

    GCs = [DNA_C, RNA_C, DNA_G, RNA_G, DNA_S, RNA_S]

    ATs = [DNA_A, RNA_A, DNA_T, RNA_U, DNA_W, RNA_W] 

    gc = sum(count(x -> x in GCs, seq))

    if !(ambiguous in ["weighted", "remove", "ignore"])
        throw(ArgumentError("ambiguous value $ambiguous not recognized"))
    end

    if ambiguous == "remove"
        len = gc + sum(count(x -> x in ATs, seq))
    else
        len = length(seq)
    end

    if ambiguous == "weighted"
        for nt in seq
            if !(nt in ATs) && !(nt in GCs)
                gc += gc_values[nt]
            end
        end
    end

    if len == 0
        return 0
    end

    return gc / len
end

gc_values = Dict{NucleicAcid, Float64}(

    #DNAs
    DNA_G => 1.000,
    DNA_C => 1.000,
    DNA_A => 0.000,
    DNA_T => 0.000,
    DNA_S => 1.000,  # Strong interaction (3 H bonds) (G or C)
    DNA_W => 0.000,  # Weak interaction (2 H bonds) (A or T)
    DNA_M => 0.500,  # Amino (A or C)
    DNA_R => 0.500,  # Purine (A or G)
    DNA_Y => 0.500,  # Pyrimidine (T or C)
    DNA_K => 0.500,  # Keto (G or T)
    DNA_V => 2 / 3,  # Not T or U (A or C or G)
    DNA_B => 2 / 3,  # Not A (C or G or T)
    DNA_H => 1 / 3,  # Not G (A or C or T)
    DNA_D => 1 / 3,  # Not C (A or G or T)
    DNA_N => 0.500,  # Any nucleotide (A or C or G or T)

    # # RNAs
    RNA_G => 1.000,
    RNA_C => 1.000,
    RNA_A => 0.000,
    RNA_U => 0.000,
    RNA_S => 1.000,  # Strong interaction (3 H bonds) (G or C)
    RNA_W => 0.000,  # Weak interaction (2 H bonds) (A or U)
    RNA_M => 0.500,  # Amino (A or C)
    RNA_R => 0.500,  # Purine (A or G)
    RNA_Y => 0.500,  # Pyrimidine (U or C)
    RNA_K => 0.500,  # Keto (G or U)
    RNA_V => 2 / 3,  # Not U or U (A or C or G)
    RNA_B => 2 / 3,  # Not A (C or G or U)
    RNA_H => 1 / 3,  # Not G (A or C or U)
    RNA_D => 1 / 3,  # Not C (A or G or U)
    RNA_N => 0.500  # Any nucleotide (A or C or G or U)
)
