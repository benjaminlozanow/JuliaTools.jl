# Tests for GC
# =====
#
#
# By Benjamin Lozano

using BioSequences
using Test

include("../GC.jl")

@testset "GC calculation tests" begin
    # Test 1: GC calculation for DNA sequence
    seq = dna"AGCT"
    gc = 0.5
    @test GC(seq) == gc

    # Test 2: GC calculation for RNA sequence
    seq = rna"AGCU"
    gc = 0.5
    @test GC(seq) == gc

    # Test 3: GC calculation with ambiguous characters
    seq = dna"AGCTN"
    gc = 0.5
    gc_ignore = 0.4
    gc_weighted = 0.5
    @test GC(seq, ambiguous="remove") == gc
    @test GC(seq, ambiguous="ignore") == gc_ignore
    @test GC(seq, ambiguous="weighted") == gc_weighted

    # Test 4: GC calculation for empty sequence
    seq = dna""
    gc = 0
    @test GC(seq) == gc

    # Test 5: Error handling for invalid ambiguous option
    seq = dna"AGCT"
    ambiguous = "invalid"
    error_msg = "ambiguous value invalid not recognized"
    @test_throws ArgumentError GC(seq, ambiguous=ambiguous)
    @test_throws ArgumentError(error_msg) GC(seq, ambiguous=ambiguous) 
end
