@testset "scaling_tree draws each generation about its scaling-law means with its own CV" begin
    # Deep enough that the last three generations hold 2048, 4096 and 8192 vessels. The CVs
    # are root first and differ there, the largest at the leaves, where σ = cv would put the
    # sample CV near 0.95 instead of 0.8. The exponents only need to be off Murray's and
    # nonzero; the sampler, not the anatomy, is under test.
    (depth, D_leaf, L_leaf, γ, β) = (13, 30.0, 600.0, 7 / 3, 1.2)
    cv_diameter = [fill(0.1, depth - 2); [0.3, 0.2, 0.25]]
    cv_length = [fill(0.3, depth - 2); [0.5, 0.3, 0.8]]
    build(seed) = scaling_tree(Xoshiro(seed), depth; leaf_diameter=D_leaf, leaf_length=L_leaf, γ, β, cv_diameter, cv_length)
    t = build(34)
    @test t.diameter == build(34).diameter && t.diameter != build(35).diameter
    @test all(i -> t.generation[i] == t.generation[t.parents[i]] + 1, 2:length(t.parents))
    @test count(==(depth), t.generation) == 2^depth

    # Standard errors of the sample mean and sample CV of n mean-one lognormal factors.
    average(x) = sum(x) / length(x)
    sample_cv(x) = sqrt(sum(abs2, x .- average(x)) / (length(x) - 1)) / average(x)
    excess_kurtosis(cv) = (s² = log1p(cv^2); exp(4s²) + 2exp(3s²) + 3exp(2s²) - 6)
    se_mean(cv, n) = cv / sqrt(n)
    se_cv(cv, n) = cv * sqrt((excess_kurtosis(cv) + 2) / (4n))
    println("  scaling_tree($depth; γ = $(round(γ; digits=3)), β = $β), target → sample:")
    for k in (depth-2):depth
        is = findall(==(k), t.generation)
        n = length(is)
        D̄ = D_leaf * 2.0^((depth - k) / γ)
        L̄ = L_leaf * (D̄ / D_leaf)^β
        (D, L) = (t.diameter[is], t.length[is])
        (cv_D, cv_L) = (cv_diameter[k+1], cv_length[k+1])
        (logD, logL) = (log.(D) .- average(log.(D)), log.(L) .- average(log.(L)))
        correlation = sum(logD .* logL) / sqrt(sum(abs2, logD) * sum(abs2, logL))
        println(
            "    level $k (n = $n): D̄ $(round(D̄; digits=2)) → $(round(average(D); digits=2)) µm, ",
            "cv_D $cv_D → $(round(sample_cv(D); digits=3)), L̄ $(round(L̄; digits=1)) → $(round(average(L); digits=1)) µm, ",
            "cv_L $cv_L → $(round(sample_cv(L); digits=3)), corr(log D, log L) $(round(correlation; digits=3))",
        )
        @test abs(average(D) / D̄ - 1) ≤ 4se_mean(cv_D, n)
        @test abs(average(L) / L̄ - 1) ≤ 4se_mean(cv_L, n)
        @test abs(sample_cv(D) - cv_D) ≤ 4se_cv(cv_D, n)
        @test abs(sample_cv(L) - cv_L) ≤ 4se_cv(cv_L, n)
        @test abs(correlation) ≤ 4 / sqrt(n)        # β scales the generation's mean, not the drawn diameter
    end

    # A Float32 tree must stay Float32 through the exponents and the draws.
    for T in TEST_TYPES
        s = scaling_tree(Xoshiro(1), 3; leaf_diameter=T(30), leaf_length=T(600), β=1.2, cv_diameter=0.1, cv_length=0.3)
        @test s isa ScalingTree{T}
    end
end
