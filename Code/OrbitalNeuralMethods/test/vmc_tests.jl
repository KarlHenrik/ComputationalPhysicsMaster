@testset "Correlated Optimization" begin
    dims = 3
    num = 5

    ham = HarmonicOscillator(dims, num, ω = 1)
    metro = Importance(equils = 1e4, samples = 1e5, step = 0.01)
    optimizer = GradientDescent(lr = 0.01, max_iter = 20, tol = 1e-6)

    for (a, α_opt, E) in zip([0, 0.00433, 0.0433], [0.5, 0.498944, 0.490226], [1.5, 1.50684, 1.56743])
        wf = Correlated(dims, num, α = 0.5, a = a)
        wf_opt, grad_results = optimize(wf, ham, metro, optimizer, nthreads = 8, verbose = false)

        @test abs(wf_opt.α - α_opt) < 1e-2
        @test (grad_results[end].E/num - E) < 1e-4
    end
end

@testset "SimpleGaussian Optimization" begin
    imp = Importance(equils = 1e3, samples = 1e4, step = 0.01)
    met = Metropolis(equils = 1e3, samples = 1e4, step = 0.1)

    optimizer = GradientDescent(lr = 0.1, max_iter = 20, tol = 1e-6)
    for metro in [met, imp]
        for dims in 1:3
            for num in [1, 5, 10]
                ham = HarmonicOscillator(dims, num, ω = 1)
                wf = SimpleGaussian(dims, num, α = 0.48)

                wf_opt, grad_results = optimize(wf, ham, metro, optimizer, nthreads = 8, verbose = false)

                @test abs(wf_opt.α - 0.5) < 1e-3
                @test (grad_results[end].E - num * dims) < 1e-4
            end
        end
    end
end