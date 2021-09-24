# Schemes for doing a type of statistics on the system
abstract type StatsScheme <: Scheme end

function run_scheme(wf, ham, metro, dimsnum, nthreads, scheme::StatsScheme)
    result = vmc(wf, ham, metro, dimsnum, nthreads, scheme)
    return result
end

struct OneBody <: StatsScheme
    start::Float64
    stop::Float64
    length::Float64
end
OneBody(dims, num; start, stop, length, dimsnum) = OneBody(start, stop, length)
createSampler(scheme::OneBody, wf, sampled_steps) = OneBodySampler(scheme.start, scheme.stop, scheme.length, sampled_steps)

struct Blocking <: StatsScheme end
createSampler(scheme::Blocking, wf, sampled_steps) = BlockingSampler(sampled_steps)